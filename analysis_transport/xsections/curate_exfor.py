"""Create the tracked dataset-level EXFOR curation audit."""

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from .channels import BY_PAIR
from .config import load


VALID_STATES = {"accepted", "excluded", "pending"}


def reaction_suffix(reaction: str):
    match = re.search(r"\)\d+-[A-Z]+-\d+(.*)$", reaction, re.IGNORECASE)
    return match.group(1) if match else ""


def _point_counts(repo, catalog_row, threshold, energy_min, energy_max):
    points = pd.read_csv(repo / catalog_row["point_file"])
    uncertainty = points[
        ["sigma_unc_minus_mb", "sigma_unc_plus_mb"]].mean(axis=1)
    in_range = points.energy_MeV.between(energy_min, energy_max)
    positive_uncertainty = uncertainty.notna() & (uncertainty > 0)
    point_fit_rule = (in_range & (points.energy_MeV > threshold)
                      & (points.sigma_mb > 0) & positive_uncertainty)
    return {
        "points_total": len(points),
        "points_in_study_range": int(in_range.sum()),
        "points_with_positive_uncertainty": int(positive_uncertainty.sum()),
        "points_passing_point_rules": int(point_fit_rule.sum()),
    }


def _point_audit(repo, catalog_row, dataset_state, threshold,
                 energy_min, energy_max):
    points = pd.read_csv(repo / catalog_row["point_file"])
    uncertainty = points[
        ["sigma_unc_minus_mb", "sigma_unc_plus_mb"]].mean(axis=1)
    rows = []
    for index, point in points.iterrows():
        if dataset_state != "accepted":
            include = False
            reason = f"dataset_{dataset_state}"
        elif not energy_min <= point.energy_MeV <= energy_max:
            include = False
            reason = "outside_energy_range"
        elif point.sigma_mb <= 0:
            include = False
            reason = "nonpositive_cross_section"
        elif point.energy_MeV <= threshold:
            include = False
            reason = "at_or_below_threshold"
        elif pd.isna(uncertainty.iloc[index]) or uncertainty.iloc[index] <= 0:
            include = False
            reason = "no_positive_reported_uncertainty"
        else:
            include = True
            reason = "accepted"
        rows.append({
            "dataset_id": catalog_row["dataset_id"],
            "point_id": int(point.point_id),
            "energy_MeV": point.energy_MeV,
            "sigma_mb": point.sigma_mb,
            "sigma_unc_mb": uncertainty.iloc[index],
            "include_in_fit": int(include),
            "reason_code": reason,
        })
    return rows


def curate(repo: Path, curation_path: Path, fit_path: Path):
    policy = load(curation_path)
    fit_config = load(fit_path)
    catalog_path = repo / "data/xsections/normalized/datasets.csv"
    catalog = pd.read_csv(catalog_path)
    exfor_catalog = catalog[catalog.library == "EXFOR"]
    overrides = policy.get("dataset_override", {})
    rows = []
    point_rows = []
    for catalog_row in exfor_catalog.to_dict("records"):
        metadata = json.loads(
            (repo / catalog_row["metadata_file"]).read_text(encoding="utf-8"))
        channel = BY_PAIR[(metadata["target"], metadata["residual"])]
        suffix = reaction_suffix(metadata["reaction"])
        checks = {
            "channel_identity": (metadata["target"], metadata["residual"]) in BY_PAIR,
            "absolute_quantity": suffix in policy["accepted_reaction_suffixes"],
            "energy_unit": metadata["original_energy_unit"] == policy["accepted_energy_unit"],
            "cross_section_unit": (
                metadata["original_cross_section_unit"]
                == policy["accepted_cross_section_unit"]),
            "quantity_label": "cross section" in metadata["quantity"].lower(),
        }
        if metadata["dataset_id"] in overrides:
            decision = overrides[metadata["dataset_id"]]
            state = decision["state"]
            reason_code = decision["reason_code"]
            verification_basis = decision["verification_basis"]
            notes = decision.get("notes", "")
        elif all(checks.values()):
            state = "accepted"
            reason_code = "absolute_exfor_cross_section"
            verification_basis = "EXFOR reaction, quantity and source units"
            notes = ""
        else:
            state = "pending"
            reason_code = "curation_check_failed"
            verification_basis = ";".join(
                name for name, passed in checks.items() if not passed)
            notes = "Requires source-level verification before fitting."
        if state not in VALID_STATES:
            raise ValueError(
                f"Invalid curation state {state} for {metadata['dataset_id']}")
        counts = _point_counts(
            repo, catalog_row,
            fit_config["threshold_MeV"][channel.channel_id],
            fit_config["energy_min_MeV"], fit_config["energy_max_MeV"])
        point_rows.extend(_point_audit(
            repo, catalog_row, state,
            fit_config["threshold_MeV"][channel.channel_id],
            fit_config["energy_min_MeV"], fit_config["energy_max_MeV"]))
        rows.append({
            "dataset_id": metadata["dataset_id"],
            "channel_id": channel.channel_id,
            "accession": metadata["accession"],
            "author": metadata["author"],
            "reference": metadata["reference"],
            "doi": metadata["doi"],
            "reaction": metadata["reaction"],
            "quantity": metadata["quantity"],
            "original_energy_unit": metadata["original_energy_unit"],
            "original_cross_section_unit": metadata["original_cross_section_unit"],
            "state": state,
            "reason_code": reason_code,
            "verification_basis": verification_basis,
            "notes": notes,
            **counts,
        })

    output = repo / "data/xsections/curation.csv"
    frame = pd.DataFrame(rows).sort_values(
        ["channel_id", "accession"]).reset_index(drop=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False, lineterminator="\n")
    point_frame = pd.DataFrame(point_rows).sort_values(
        ["dataset_id", "point_id"]).reset_index(drop=True)
    point_frame.to_csv(repo / "data/xsections/point_curation.csv",
                       index=False, lineterminator="\n")
    metadata = {
        "schema_version": 1,
        "dataset_count": len(frame),
        "point_count": len(point_frame),
        "fit_point_count": int(point_frame.include_in_fit.sum()),
        "state_counts": frame.state.value_counts().sort_index().to_dict(),
        "curation_config": str(curation_path.relative_to(repo)),
        "curation_config_sha256": hashlib.sha256(
            curation_path.read_bytes()).hexdigest(),
        "fit_config_sha256": hashlib.sha256(fit_path.read_bytes()).hexdigest(),
        "normalized_catalog_sha256": hashlib.sha256(
            catalog_path.read_bytes()).hexdigest(),
        "modifier_reference": (
            "IAEA EXFOR Formats/LEXFOR, reaction SF8 general quantity modifier A"
        ),
    }
    (repo / "data/xsections/curation_meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--curation-config", type=Path,
                        default=Path("config/xsection_curation.toml"))
    parser.add_argument("--fit-config", type=Path,
                        default=Path("config/xsection_fit.toml"))
    args = parser.parse_args()
    frame = curate(args.repo.resolve(), args.curation_config.resolve(),
                   args.fit_config.resolve())
    counts = frame.state.value_counts().sort_index().to_dict()
    print(f"Curated {len(frame)} EXFOR series: {counts}")


if __name__ == "__main__":
    main()
