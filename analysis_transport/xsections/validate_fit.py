"""Validate the EXFOR-fit products used by the standalone report."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .channels import CHANNELS
from .config import load


def validate(repo: Path, config_path: Path):
    config = load(config_path)
    problems = []
    fit_dir = repo / "data/xsections/fits"
    summary_path = fit_dir / "fit_summary.csv"
    if not summary_path.exists():
        raise ValueError("missing fit summary")
    summary = pd.read_csv(summary_path)
    curation_path = repo / "data/xsections/curation.csv"
    point_curation_path = repo / "data/xsections/point_curation.csv"
    curation = pd.read_csv(curation_path)
    point_curation = pd.read_csv(point_curation_path)
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    exfor_ids = set(catalog.loc[catalog.library == "EXFOR", "dataset_id"])
    if set(curation.dataset_id) != exfor_ids:
        problems.append("curation does not cover every EXFOR dataset exactly once")
    if not set(curation.state).issubset({"accepted", "excluded", "pending"}):
        problems.append("invalid dataset curation state")
    if curation.dataset_id.duplicated().any():
        problems.append("duplicate dataset curation decision")
    accepted = curation[curation.state == "accepted"]
    if accepted.reaction.str.contains(",A,", regex=False).any():
        problems.append("abundance-weighted EXFOR quantity accepted")
    b0077 = curation[curation.accession.isin(["B0077.002", "B0077.003"])]
    if len(b0077) != 2 or not (b0077.state == "pending").all():
        problems.append("B0077 curation decision changed")
    b0095 = curation[curation.accession == "B0095.002"]
    if len(b0095) != 1 or b0095.state.iloc[0] != "pending":
        problems.append("B0095 must remain pending")
    masuda = curation[curation.accession.isin(
        ["E2568.002", "E2568.003", "E2568.004"])]
    if len(masuda) != 3 or not (masuda.state == "pending").all():
        problems.append("Masuda series must remain pending shape comparisons")
    if set(summary.channel_id) != {channel.channel_id for channel in CHANNELS}:
        problems.append("fit summary has incomplete channel coverage")
    for channel in CHANNELS:
        stem = fit_dir / channel.channel_id
        curve = pd.read_csv(Path(f"{stem}_curve.csv"))
        replicas = pd.read_csv(Path(f"{stem}_replicas.csv"))
        representatives = pd.read_csv(Path(f"{stem}_representatives.csv"))
        if len(curve) != int(config["fit_grid_points"]):
            problems.append(f"wrong dense-grid size for {channel.channel_id}")
        numeric = curve.select_dtypes(include=[np.number]).to_numpy()
        if not np.isfinite(numeric).all():
            problems.append(f"non-finite curve value for {channel.channel_id}")
        threshold = config["threshold_MeV"][channel.channel_id]
        below = curve.energy_MeV <= threshold
        if not (curve.loc[below, ["sigma_nominal_mb", "sigma_lower_16_mb",
                                  "sigma_upper_84_mb"]] == 0).all().all():
            problems.append(f"nonzero subthreshold fit for {channel.channel_id}")
        if (curve.sigma_nominal_mb < 0).any():
            problems.append(f"negative fit for {channel.channel_id}")
        if (curve.sigma_lower_16_mb > curve.sigma_nominal_mb).any():
            problems.append(f"lower quantile exceeds median for {channel.channel_id}")
        if (curve.sigma_upper_84_mb < curve.sigma_nominal_mb).any():
            problems.append(f"upper quantile below median for {channel.channel_id}")
        if len(replicas) != int(config["n_replicas"]):
            problems.append(f"wrong replica count for {channel.channel_id}")
        expected_representatives = len(config["representative_distance_quantiles"])
        if representatives.representative_rank.nunique() != expected_representatives:
            problems.append(f"wrong representative count for {channel.channel_id}")
        figure = repo / "docs/figures/xsection_fit" / f"{channel.channel_id}.pdf"
        if not figure.exists() or figure.stat().st_size == 0:
            problems.append(f"missing figure for {channel.channel_id}")

    sensitivity = pd.read_csv(fit_dir / "sensitivity_summary.csv")
    expected_sensitivity = set(config.get("sensitivity_include_dataset_ids", []))
    actual_sensitivity = set(
        ";".join(sensitivity.included_pending_dataset_ids).split(";"))
    if actual_sensitivity != expected_sensitivity:
        problems.append("sensitivity dataset coverage mismatch")
    c12_sensitivity = fit_dir / "p_C12_x_C11_sensitivity.csv"
    if not c12_sensitivity.exists() or c12_sensitivity.stat().st_size == 0:
        problems.append("missing B0095 sensitivity curve")

    metadata = json.loads((fit_dir / "fit_meta.json").read_text())
    expected_hash = hashlib.sha256(config_path.read_bytes()).hexdigest()
    if metadata.get("config_sha256") != expected_hash:
        problems.append("fit configuration hash mismatch")
    if metadata.get("curation_sha256") != hashlib.sha256(
            curation_path.read_bytes()).hexdigest():
        problems.append("fit curation hash mismatch")
    for channel in CHANNELS:
        expected_points = int(point_curation[
            point_curation.dataset_id.isin(
                curation.loc[curation.channel_id == channel.channel_id,
                             "dataset_id"])
            & (point_curation.include_in_fit == 1)
        ].shape[0])
        actual_points = int(summary.loc[
            summary.channel_id == channel.channel_id, "points_used"].iloc[0])
        if expected_points != actual_points:
            problems.append(f"point-audit mismatch for {channel.channel_id}")
    comparison = pd.read_csv(
        repo / "docs/generated/xsection_fit/evaluation_comparison.csv")
    expected_pairs = {(channel.channel_id, evaluation)
                      for channel in CHANNELS
                      for evaluation in ("JENDL-4.0/HE", "LANL ENDF/B-VII.1")}
    actual_pairs = set(zip(comparison.channel_id, comparison.evaluation))
    if actual_pairs != expected_pairs:
        problems.append("incomplete JENDL/LANL comparison coverage")
    if problems:
        raise ValueError("; ".join(problems))
    return len(CHANNELS), int(config["n_replicas"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--config", type=Path,
                        default=Path("config/xsection_fit.toml"))
    args = parser.parse_args()
    channels, replicas = validate(args.repo.resolve(), args.config.resolve())
    print(f"Validated {replicas} replicas and one figure for each of {channels} channels")


if __name__ == "__main__":
    main()
