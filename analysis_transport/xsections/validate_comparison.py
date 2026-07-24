"""Validate normalized inputs and generated cross-section comparison products."""

import argparse
import json
import hashlib
from pathlib import Path

import pandas as pd

from .channels import CHANNELS


def validate(repo: Path):
    problems = []
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    if catalog.dataset_id.duplicated().any(): problems.append("duplicate dataset identifiers")
    for row in catalog.to_dict("records"):
        points = pd.read_csv(repo / row["point_file"])
        metadata_path = repo / row["metadata_file"]
        if not metadata_path.exists(): problems.append(f"missing metadata for {row['dataset_id']}")
        else:
            source_file = json.loads(metadata_path.read_text()).get("source_file", "")
            if not source_file or not (repo / source_file).exists():
                problems.append(f"missing raw source for {row['dataset_id']}")
        if points.point_id.duplicated().any(): problems.append(f"duplicate point ids in {row['dataset_id']}")
        if not points.energy_MeV.between(0, 1.0e6).all(): problems.append(f"invalid energy in {row['dataset_id']}")
        if not points.sigma_mb.ge(0).all(): problems.append(f"negative cross section in {row['dataset_id']}")
    for channel in CHANNELS:
        figure = repo / "docs/figures/xsections" / f"{channel.channel_id}.pdf"
        if not figure.exists() or figure.stat().st_size == 0: problems.append(f"missing figure {figure}")
    metadata = json.loads((repo / "docs/generated/xsections/comparison_meta.json").read_text())
    if metadata["normalized_dataset_count"] != len(catalog): problems.append("dataset count mismatch")
    sources = pd.read_csv(repo / "data/xsections/sources.csv").fillna("")
    for source in sources.to_dict("records"):
        if not source["sha256"]:
            continue
        local_path = repo / "data/xsections" / source["local_path"]
        if hashlib.sha256(local_path.read_bytes()).hexdigest() != source["sha256"]:
            problems.append(f"checksum mismatch for {source['source_id']}")
    g4_counts = repo / "data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP/thin_target_counts.csv"
    if g4_counts.exists():
        counts = pd.read_csv(g4_counts)
        if (counts.n_protons <= 0).any(): problems.append("invalid Geant4 proton count")
        if (counts.mean_continuous_loss_MeV >= 0.5).any(): problems.append("Geant4 energy-loss criterion failed")
        if ((counts.n_inelastic / counts.n_protons) >= 1e-3).any(): problems.append("Geant4 interaction criterion failed")
        run_meta = json.loads((g4_counts.parent / "run_meta.json").read_text())
        config_path = repo / "config/xsections_scan.toml"
        current_hash = hashlib.sha256(config_path.read_bytes()).hexdigest()
        if run_meta.get("config_sha256") != current_hash: problems.append("Geant4 scan configuration hash mismatch")
        if set(counts.target) != {"C12", "N14", "O16"}: problems.append("incomplete Geant4 target coverage")
        if counts.groupby("target").size().to_dict() != {"C12": 36, "N14": 36, "O16": 36}:
            problems.append("incomplete Geant4 energy coverage")
    if problems:
        raise ValueError("; ".join(problems))
    return len(catalog)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    count = validate(args.repo.resolve())
    print(f"Validated {count} normalized datasets and five comparison figures")


if __name__ == "__main__":
    main()
