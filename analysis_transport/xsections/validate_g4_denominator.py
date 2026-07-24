"""Validate the direct Geant4 denominator pilot and support audit."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from .channels import CHANNELS
from .config import load


def validate(repo: Path, config_path: Path, executable: Path):
    config = load(config_path)
    root = (repo / "data/xsections/g4"
            / f"geant4-{config['geant4_version']}" / config["physics_list"])
    counts = pd.read_csv(root / "denominator_pilot_counts.csv")
    metadata = json.loads(
        (root / "denominator_pilot_meta.json").read_text(encoding="utf-8"))
    problems = []
    expected_rows = len(config["targets"]) * len(config["energies_MeV"])
    if len(counts) != expected_rows:
        problems.append("incomplete target-energy pilot grid")
    if counts.groupby("target").size().to_dict() != {
            target: len(config["energies_MeV"]) for target in config["targets"]}:
        problems.append("unequal target coverage")
    if not (counts.n_interactions
            == int(config["denominator_pilot_interactions"])).all():
        problems.append("wrong pilot interaction count")
    if not (counts.sigma_inelastic_mb > 0).all():
        problems.append("nonpositive inelastic cross section")
    if set(counts.physics_list) != {config["physics_list"]}:
        problems.append("physics-list mismatch")
    if not counts.cross_section_data_sets.str.len().gt(0).all():
        problems.append("missing cross-section data-set provenance")
    if not counts.model_counts.str.len().gt(0).all():
        problems.append("missing final-state model provenance")
    if metadata.get("config_sha256") != hashlib.sha256(
            config_path.read_bytes()).hexdigest():
        problems.append("denominator configuration hash mismatch")
    if metadata.get("executable_sha256") != hashlib.sha256(
            executable.read_bytes()).hexdigest():
        problems.append("denominator executable hash mismatch")

    pilot = root / "denominator_pilot"
    summary = pd.read_csv(pilot / "support_summary.csv")
    if set(summary.channel_id) != {channel.channel_id for channel in CHANNELS}:
        problems.append("incomplete support summary")
    for channel in CHANNELS:
        curve = pd.read_csv(pilot / f"{channel.channel_id}_curve.csv")
        if len(curve) != len(config["energies_MeV"]):
            problems.append(f"incomplete denominator curve for {channel.channel_id}")
        if (curve.sigma_mb < 0).any():
            problems.append(f"negative denominator for {channel.channel_id}")
    c12 = summary[summary.channel_id == "p_C12_x_C11"].iloc[0]
    if int(c12.zero_g4_support_in_production_region) != 8:
        problems.append("C-12 support result changed")
    if str(c12.zero_support_energies_MeV) != "20;22;24;26;28;30;35;40":
        problems.append("C-12 zero-support energies changed")
    figure = repo / "docs/figures/xsection_denominator/pilot_support.pdf"
    if not figure.exists() or figure.stat().st_size == 0:
        problems.append("missing denominator support figure")
    if problems:
        raise ValueError("; ".join(problems))
    return expected_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--config", type=Path,
                        default=Path("config/xsections_denominator.toml"))
    parser.add_argument("--executable", type=Path,
                        default=Path("xsections_g4/build/denominator_sampler"))
    args = parser.parse_args()
    rows = validate(args.repo.resolve(), args.config.resolve(),
                    args.executable.resolve())
    print(f"Validated {rows} direct Geant4 denominator pilot rows")


if __name__ == "__main__":
    main()
