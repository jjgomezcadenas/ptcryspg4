"""Export the nominal fitted curves as one compact table for the Geant4 sampler.

`data/xsections/fits/sampling_curves.csv` holds the five channels in long
form (channel, target, residual, threshold, energy grid, cross section), so
the sampling application parses a single flat file whose provenance is the
fit directory.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pandas as pd

from .channels import CHANNELS


def export(repo: Path) -> Path:
    fits = repo / "data/xsections/fits"
    thresholds = json.loads((fits / "fit_meta.json").read_text())["threshold_MeV"]
    output = fits / "sampling_curves.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["channel_id", "target", "residual", "threshold_MeV",
                         "energy_MeV", "sigma_mb"])
        for channel in CHANNELS:
            curve = pd.read_csv(fits / f"{channel.channel_id}_curve.csv")
            threshold = thresholds[channel.channel_id]
            for row in curve.itertuples(index=False):
                writer.writerow([
                    channel.channel_id, channel.target, channel.residual,
                    threshold, row.energy_MeV, row.sigma_nominal_mb,
                ])
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    output = export(args.repo.resolve())
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
