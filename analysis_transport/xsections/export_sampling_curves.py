"""Export the nominal fitted curves as one compact table for the Geant4 sampler.

`data/xsections/fits/sampling_curves.csv` holds the five channels in long
form (channel, target, residual, threshold, energy grid, cross section,
envelope), so the sampling application parses a single flat file whose
provenance is the fit directory.  The envelope column is the bank sampling
curve: the nominal scaled by the channel's replica cover factor (the largest
replica-to-nominal ratio over the production-bearing energies, plus a 10%
margin), so the bank keep-probability covers every replica.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .channels import CHANNELS

REPLICA_COLUMN = re.compile(r"^sigma_([0-9]+(?:\.[0-9]+)?)_MeV_mb$")


def cover_factor(fits: Path, channel_id: str) -> float:
    """Largest replica/nominal ratio over production-bearing energies, x1.1."""
    replicas = pd.read_csv(fits / f"{channel_id}_replicas.csv")
    table = pd.read_csv(fits / f"{channel_id}_table.csv")
    columns = [c for c in replicas.columns if REPLICA_COLUMN.match(c)]
    energies = np.asarray([float(REPLICA_COLUMN.match(c).group(1))
                           for c in columns])
    order = np.argsort(energies)
    replica_max = replicas[columns].to_numpy(dtype=float)[:, order].max(axis=0)
    nominal = np.interp(energies[order], table.energy_MeV,
                        table.sigma_nominal_mb)
    bearing = nominal >= 0.01 * nominal.max()
    return 1.1 * float(np.max(replica_max[bearing] / nominal[bearing]))


def export(repo: Path) -> Path:
    fits = repo / "data/xsections/fits"
    thresholds = json.loads((fits / "fit_meta.json").read_text())["threshold_MeV"]
    output = fits / "sampling_curves.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["channel_id", "target", "residual", "threshold_MeV",
                         "energy_MeV", "sigma_mb", "sigma_env_mb"])
        for channel in CHANNELS:
            curve = pd.read_csv(fits / f"{channel.channel_id}_curve.csv")
            threshold = thresholds[channel.channel_id]
            factor = cover_factor(fits, channel.channel_id)
            for row in curve.itertuples(index=False):
                writer.writerow([
                    channel.channel_id, channel.target, channel.residual,
                    threshold, row.energy_MeV, row.sigma_nominal_mb,
                    factor * row.sigma_nominal_mb,
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
