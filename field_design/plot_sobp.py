#!/usr/bin/env python3
"""Plot the realized Geant4 SOBP depth-dose of a run.

depth_dose.csv (written by proton_transport every run) is the *realized* field:
energy deposit per depth bin, total and primary-proton-only. This plots it
versus depth from the entrance face, marks the target window read from the
run's run_meta.csv, and quotes the plateau uniformity graded by sobp_metrics
(on the central-axis core dose — the acceptance number check_run gates on).

Usage:
    python field_design/plot_sobp.py <run_dir> [--prox CM] [--dist CM]
"""

import argparse
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "analysis_transport"))
from sobp_metrics import metrics as plateau_metrics  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", help="a Stage-A run directory")
    ap.add_argument("--prox", type=float, default=None,
                    help="target proximal depth [cm] (default: the run's)")
    ap.add_argument("--dist", type=float, default=None,
                    help="target distal depth [cm] (default: the run's)")
    args = ap.parse_args()

    d = pd.read_csv(os.path.join(args.data_dir, "depth_dose.csv"))
    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    half_len_mm = 0.5 * float(meta["phantom_length_mm"])
    prox = args.prox if args.prox is not None else float(meta["target_prox_depth_mm"]) / 10.0
    dist = args.dist if args.dist is not None else float(meta["target_dist_depth_mm"]) / 10.0

    # Plateau uniformity from sobp_metrics (central-axis core dose), for an
    # SOBP run whose depth_dose carries the core tally.
    is_sobp = os.path.exists(os.path.join(args.data_dir, "sobp_layers.csv"))
    if is_sobp:
        try:
            grade = (f"plateau uniformity "
                     f"{plateau_metrics(args.data_dir)['uniformity_pct']:.1f}% (core)")
        except ValueError:
            grade = "plateau uniformity n/a (no central-axis tally)"
    else:
        grade = "pencil run"

    # Depth from the entrance face (beam enters at z = -half_len).
    depth = (d["z_mm"] + half_len_mm) / 10.0  # cm
    tot = d["edep_total_MeV"].to_numpy()
    prim = d["edep_primary_MeV"].to_numpy()

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(depth, tot, "k", label="total")
    ax.plot(depth, prim, "C0--", label="primary proton")
    ax.axvspan(prox, dist, color="C2", alpha=0.15, label="target")
    ax.set(xlabel="depth from entrance [cm]", ylabel="energy deposit [MeV / mm]",
           title=f"G4 depth-dose  ({meta['phantom_material']}, "
                 f"{int(meta['n_protons'])} p)  -  {grade}")
    ax.legend()
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    figdir = os.path.join(args.data_dir, "figures")
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "sobp_g4.png")
    fig.savefig(out, dpi=120)
    print(f"target window [{prox:g}, {dist:g}] cm; {grade}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
