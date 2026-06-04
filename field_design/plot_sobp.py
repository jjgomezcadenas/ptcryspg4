#!/usr/bin/env python3
"""Plot the actual Geant4 SOBP depth-dose from data/depth_dose.csv.

depth_dose.csv (written by proton_transport every run) is the *realized* field:
energy deposit per 1 mm depth bin, total and primary-proton-only. This plots it
versus depth from the entrance face, marks the target, and reports the plateau
flatness — the authoritative check of the SOBP design.

Usage:
    python field_design/plot_sobp.py [data_dir] [--prox CM] [--dist CM]
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?",
                    default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--prox", type=float, default=5.5, help="target proximal depth [cm]")
    ap.add_argument("--dist", type=float, default=10.5, help="target distal depth [cm]")
    args = ap.parse_args()

    d = pd.read_csv(os.path.join(args.data_dir, "depth_dose.csv"))
    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    half_len_mm = 0.5 * float(meta["phantom_length_mm"])

    # Depth from the entrance face (beam enters at z = -half_len).
    depth = (d["z_mm"] + half_len_mm) / 10.0  # cm
    tot = d["edep_total_MeV"].to_numpy()
    prim = d["edep_primary_MeV"].to_numpy()

    # Plateau flatness over the target, excluding the last 0.5 cm distal falloff.
    flat = (depth >= args.prox) & (depth <= args.dist - 0.5)
    p = tot[flat]
    flatness = (p.max() - p.min()) / p.mean() * 100 if len(p) else float("nan")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(depth, tot, "k", label="total")
    ax.plot(depth, prim, "C0--", label="primary proton")
    ax.axvspan(args.prox, args.dist, color="C2", alpha=0.15, label="target")
    ax.set(xlabel="depth from entrance [cm]", ylabel="energy deposit [MeV / mm]",
           title=f"G4 SOBP depth-dose  ({meta['phantom_material']}, "
                 f"{int(meta['n_protons'])} p)  -  plateau flatness {flatness:.1f}%")
    ax.legend()
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    out = os.path.join(args.data_dir, "sobp_g4.png")
    fig.savefig(out, dpi=120)
    print(f"plateau flatness over [{args.prox}, {args.dist - 0.5}] cm: {flatness:.1f}%")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
