#!/usr/bin/env python3
"""Render the scenario coordinate frame: phantom, target box, and the source.

A z-y cross-section (the beam axis is z) drawn from the actual run_meta.csv
geometry, with a subsample of the emitters.csv annihilation points overlaid so
the source cloud is shown sitting inside the medium. This is the figure a
downstream consumer needs to co-register the source with the phantom: the
phantom is centred at the origin, its axis is +z (the beam direction), and the
beam enters at the z = -L/2 face. Used in docs/scenario_format.tex.

Usage:
    python analysis_transport/plot_geometry.py [data_dir] [--out PNG] [--n N]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?", default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--out", default=None, help="output PNG (default data/scenario_geometry.png)")
    ap.add_argument("--n", type=int, default=15000, help="annihilation points to scatter")
    args = ap.parse_args()

    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    emit = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))

    half_z = 0.5 * float(meta["phantom_length_mm"])
    R = 0.5 * float(meta["phantom_diameter_mm"])
    tr = float(meta["target_radius_mm"])
    z_lo = float(meta["target_prox_depth_mm"]) - half_z
    z_hi = float(meta["target_dist_depth_mm"]) - half_z

    # Subsample the annihilation points for a legible scatter (reproducible).
    if len(emit) > args.n:
        emit = emit.sample(n=args.n, random_state=0)

    fig, ax = plt.subplots(figsize=(8, 5))

    # Source cloud, coloured by isotope (z-y projection of the annihilation points).
    for iid in sorted(ISOTOPES):
        sel = emit["isotope_id"] == iid
        if sel.any():
            ax.scatter(emit.loc[sel, "anh_z_mm"], emit.loc[sel, "anh_y_mm"],
                       s=2, alpha=0.25, linewidths=0, label=ISOTOPES[iid].name)

    # Phantom cross-section (cylinder axis along z) and target box.
    ax.add_patch(Rectangle((-half_z, -R), 2 * half_z, 2 * R, fill=False,
                           edgecolor="k", lw=1.8, label="phantom"))
    ax.add_patch(Rectangle((z_lo, -tr), z_hi - z_lo, 2 * tr, fill=False,
                           edgecolor="C3", lw=1.6, ls="--", label="target box"))

    # Origin, entrance face, and the beam arrow (+z).
    ax.plot(0, 0, "k+", ms=12, mew=2)
    ax.annotate("origin", (0, 0), textcoords="offset points", xytext=(6, 6))
    ax.axvline(-half_z, color="0.5", lw=0.8, ls=":")
    ax.annotate("entrance\n(z = -L/2)", (-half_z, R), textcoords="offset points",
                xytext=(4, -14), color="0.4", fontsize=9)
    ax.annotate("", xy=(-half_z, 0), xytext=(-half_z - 18, 0),
                arrowprops=dict(arrowstyle="-|>", color="C0", lw=2))
    ax.annotate("beam +z", (-half_z - 18, 0), textcoords="offset points",
                xytext=(-2, 6), color="C0")

    ax.set(xlabel="z [mm]  (beam axis)", ylabel="y [mm]",
           title=f"Scenario frame: {meta['phantom_material']} cylinder "
                 f"({2*R:.0f} mm x {2*half_z:.0f} mm), source overlaid")
    ax.set_aspect("equal")
    ax.set_xlim(-half_z - 28, half_z + 8)
    ax.set_ylim(-R - 8, R + 8)
    leg = ax.legend(loc="upper right", markerscale=4, framealpha=0.9, fontsize=8)
    for h in leg.legend_handles:
        if hasattr(h, "set_alpha"):
            h.set_alpha(1.0)
    fig.tight_layout()

    out = args.out or os.path.join(args.data_dir, "scenario_geometry.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
