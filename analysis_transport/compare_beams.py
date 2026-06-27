#!/usr/bin/env python3
"""Compare central-axis depth dose across runs: single pencil vs SOBP, per phantom.

Overlays the on-axis `dose_core` of several runs (each normalised to its own peak,
so different statistics are comparable as shapes), one panel per geometry. A run
that carries an SOBP field (sobp_layers_meta) is labelled SOBP, else pencil. Shows
why a spread-out field is needed (a pencil peaks at one depth; the SOBP fills the
target) and, across the head panels, the bone's effect.

Usage:
    python analysis_transport/compare_beams.py RUN_DIR [RUN_DIR ...] --out PNG
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def load(run_dir):
    d = pd.read_csv(os.path.join(run_dir, "depth_dose.csv"))
    m = pd.read_csv(os.path.join(run_dir, "run_meta.csv")).iloc[0]
    half = 0.5 * float(m["phantom_length_mm"])
    depth = (d["z_mm"].to_numpy() + half) / 10.0  # cm from entrance
    dose = d["dose_core_Gy"].to_numpy()
    is_sobp = os.path.exists(os.path.join(run_dir, "sobp_layers_meta.csv"))
    return dict(geometry=str(m["geometry"]), depth=depth, dose=dose,
                beam="SOBP" if is_sobp else "pencil",
                prox=float(m["target_prox_depth_mm"]) / 10.0,
                dist=float(m["target_dist_depth_mm"]) / 10.0)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    runs = [load(r) for r in args.run_dirs]
    geoms = list(dict.fromkeys(r["geometry"] for r in runs))  # preserve order
    fig, axes = plt.subplots(1, len(geoms), figsize=(6.2 * len(geoms), 4.6),
                             squeeze=False)
    style = {"pencil": dict(color="C3", lw=1.4, ls="--"),
             "SOBP": dict(color="C0", lw=1.7)}
    for ax, g in zip(axes[0], geoms):
        prox = dist = None
        for r in [x for x in runs if x["geometry"] == g]:
            ax.plot(r["depth"], r["dose"] / r["dose"].max(),
                    label=r["beam"], **style[r["beam"]])
            prox, dist = r["prox"], r["dist"]
        ax.axvspan(prox, dist, color="C2", alpha=0.15, label="target")
        ax.set(xlabel="depth from entrance [cm]", ylabel="dose (norm. to peak)",
               title=g, xlim=(0, dist + 3), ylim=(0, 1.08))
        ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
