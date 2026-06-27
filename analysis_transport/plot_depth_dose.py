#!/usr/bin/env python3
"""Central-axis depth dose from depth_dose.csv (Step 1): the clean on-axis core
vs the contaminated full-plane tally.

Left: dose_core_Gy along the beam (dose to the medium on the axis) with the
target window shaded — the physical central-axis depth dose graded for the SOBP
(R80, plateau flatness come in Step 3). Right: the shape comparison, edep_total
(whole transverse plane) vs edep_core (r <= 5 mm), each normalised to its own
peak — through the head the full-plane curve carries the ellipsoid's varying
cross-section and the bone shells, the core does not. Writes
<run_dir>/figures/depth_dose.png.

Usage:
    python analysis_transport/plot_depth_dose.py <run_dir> [--out PNG]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a Stage-A run directory")
    ap.add_argument("--out", default=None,
                    help="output PNG (default <run_dir>/figures/depth_dose.png)")
    args = ap.parse_args()

    d = pd.read_csv(os.path.join(args.run_dir, "depth_dose.csv"))
    meta = pd.read_csv(os.path.join(args.run_dir, "run_meta.csv")).iloc[0]
    half = 0.5 * float(meta["phantom_length_mm"])
    depth = (d["z_mm"] + half) / 10.0  # cm from the entrance face
    prox = float(meta["target_prox_depth_mm"]) / 10.0
    dist = float(meta["target_dist_depth_mm"]) / 10.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: the physical on-axis dose.
    ax1.plot(depth, d["dose_core_Gy"], color="C0", lw=1.6)
    ax1.axvspan(prox, dist, color="C2", alpha=0.15, label="target")
    ax1.set(xlabel="depth from entrance [cm]", ylabel="dose_core [Gy]",
            title="Central-axis dose (on-axis core, r ≤ 5 mm)")
    ax1.set_xlim(0, 2 * half / 10.0)
    ax1.set_ylim(bottom=0)
    ax1.legend()

    # Right: shape comparison, each normalised to its own peak.
    for col, c, lab in [("edep_total_MeV", "C3", "full plane (edep_total)"),
                        ("edep_core_MeV", "C0", "on-axis core (edep_core)")]:
        peak = d[col].max() or 1.0
        ax2.plot(depth, d[col] / peak, color=c, lw=1.5, label=lab)
    ax2.axvspan(prox, dist, color="C2", alpha=0.15)
    ax2.set(xlabel="depth from entrance [cm]", ylabel="normalised to own peak",
            title="Shape: full-plane vs on-axis core")
    ax2.set_xlim(0, 2 * half / 10.0)
    ax2.set_ylim(bottom=0)
    ax2.legend()

    fig.suptitle(f"Depth dose: {meta['geometry']} ({meta['phantom_material']}), "
                 f"{int(meta['n_protons']):g} protons", fontsize=12, weight="bold")
    fig.tight_layout()

    if args.out:
        out = args.out
    else:
        figdir = os.path.join(args.run_dir, "figures")
        os.makedirs(figdir, exist_ok=True)
        out = os.path.join(figdir, "depth_dose.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
