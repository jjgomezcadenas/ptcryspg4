#!/usr/bin/env python3
"""Overlay of depth dose and beta+ activity along the beam, in arbitrary units.

Per run: one axes, two curves vs depth from the entrance face -- the central-axis
dose (dose_core_Gy, the Bragg curve / SOBP plateau) and the beta+ activity profile
(a histogram of annihilation points along the beam). Both are normalised to their
own peak, so the reader sees directly how the activity compares with the dose: it
is roughly flat along the proton track and falls to zero just before the distal
dose edge, where the protons drop below the nuclear-reaction thresholds. That
distal fall-off, A(x) -> 0, is the estimator of the beam's distal reach.

Pass one run for a single overlay; pass several (e.g. a pencil and an SOBP run) to
get them side by side in one figure for comparison. Writes
<run_dir>/figures/dose_activity.png (single run) or the path given by --out.

Usage:
    python analysis_transport/plot_dose_activity.py <run_dir> [<run_dir> ...]
        [--out PNG] [--bin-mm MM]
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def draw_overlay(ax, run_dir, bin_mm=2.0):
    """Draw the dose / beta+ activity overlay for one run on the given axes."""
    dd = pd.read_csv(os.path.join(run_dir, "depth_dose.csv"))
    emit = pd.read_csv(os.path.join(run_dir, "emitters.csv"))
    meta = pd.read_csv(os.path.join(run_dir, "run_meta.csv")).iloc[0]
    half = 0.5 * float(meta["phantom_length_mm"])

    # Depth from the entrance face (cm), for the dose curve and the emitters.
    depth_dose = (dd["z_mm"].to_numpy() + half) / 10.0
    dose = dd["dose_core_Gy"].to_numpy()

    # beta+ activity profile: annihilation points binned along the beam.
    edges_mm = np.arange(0.0, 2.0 * half + bin_mm, bin_mm)
    counts, edges = np.histogram(emit["anh_z_mm"].to_numpy() + half, bins=edges_mm)
    centres = 0.5 * (edges[:-1] + edges[1:]) / 10.0  # cm

    # Arbitrary units: each curve to its own peak so both fit the same frame.
    dose_n = dose / (dose.max() or 1.0)
    act_n = counts / (counts.max() or 1.0)

    # Show the beam up to a little past the distal dose edge (where both have died).
    above = np.where(dose_n > 0.05)[0]
    distal = depth_dose[above[-1]] if len(above) else depth_dose[-1]
    xmax = min(2.0 * half / 10.0, distal + 2.0)

    # Beam label: an SOBP run carries a layer table; otherwise it is a pencil.
    is_sobp = os.path.exists(os.path.join(run_dir, "sobp_layers.csv"))
    beam = "SOBP" if is_sobp else f"{float(meta['beam_energy_MeV']):g} MeV pencil"

    ax.plot(depth_dose, dose_n, color="C3", lw=1.8, label="dose")
    ax.step(centres, act_n, where="mid", color="C0", lw=1.6,
            label=r"$\beta^+$ activity")
    ax.set(xlabel="depth from entrance [cm]", ylabel="arbitrary units (peak = 1)",
           xlim=(0, xmax), ylim=(0, 1.08))
    ax.set_title(f"{beam} in {meta['phantom_material']}", fontsize=11)
    ax.legend()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+", help="one or more Stage-A run directories")
    ap.add_argument("--out", default=None,
                    help="output PNG (default <first run>/figures/dose_activity.png)")
    ap.add_argument("--bin-mm", type=float, default=2.0,
                    help="activity histogram bin width [mm] (default 2)")
    args = ap.parse_args()

    n = len(args.run_dirs)
    fig, axes = plt.subplots(1, n, figsize=(8.5 * n if n > 1 else 8.5, 5.0),
                             squeeze=False)
    for ax, run_dir in zip(axes[0], args.run_dirs):
        draw_overlay(ax, run_dir, args.bin_mm)
    fig.tight_layout()

    if args.out:
        out = args.out
    else:
        figdir = os.path.join(args.run_dirs[0], "figures")
        os.makedirs(figdir, exist_ok=True)
        out = os.path.join(figdir, "dose_activity.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
