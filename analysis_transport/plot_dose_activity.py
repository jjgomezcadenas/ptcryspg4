#!/usr/bin/env python3
"""Overlay of depth dose and beta+ activity along the beam, in arbitrary units.

One axes, two curves vs depth from the entrance face: the central-axis dose
(dose_core_Gy, the Bragg curve) and the beta+ activity profile (a histogram of
annihilation points along the beam). Both are normalised to their own peak, so
the reader can see directly how the activity compares with the dose: the activity
is roughly flat along the proton track and falls away at the Bragg peak, where the
protons drop below the nuclear-reaction thresholds, while the dose peaks there.
Writes <run_dir>/figures/dose_activity.png.

Usage:
    python analysis_transport/plot_dose_activity.py <run_dir> [--out PNG] [--bin-mm MM]
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a Stage-A run directory")
    ap.add_argument("--out", default=None,
                    help="output PNG (default <run_dir>/figures/dose_activity.png)")
    ap.add_argument("--bin-mm", type=float, default=2.0,
                    help="activity histogram bin width [mm] (default 2)")
    args = ap.parse_args()

    dd = pd.read_csv(os.path.join(args.run_dir, "depth_dose.csv"))
    emit = pd.read_csv(os.path.join(args.run_dir, "emitters.csv"))
    meta = pd.read_csv(os.path.join(args.run_dir, "run_meta.csv")).iloc[0]
    half = 0.5 * float(meta["phantom_length_mm"])

    # Depth from the entrance face (cm), for both the dose curve and the emitters.
    depth_dose = (dd["z_mm"].to_numpy() + half) / 10.0
    dose = dd["dose_core_Gy"].to_numpy()
    depth_anh = (emit["anh_z_mm"].to_numpy() + half) / 10.0

    # beta+ activity profile: annihilation points binned along the beam.
    edges_mm = np.arange(0.0, 2.0 * half + args.bin_mm, args.bin_mm)
    counts, edges = np.histogram(depth_anh * 10.0, bins=edges_mm)
    centres = 0.5 * (edges[:-1] + edges[1:]) / 10.0  # cm

    # Arbitrary units: each curve to its own peak so both fit the same frame.
    dose_n = dose / (dose.max() or 1.0)
    act_n = counts / (counts.max() or 1.0)

    # Show the beam up to a little past the Bragg peak (where both have died away).
    peak_depth = depth_dose[int(np.argmax(dose))]
    xmax = min(2.0 * half / 10.0, peak_depth + 2.5)

    # Beam label: an SOBP run carries a layer table; otherwise it is a pencil.
    is_sobp = os.path.exists(os.path.join(args.run_dir, "sobp_layers.csv"))
    beam = "SOBP" if is_sobp else f"{float(meta['beam_energy_MeV']):g} MeV pencil"

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.plot(depth_dose, dose_n, color="C3", lw=1.8, label="dose")
    ax.step(centres, act_n, where="mid", color="C0", lw=1.6,
            label=r"$\beta^+$ activity")
    ax.set(xlabel="depth from entrance [cm]", ylabel="arbitrary units (peak = 1)",
           xlim=(0, xmax), ylim=(0, 1.08))
    ax.set_title(f"Dose vs $\\beta^+$ activity: {meta['geometry']} "
                 f"({meta['phantom_material']}), {beam}, "
                 f"{int(meta['n_protons']):g} protons", fontsize=11)
    ax.legend()
    fig.tight_layout()

    if args.out:
        out = args.out
    else:
        figdir = os.path.join(args.run_dir, "figures")
        os.makedirs(figdir, exist_ok=True)
        out = os.path.join(figdir, "dose_activity.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
