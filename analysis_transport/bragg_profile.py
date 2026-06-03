#!/usr/bin/env python3
"""Terminal quick-look at the Stage-A proton Bragg depth-dose profile.

Reads data/depth_dose.csv (written by the proton_transport app) and prints the
energy-deposit-vs-depth curve as an ASCII bar chart, with the Bragg-peak
location and the primary/secondary split. A fast sanity check that the beam
ranges out where expected and the depth profile is smooth -- complements the
PNG plots from validate_transport.py.

Usage:
    python analysis_transport/bragg_profile.py [data_dir] [--step MM] [--width N]
"""

import argparse
import os

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("data_dir", nargs="?",
                    default=os.path.join(_HERE, "..", "data"),
                    help="directory holding depth_dose.csv")
    ap.add_argument("--step", type=float, default=2.0,
                    help="depth spacing between printed rows [mm]")
    ap.add_argument("--width", type=int, default=48,
                    help="bar width in characters at the peak")
    args = ap.parse_args()

    d = pd.read_csv(os.path.join(args.data_dir, "depth_dose.csv"))
    tot, prim = d["edep_total_MeV"], d["edep_primary_MeV"]

    peak_i = tot.idxmax()
    peak_z = d["z_mm"][peak_i]
    peak_v = tot[peak_i]
    sum_tot, sum_prim = tot.sum(), prim.sum()
    prim_frac = sum_prim / sum_tot if sum_tot else 0.0

    # Distal edge: deepest bin still above 50% of the peak (range proxy).
    half = peak_v * 0.5
    distal = d["z_mm"][tot >= half].max()

    print(f"\nBragg depth-dose profile  ({os.path.abspath(args.data_dir)})")
    print(f"  bins            : {len(d)}  over z = "
          f"[{d['z_mm'].min():.1f}, {d['z_mm'].max():.1f}] mm")
    print(f"  Bragg peak      : z = {peak_z:.1f} mm   ({peak_v:.1f} MeV/bin)")
    print(f"  distal edge @50%: z = {distal:.1f} mm")
    print(f"  energy split    : {100 * prim_frac:.2f}% primary proton, "
          f"{100 * (1 - prim_frac):.2f}% secondaries\n")

    bin_w = d["z_mm"].iloc[1] - d["z_mm"].iloc[0]
    stride = max(1, int(round(args.step / bin_w)))
    scale = args.width / peak_v if peak_v else 0.0

    for _, r in d[::stride].iterrows():
        # Skip the near-zero noise past the distal fall-off.
        if r["edep_total_MeV"] < 0.01 * peak_v and r["z_mm"] > peak_z:
            continue
        bar = "#" * int(scale * r["edep_total_MeV"])
        mark = "  <-- peak" if abs(r["z_mm"] - peak_z) < bin_w * stride / 2 else ""
        print(f"  z={r['z_mm']:7.1f}  {r['edep_total_MeV']:9.1f} MeV  "
              f"{bar}{mark}")
    print()


if __name__ == "__main__":
    main()
