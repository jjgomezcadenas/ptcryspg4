#!/usr/bin/env python3
"""SOBP plateau metrics on a run's central-axis depth dose: R80 + uniformity.

Reads <run_dir>/depth_dose.csv (the dose_core_Gy column from Step 1) and
run_meta.csv, and reports the standard acceptance numbers for a spread-out Bragg
peak through the phantom:

  * plateau uniformity over the target window (excluding the inherent distal
    falloff) -- (max-min)/mean of dose_core;
  * R80 -- the distal depth at 80% of the plateau dose (nearly independent of the
    energy spread, so the robust range measure), and its offset from the target's
    distal edge (the bone-offset check: a WEPL-designed head field lands R80 at
    the target edge, a water-tuned one falls ~1 cm short);
  * the proximal 80% rise and the resulting 80%-80% modulation width.

Graded on dose_core (the on-axis core, edep proportional to dose in constant-
density material), not the cross-section-contaminated full-plane edep. With
acceptance flags (--max-uniformity, --r80-tol) it exits non-zero on a fail, so it
can gate a run; without them it just reports.

Usage:
    python analysis_transport/sobp_metrics.py <run_dir>
        [--max-uniformity PCT] [--r80-tol MM]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

FALLOFF_CM = 0.5  # excluded from the plateau: any SOBP must drop to zero here


def crossing(depth, dose, level, lo, hi):
    """Depth where dose first crosses `level` within [lo, hi], interpolated;
    None if it does not cross in that span."""
    m = (depth >= lo) & (depth <= hi)
    d, y = depth[m], dose[m]
    below = y < level
    flips = np.where(below[:-1] != below[1:])[0]
    if len(flips) == 0:
        return None
    i = flips[0] if y[0] >= level else flips[-1]  # proximal rise vs distal fall
    return float(np.interp(level, [y[i], y[i + 1]][:: 1 if y[i] < y[i + 1] else -1],
                           [d[i], d[i + 1]][:: 1 if y[i] < y[i + 1] else -1]))


def metrics(run_dir):
    d = pd.read_csv(os.path.join(run_dir, "depth_dose.csv"))
    if "dose_core_Gy" not in d.columns:
        sys.exit(f"{run_dir}/depth_dose.csv has no dose_core_Gy "
                 "(re-run Stage A; needs the Step-1 central-axis tally)")
    m = pd.read_csv(os.path.join(run_dir, "run_meta.csv")).iloc[0]
    half = 0.5 * float(m["phantom_length_mm"])
    depth = (d["z_mm"].to_numpy() + half) / 10.0          # cm from the entrance
    dose = d["dose_core_Gy"].to_numpy()
    prox = float(m["target_prox_depth_mm"]) / 10.0
    dist = float(m["target_dist_depth_mm"]) / 10.0

    plat = dose[(depth >= prox) & (depth <= dist - FALLOFF_CM)]
    plateau = float(plat.mean())
    uniformity = 100.0 * (plat.max() - plat.min()) / plateau

    mid = 0.5 * (prox + dist)
    r80_distal = crossing(depth, dose, 0.8 * plateau, mid, depth.max())
    r80_prox = crossing(depth, dose, 0.8 * plateau, 0.0, mid)
    return dict(geometry=str(m["geometry"]), plateau_Gy=plateau,
                uniformity_pct=uniformity, prox=prox, dist=dist,
                r80_distal=r80_distal, r80_prox=r80_prox)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a Stage-A SOBP run directory")
    ap.add_argument("--max-uniformity", type=float, default=None,
                    help="fail if plateau uniformity exceeds this %% (acceptance)")
    ap.add_argument("--r80-tol", type=float, default=None,
                    help="fail if |R80 - target distal edge| exceeds this mm")
    args = ap.parse_args()

    r = metrics(args.run_dir)
    print(f"plateau metrics: {args.run_dir}  ({r['geometry']})")
    print(f"  target window        : {r['prox']:.2f} - {r['dist']:.2f} cm")
    print(f"  plateau dose (core)  : {r['plateau_Gy']:.3e} Gy")
    print(f"  plateau uniformity   : {r['uniformity_pct']:.1f} %  "
          f"(over {r['prox']:.1f}-{r['dist']-FALLOFF_CM:.1f} cm)")
    if r["r80_prox"] is not None and r["r80_distal"] is not None:
        print(f"  80% modulation width : {r['r80_prox']:.2f} - {r['r80_distal']:.2f} cm "
              f"({r['r80_distal']-r['r80_prox']:.2f} cm)")
    fails = []
    if r["r80_distal"] is not None:
        off_mm = (r["r80_distal"] - r["dist"]) * 10.0
        print(f"  R80 (distal range)   : {r['r80_distal']:.2f} cm  "
              f"({off_mm:+.0f} mm vs target distal edge {r['dist']:.1f} cm)")
        if args.r80_tol is not None and abs(off_mm) > args.r80_tol:
            fails.append(f"R80 off by {off_mm:+.0f} mm (> {args.r80_tol:g})")
    else:
        print("  R80 (distal range)   : no distal 80% crossing found")
    if args.max_uniformity is not None and r["uniformity_pct"] > args.max_uniformity:
        fails.append(f"uniformity {r['uniformity_pct']:.1f}% (> {args.max_uniformity:g}%)")

    if fails:
        print("\nFAIL: " + "; ".join(fails))
        sys.exit(1)
    if args.max_uniformity is not None or args.r80_tol is not None:
        print("\nacceptance: PASS")


if __name__ == "__main__":
    main()
