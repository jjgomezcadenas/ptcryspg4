#!/usr/bin/env python3
"""Distal emitter-pool density: is a run dense enough for a sigma_R study?

A downstream sigma_R(distal) study draws its annihilation points, with
replacement, from emitters.csv -- so the spread it measures over realizations
is limited by how many DISTINCT emitters populate the distal falloff, the only
region the range fit uses. Run statistics do not change the measured decays N_j
(those are dose-normalized), only this pool granularity.

This reports, for the isotope that sets the range floor (O15 by default) and for
all isotopes, the emitter count within +-HALF cm of the activity distal edge,
and the resulting pool edge-position uncertainty ~ window / sqrt(n). With
--min-o15 it exits non-zero when the O15 count is below the target, so it can
gate a scenario before it is frozen.

The distal edge is R80 of the central-axis dose (sobp_metrics), the same depth
check_run gates on; --edge overrides it for a pencil run or a custom reference.

Usage:
    python analysis_transport/distal_pool.py <run_dir>
        [--half CM] [--isotope ID] [--min-o15 N]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402
from sobp_metrics import metrics as plateau_metrics  # noqa: E402


def distal_edge_cm(run_dir):
    """R80 of the central-axis dose [cm from entrance] -- the activity distal
    edge the range fit targets. Raises ValueError if the run cannot be graded."""
    r = plateau_metrics(run_dir)
    if r["r80_distal"] is None:
        raise ValueError(f"{run_dir}: no distal R80 crossing (not an SOBP run?)")
    return r["r80_distal"]


def pool(run_dir, edge_cm, half_cm, iso_id):
    """Emitter counts in the distal window [edge-half, edge+half] cm."""
    m = pd.read_csv(os.path.join(run_dir, "run_meta.csv")).iloc[0]
    e = pd.read_csv(os.path.join(run_dir, "emitters.csv"),
                    usecols=["isotope_id", "anh_z_mm"])
    half_z = 0.5 * float(m["phantom_length_mm"])
    depth = (e["anh_z_mm"].to_numpy() + half_z) / 10.0        # cm from entrance
    win = (depth >= edge_cm - half_cm) & (depth <= edge_cm + half_cm)
    iso = e["isotope_id"].to_numpy()
    per = {iid: int(((iso == iid) & win).sum()) for iid in sorted(ISOTOPES)}
    return per, int(win.sum()), int(((iso == iso_id) & win).sum())


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a Stage-A SOBP run directory")
    ap.add_argument("--half", type=float, default=0.5,
                    help="half-width of the distal window [cm] (default 0.5)")
    ap.add_argument("--isotope", type=int, default=0,
                    help="isotope_id that sets the range floor (default 0 = O15)")
    ap.add_argument("--edge", type=float, default=None,
                    help="distal edge [cm from entrance]; default R80 of the run")
    ap.add_argument("--min-o15", type=int, default=None,
                    help="fail if the target isotope's distal count is below this")
    args = ap.parse_args()

    try:
        edge = args.edge if args.edge is not None else distal_edge_cm(args.run_dir)
    except ValueError as ex:
        sys.exit(str(ex))

    per, total, target = pool(args.run_dir, edge, args.half, args.isotope)
    name = ISOTOPES[args.isotope].name
    # Pool edge-position uncertainty: the empirical falloff position averages
    # ~n points over the window, so its standard error scales as window/sqrt(n).
    sigma_edge_mm = (2 * args.half * 10.0) / np.sqrt(target) if target else float("inf")

    print(f"distal pool: {args.run_dir}")
    print(f"  distal edge (R80)    : {edge:.2f} cm from entrance")
    print(f"  window               : {edge - args.half:.2f} - {edge + args.half:.2f} cm "
          f"(+-{args.half:g} cm)")
    print(f"  emitters in window   : {total}  (all isotopes)")
    for iid in sorted(ISOTOPES):
        print(f"    {ISOTOPES[iid].name:>4} : {per[iid]}")
    print(f"  {name} pool edge sigma : {sigma_edge_mm:.3f} mm  "
          f"(~ window / sqrt(N_{name}))")

    if args.min_o15 is not None:
        ok = target >= args.min_o15
        print(f"\n{'PASS' if ok else 'FAIL'}: {name} distal pool {target} "
              f"{'>=' if ok else '<'} {args.min_o15}")
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
