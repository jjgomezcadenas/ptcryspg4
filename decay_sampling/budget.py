#!/usr/bin/env python3
"""Handoff budget (deterministic): measured decays N_j per isotope (Eq. 1).

Reads the Stage-A output, scales production to a clinical dose
(P_j(D)=count_j·D/target_dose), and applies the three-factor survival of Eq. 1 at
the operating point to get the expected measured decays N_j. This is the
detector-independent, RNG-free source budget — the thin quantity handed across the
A|B seam to the downstream detector study. The stochastic Poisson realizations
and the σ(range) figure of merit live downstream (budget_gen.py here for now;
moves there later). Writes:
  data/sampling_budget_<scenario>.csv       (isotope_id, N_expected)
  data/sampling_budget_<scenario>_meta.csv  (operating point, source)
See docs/handoff.tex.

Usage:
    python decay_sampling/budget.py [data_dir] [--scenario NAME] [--dose GY]
        [--t-irr S] [--t-del S] [--t-meas S]
"""

import argparse
import math
import os
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402


def survival(lam, t_irr, t_del, t_meas):
    """Eq. 1 factors (build-up, transport, window) for decay constant lam."""
    build = (1.0 - math.exp(-lam * t_irr)) / (lam * t_irr)
    transport = math.exp(-lam * t_del)
    window = 1.0 - math.exp(-lam * t_meas)
    return build, transport, window


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?", default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--scenario", default="inroom")
    ap.add_argument("--dose", type=float, default=1.0, help="delivered dose D [Gy]")
    ap.add_argument("--t-irr", type=float, default=60.0)
    ap.add_argument("--t-del", type=float, default=120.0)
    ap.add_argument("--t-meas", type=float, default=1200.0)
    args = ap.parse_args()

    emit = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))
    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    t_dose = float(meta["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()

    # Expected measured decays per isotope: N_j = P_j(D) · survival.
    rows = []
    print(f"\nscenario '{args.scenario}'  ({meta['phantom_material']}, "
          f"{args.dose:g} Gy; t_irr={args.t_irr:g}, t_del={args.t_del:g}, "
          f"t_meas={args.t_meas:g} s)")
    print(f"{'iso':>5} {'build':>6} {'transp':>7} {'window':>7} "
          f"{'P_j(D)':>10} {'N_j':>10}")
    print("-" * 50)
    for iid in sorted(ISOTOPES):
        lam = ISOTOPES[iid].lam
        pj = counts.get(iid, 0) * args.dose / t_dose
        b, tr, w = survival(lam, args.t_irr, args.t_del, args.t_meas)
        n_exp = pj * b * tr * w
        rows.append((iid, n_exp))
        print(f"{ISOTOPES[iid].name:>5} {b:>6.3f} {tr:>7.3f} {w:>7.3f} "
              f"{pj:>10.3e} {n_exp:>10.3e}")
    total = sum(r[1] for r in rows)
    print("-" * 50)
    print(f"{'total':>5} {'':>22} {'':>10} {total:>10.3e}")
    n_exp_by_id = dict(rows)
    o15, c11 = n_exp_by_id.get(0, 0.0), n_exp_by_id.get(1, 0.0)
    if c11:
        print(f"\nmeasured 15O/11C = {o15 / c11:.2f}")

    budget = pd.DataFrame(rows, columns=["isotope_id", "N_expected"])
    bpath = os.path.join(args.data_dir, f"sampling_budget_{args.scenario}.csv")
    budget.to_csv(bpath, index=False, float_format="%.6e")

    meta_out = {
        "scenario": args.scenario,
        "source_file": "emitters.csv",
        "dose_Gy": args.dose,
        "t_irr_s": args.t_irr, "t_del_s": args.t_del, "t_meas_s": args.t_meas,
        "target_dose_Gy": t_dose,
    }
    mpath = os.path.join(args.data_dir, f"sampling_budget_{args.scenario}_meta.csv")
    pd.DataFrame([meta_out]).to_csv(mpath, index=False)

    print(f"\nwrote budget -> {bpath}")
    print(f"      meta   -> {mpath}")


if __name__ == "__main__":
    main()
