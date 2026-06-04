#!/usr/bin/env python3
"""Handoff sampling budget: measured decays N_j (Eq. 1) and Poisson realizations.

Reads the Stage-A output, scales production to a clinical dose
(P_j(D)=count_j·D/target_dose), applies the three-factor survival of Eq. 1 at the
operating point to get the expected measured decays N_j, and draws Z Poisson
realizations M_j^(z) ~ Poisson(N_j). Writes:
  data/sampling_budget_<scenario>.csv       (realization, isotope_id,
                                              N_expected, N_poisson)
  data/sampling_budget_<scenario>_meta.csv  (operating point, source, seed)
Stage B draws M_j^(z) annihilation points from emitters.csv with seed
(master_seed + realization), so every detector sees the identical source.
See docs/handoff.tex.

Usage:
    python decay_sampling/budget.py [data_dir] [--scenario NAME] [--dose GY]
        [--t-irr S] [--t-del S] [--t-meas S] [--realizations Z] [--seed S]
"""

import argparse
import math
import os
import sys

import numpy as np
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
    ap.add_argument("--realizations", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    emit = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))
    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    t_dose = float(meta["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()

    # Expected measured decays per isotope: N_j = P_j(D) · survival.
    n_exp = {}
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
        n_exp[iid] = pj * b * tr * w
        print(f"{ISOTOPES[iid].name:>5} {b:>6.3f} {tr:>7.3f} {w:>7.3f} "
              f"{pj:>10.3e} {n_exp[iid]:>10.3e}")
    total = sum(n_exp.values())
    print("-" * 50)
    print(f"{'total':>5} {'':>22} {'':>10} {total:>10.3e}")
    o15, c11 = n_exp.get(0, 0.0), n_exp.get(1, 0.0)
    if c11:
        print(f"\nmeasured 15O/11C = {o15 / c11:.2f}")

    # Z Poisson realizations.
    rng = np.random.default_rng(args.seed)
    rows = []
    for z in range(args.realizations):
        for iid in sorted(ISOTOPES):
            m = int(rng.poisson(n_exp[iid]))
            rows.append((z, iid, n_exp[iid], m))
    budget = pd.DataFrame(rows, columns=["realization", "isotope_id",
                                         "N_expected", "N_poisson"])
    bpath = os.path.join(args.data_dir, f"sampling_budget_{args.scenario}.csv")
    budget.to_csv(bpath, index=False, float_format="%.6e")

    meta_out = {
        "scenario": args.scenario,
        "source_file": "emitters.csv",
        "dose_Gy": args.dose,
        "t_irr_s": args.t_irr, "t_del_s": args.t_del, "t_meas_s": args.t_meas,
        "n_realizations": args.realizations,
        "master_seed": args.seed,
        "target_dose_Gy": t_dose,
    }
    mpath = os.path.join(args.data_dir, f"sampling_budget_{args.scenario}_meta.csv")
    pd.DataFrame([meta_out]).to_csv(mpath, index=False)

    print(f"\nwrote {args.realizations} realizations -> {bpath}")
    print(f"      meta -> {mpath}")


if __name__ == "__main__":
    main()
