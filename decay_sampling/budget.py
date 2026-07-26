#!/usr/bin/env python3
"""Handoff budget (deterministic): measured decays N_j per isotope.

Reads the Stage-A output, scales production to a clinical dose
(P_j(D)=count_j·D/target_dose), and applies the three-factor survival ---
build-up during irradiation, decay over the beam-off delay, decay inside the
measurement window; derived in ptcrysp_physics.pdf, "Decay kinetics" --- at
the operating point to get the expected measured decays N_j. This is the
detector-independent, RNG-free source budget — the thin quantity handed to the
downstream detector study. The stochastic Poisson realizations
and the σ(range) figure of merit live downstream (budget_gen.py here for now;
moves there later). Writes:
  data/sampling_budget_<scenario>.csv       (isotope_id, N_expected)
  data/sampling_budget_<scenario>_meta.csv  (operating point, source)
See latex/ptcrysp_physics.tex.

Usage:
    python decay_sampling/budget.py <run_dir> [--scenario NAME] [--dose GY]
        [--scenario-config FILE]
"""

import argparse
import os
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _REPO)
from common.isotopes import ISOTOPES  # noqa: E402
from decay_sampling.scenarios import (  # noqa: E402
    DEFAULT_SCENARIO_CONFIG, resolve_scenario, survival,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", help="a run directory, e.g. data/runs/cylinder_sobp_1e7")
    ap.add_argument("--scenario", default="d120s300")
    ap.add_argument("--scenario-config", default=str(DEFAULT_SCENARIO_CONFIG))
    ap.add_argument("--dose", type=float, default=1.0, help="delivered dose D [Gy]")
    args = ap.parse_args()
    scenario = resolve_scenario(args.scenario, args.scenario_config)

    emit = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"),
                       usecols=["isotope_id"])  # counts only — skip the positions
    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    t_dose = float(meta["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()

    # Expected measured decays per isotope: N_j = P_j(D) · survival.
    rows = []
    print(f"\nscenario '{args.scenario}'  ({meta['phantom_material']}, "
          f"{args.dose:g} Gy; t_irr={scenario.t_irr_s:g}, "
          f"t_del={scenario.t_del_s:g}, t_meas={scenario.t_meas_s:g} s)")
    print(f"{'iso':>5} {'build':>6} {'transp':>7} {'window':>7} "
          f"{'P_j(D)':>10} {'N_j':>10}")
    print("-" * 50)
    for iid in sorted(ISOTOPES):
        lam = ISOTOPES[iid].lam
        pj = counts.get(iid, 0) * args.dose / t_dose
        b, tr, w = scenario.factors(lam)
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
        "scenario_config": str(scenario.config_path),
        "scenario_config_sha256": scenario.config_sha256,
        "source_file": "emitters.csv",
        "dose_Gy": args.dose,
        "t_irr_s": scenario.t_irr_s,
        "t_del_s": scenario.t_del_s,
        "t_meas_s": scenario.t_meas_s,
        "target_dose_Gy": t_dose,
    }
    mpath = os.path.join(args.data_dir, f"sampling_budget_{args.scenario}_meta.csv")
    pd.DataFrame([meta_out]).to_csv(mpath, index=False)

    print(f"\nwrote budget -> {bpath}")
    print(f"      meta   -> {mpath}")


if __name__ == "__main__":
    main()
