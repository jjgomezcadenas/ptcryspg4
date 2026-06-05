#!/usr/bin/env python3
"""Handoff realizations (stochastic): Z Poisson draws M_j ~ Poisson(N_j).

Reads the deterministic budget written by budget.py (N_expected per isotope) and
draws Z independent Poisson realizations M_j^(z) ~ Poisson(N_expected). Each
realization is one "experiment": PTCryspMC.jl samples M_j^(z) annihilation points
from emitters.csv (seed master_seed + realization, so every detector sees the
identical source), reconstructs, and fits the range; the spread of the fitted
ranges over the Z realizations is σ(range). See docs/handoff.tex.

This is the stochastic, detector-side half of the old budget.py. It reads exactly
the file PTCryspMC.jl/py/ will read, so relocating it there is a pure move.
Writes:
  data/sampling_realizations_<scenario>.csv       (realization, isotope_id, N_poisson)
  data/sampling_realizations_<scenario>_meta.csv  (n_realizations, master_seed, source)

Usage:
    python decay_sampling/budget_gen.py [data_dir] [--scenario NAME]
        [--realizations Z] [--seed S]
"""

import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    _here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("data_dir", nargs="?", default=os.path.join(_here, "..", "data"))
    ap.add_argument("--scenario", default="inroom")
    ap.add_argument("--realizations", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    bpath = os.path.join(args.data_dir, f"sampling_budget_{args.scenario}.csv")
    budget = pd.read_csv(bpath)
    n_exp = dict(zip(budget["isotope_id"], budget["N_expected"]))

    rng = np.random.default_rng(args.seed)
    rows = []
    for z in range(args.realizations):
        for iid in sorted(n_exp):
            m = int(rng.poisson(n_exp[iid]))
            rows.append((z, iid, m))
    real = pd.DataFrame(rows, columns=["realization", "isotope_id", "N_poisson"])
    rpath = os.path.join(args.data_dir, f"sampling_realizations_{args.scenario}.csv")
    real.to_csv(rpath, index=False)

    meta_out = {
        "scenario": args.scenario,
        "source_budget": f"sampling_budget_{args.scenario}.csv",
        "n_realizations": args.realizations,
        "master_seed": args.seed,
    }
    mpath = os.path.join(args.data_dir,
                         f"sampling_realizations_{args.scenario}_meta.csv")
    pd.DataFrame([meta_out]).to_csv(mpath, index=False)

    print(f"scenario '{args.scenario}': {args.realizations} realizations "
          f"x {len(n_exp)} isotopes")
    print(f"wrote realizations -> {rpath}")
    print(f"      meta         -> {mpath}")


if __name__ == "__main__":
    main()
