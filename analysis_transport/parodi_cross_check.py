#!/usr/bin/env python3
"""Cross-check Stage-A absolute β+ yields against Parodi et al. 2008 (Table 2).

For the Parodi standard scenario (proton SOBP delivering 1 Gy to a target in a
head), the absolute integral production is fixed by the dose normalization:

    P_j(1 Gy) = count_j / target_dose_Gy          (see CLAUDE.md)

This prints our P_j against Parodi's published head-field yields per Gy.
Agreement to ~2x for the dominant isotopes is the expected "gross" check —
G4 cross-section bias plus our homogeneous brain vs the real heterogeneous head
(brain is more oxygen-rich, so it overproduces 15O).

Usage:
    python analysis_transport/parodi_cross_check.py [data_dir]
"""

import os
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402

# Parodi et al., Int. J. Radiat. Oncol. Biol. Phys. 71, 945 (2008), Table 2:
# integral beta+ yield per 1 Gy proton field on the HEAD (skull-base).
# Reference literature constants, keyed by isotope_id (see common/SCHEMA.md).
PARODI_HEAD_PER_GY = {
    0: 9.4e7,  # 15O
    1: 7.7e7,  # 11C
    2: 9.4e6,  # 13N
    3: 1.6e6,  # 10C
    4: 5.5e5,  # 14O
}


def main(data_dir: str) -> None:
    emit = pd.read_csv(os.path.join(data_dir, "emitters.csv"))
    meta = pd.read_csv(os.path.join(data_dir, "run_meta.csv")).iloc[0]
    t_dose = float(meta["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()

    print(f"\nmaterial: {meta['phantom_material']}   "
          f"target_dose: {t_dose:.3e} Gy   "
          f"N_p(1 Gy): {meta['Np_per_Gy']:.3e}")
    if meta["phantom_material"] != "G4_BRAIN_ICRP":
        print("  (NOTE: Parodi Table 2 is a brain/head field; cross-check is "
              "only meaningful for G4_BRAIN_ICRP)")

    print(f"\n{'iso':>5} {'ours/Gy':>12} {'Parodi/Gy':>12} {'ratio':>7}")
    print("-" * 40)
    tot_ours = tot_parodi = 0.0
    for iid in sorted(ISOTOPES):
        name = ISOTOPES[iid].name
        ours = counts.get(iid, 0) / t_dose          # P_j(1 Gy)
        parodi = PARODI_HEAD_PER_GY.get(iid, float("nan"))
        ratio = ours / parodi if parodi else float("nan")
        print(f"{name:>5} {ours:>12.3e} {parodi:>12.3e} {ratio:>6.1f}x")
        tot_ours += ours
        tot_parodi += parodi
    print("-" * 40)
    print(f"{'total':>5} {tot_ours:>12.3e} {tot_parodi:>12.3e} "
          f"{tot_ours / tot_parodi:>6.1f}x")
    print("\n(~2x on 15O/11C/total is the expected gross agreement.)")


if __name__ == "__main__":
    ddir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "..", "data")
    main(ddir)
