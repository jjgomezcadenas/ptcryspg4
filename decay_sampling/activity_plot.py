#!/usr/bin/env python3
"""Activity vs time after irradiation, per isotope (cf. CRYSP-paper Fig. 3).

A_j(t) = A_peak_j · exp(-lambda_j·t), where t is the time after the end of the
beam and A_peak_j = (P_j/t_irr)(1 - exp(-lambda_j·t_irr)) is the activity at beam
end. The absolute production is P_j(D) = count_j · D / target_dose_Gy from the
Stage-A run (docs/handoff.tex). The integral of A_j over the acquisition window
equals the measured decays N_j of Eq. 1.

Usage:
    python decay_sampling/activity_plot.py [data_dir] [--dose GY] [--t-irr S]
        [--t-del S] [--t-meas S] [--t-max MIN]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?", default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--dose", type=float, default=1.0, help="delivered dose D [Gy]")
    ap.add_argument("--t-irr", type=float, default=60.0, help="irradiation time [s]")
    ap.add_argument("--t-del", type=float, default=120.0, help="transport delay [s]")
    ap.add_argument("--t-meas", type=float, default=1200.0, help="acquisition window [s]")
    ap.add_argument("--t-max", type=float, default=45.0, help="plot range [min]")
    args = ap.parse_args()

    emit = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))
    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    t_dose = float(meta["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()

    t = np.linspace(0.0, args.t_max * 60.0, 2000)  # s after end of beam
    fig, ax = plt.subplots(figsize=(8, 5))
    total = np.zeros_like(t)
    for iid in sorted(ISOTOPES):
        iso = ISOTOPES[iid]
        lam = iso.lam
        Pj = counts.get(iid, 0) * args.dose / t_dose            # produced decays
        a_peak = (Pj / args.t_irr) * (1.0 - np.exp(-lam * args.t_irr))  # Bq
        a_mbq = a_peak * np.exp(-lam * t) / 1e6
        total += a_mbq
        ax.plot(t / 60.0, a_mbq, label=iso.name)
    ax.plot(t / 60.0, total, "k", lw=2.0, label="total")

    # Acquisition window (starts t_del after beam end), and the in-room start.
    ax.axvspan(args.t_del / 60.0, (args.t_del + args.t_meas) / 60.0,
               color="gray", alpha=0.12, label="acquisition")
    ax.axvline(args.t_del / 60.0, color="gray", ls="--", lw=1)

    ax.set(xlabel="time after irradiation [min]", ylabel="activity [MBq]",
           title=f"Induced activity ({meta['phantom_material']}, {args.dose:g} Gy)")
    ax.set_ylim(bottom=0)
    ax.legend()
    fig.tight_layout()
    out = os.path.join(args.data_dir, "activity.png")
    fig.savefig(out, dpi=120)

    # Quick numeric readout at a few times (cf. the paper's table).
    print(f"\nactivity [MBq] ({meta['phantom_material']}, {args.dose:g} Gy):")
    print(f"{'t[min]':>7} " + " ".join(f"{ISOTOPES[i].name:>7}" for i in sorted(ISOTOPES))
          + f"{'total':>8}")
    for tmin in (0, 2, 10, 20):
        row = []
        tot = 0.0
        for iid in sorted(ISOTOPES):
            lam = ISOTOPES[iid].lam
            Pj = counts.get(iid, 0) * args.dose / t_dose
            a0 = (Pj / args.t_irr) * (1 - np.exp(-lam * args.t_irr)) / 1e6
            a = a0 * np.exp(-lam * tmin * 60)
            row.append(a)
            tot += a
        print(f"{tmin:>7} " + " ".join(f"{v:>7.3f}" for v in row) + f"{tot:>8.3f}")
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
