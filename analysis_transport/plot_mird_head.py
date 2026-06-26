#!/usr/bin/env python3
"""Visualize the MIRD-head run: phantom, beam, and emitter paths.

Reads the medium geometry from phantom_regions.csv (Phase 2) -- world-frame,
axis-aligned ellipsoids, priority-ordered -- so nothing about the head is
hardcoded here. Renders, from a head run's emitters.csv:
  (A) the labelled phantom + beam: scalp / skull / brain cross-sections, filled;
  (B) the emitter trail: production points coloured by isotope over the outlines;
  (C) the isotope mix per region -- bone shifts the local O15/C11 balance;
  (D) the production profile along the beam.

Usage:
    python analysis_transport/plot_mird_head.py [data_dir] [--out PNG]
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))

NAMES = {0: "O15", 1: "C11", 2: "N13", 3: "C10", 4: "O14"}
REG_COLOURS = {"brain": "#9aa8e0", "skull": "#e8e8e0", "scalp": "#d9b38c",
               "head": "#9aa8e0"}  # uniform-head single region


def classify(X, Y, Z, regions):
    """Region of each world point: the first (lowest-priority) ellipsoid that
    contains it, else 'air'. Ellipsoids are axis-aligned (Euler angles 0)."""
    reg = np.full(len(X), "air", dtype=object)
    assigned = np.zeros(len(X), dtype=bool)
    for _, r in regions.sort_values("priority").iterrows():
        inside = (((X - r.cx_mm) / r.a_mm) ** 2 + ((Y - r.cy_mm) / r.b_mm) ** 2
                  + ((Z - r.cz_mm) / r.c_mm) ** 2) <= 1.0
        take = inside & ~assigned
        reg[take] = r.region
        assigned |= take
    return reg


def _draw_regions(ax, regions, fill):
    """Draw each ellipsoid region in the z-x plane (z horizontal). Outer first so
    inner regions are on top; the brain on top leaves the bone shell visible."""
    for _, r in regions.sort_values("priority", ascending=False).iterrows():
        if r.solid != "ellipsoid":
            continue
        col = REG_COLOURS.get(r.region, "0.7")
        # z-x plane: centre (cz, cx); width 2*c along z, height 2*a along x.
        kw = (dict(facecolor=col, edgecolor="0.3", lw=1.2)
              if fill else dict(fill=False, edgecolor="0.6", lw=1.0))
        ax.add_patch(Ellipse((r.cz_mm, r.cx_mm), 2 * r.c_mm, 2 * r.a_mm,
                             zorder=1, **kw))


def _beam_arrow(ax):
    ax.annotate("", xy=(-74, 0), xytext=(-104, 0),
                arrowprops=dict(arrowstyle="-|>", color="C2", lw=2.2), zorder=5)
    ax.annotate("beam +z", (-104, 0), textcoords="offset points", xytext=(-2, 6),
                color="C2", fontsize=9)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?", default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--out", default=None, help="output PNG (default data/mird_head.png)")
    args = ap.parse_args()

    e = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))
    regions = pd.read_csv(os.path.join(args.data_dir, "phantom_regions.csv"))
    reg = classify(e["anh_x_mm"].to_numpy(), e["anh_y_mm"].to_numpy(),
                   e["anh_z_mm"].to_numpy(), regions)
    reg_names = list(regions.sort_values("priority")["region"])
    isos = [0, 1, 2, 3, 4]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # (A) labelled phantom + beam (filled regions, no emitters).
    ax = axes[0, 0]
    _draw_regions(ax, regions, fill=True)
    _beam_arrow(ax)
    # Leader-line labels to a clear right margin (the regions are ~concentric, so
    # inline labels would overlap). Target a representative point per region.
    label_xy = {"brain": ((92, 25), (18, 6)), "skull": ((92, -18), (0, -78)),
                "scalp": ((92, -58), (6, -90))}
    for _, r in regions.iterrows():
        txt, pt = label_xy.get(r.region, ((r.cz_mm, r.cx_mm), (r.cz_mm, r.cx_mm)))
        ax.annotate(r.region, xy=pt, xytext=txt, fontsize=9, fontweight="bold",
                    va="center", arrowprops=dict(arrowstyle="-", color="0.3", lw=0.8))
    ax.set(xlabel="z [mm] (beam)", ylabel="x [mm] (L-R)",
           title="(A) phantom + lateral beam")
    ax.set_aspect("equal"); ax.set_xlim(-110, 125); ax.set_ylim(-105, 105)

    # (B) emitter trail by isotope over the region outlines.
    ax = axes[0, 1]
    _draw_regions(ax, regions, fill=False)
    for iid in isos:
        m = e["isotope_id"] == iid
        ax.scatter(e.loc[m, "prod_z_mm"], e.loc[m, "prod_x_mm"], s=6, alpha=0.6,
                   linewidths=0, label=NAMES[iid], zorder=3)
    _beam_arrow(ax)
    ax.set(xlabel="z [mm] (beam)", ylabel="x [mm] (L-R)",
           title="(B) emitter trail by isotope")
    ax.set_aspect("equal"); ax.set_xlim(-110, 95); ax.set_ylim(-105, 105)
    ax.legend(loc="lower right", fontsize=8)

    # (C) isotope mix per region -- the heterogeneity payoff.
    ax = axes[1, 0]
    width = 0.15
    for k, iid in enumerate(isos):
        counts = [int(((reg == rn) & (e["isotope_id"] == iid)).sum()) for rn in reg_names]
        ax.bar(np.arange(len(reg_names)) + (k - 2) * width, counts, width, label=NAMES[iid])
    ax.set_xticks(range(len(reg_names))); ax.set_xticklabels(reg_names)
    ax.set(ylabel="emitters", title="(C) isotope mix per region")
    ax.legend(fontsize=8)

    # (D) production profile along the beam.
    ax = axes[1, 1]
    for iid in isos:
        m = e["isotope_id"] == iid
        ax.hist(e.loc[m, "prod_z_mm"], bins=40, range=(-75, 75), histtype="step",
                label=NAMES[iid])
    ax.set(xlabel="production z [mm] (beam)", ylabel="emitters / bin",
           title="(D) production profile along the beam")
    ax.legend(fontsize=8)

    summ = ", ".join(f"{rn} {int((reg==rn).sum())}" for rn in reg_names)
    fig.suptitle(f"MIRD head — {len(e)} emitters ({summ})", fontsize=12)
    fig.tight_layout()
    out = args.out or os.path.join(args.data_dir, "mird_head.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")
    for rn in reg_names:
        m = reg == rn
        o15 = int(((reg == rn) & (e["isotope_id"] == 0)).sum())
        c11 = int(((reg == rn) & (e["isotope_id"] == 1)).sum())
        ratio = f"{o15/c11:.2f}" if c11 else "n/a"
        print(f"  {rn:>5}: {int(m.sum()):5d} emitters   O15/C11 = {ratio}")


if __name__ == "__main__":
    main()
