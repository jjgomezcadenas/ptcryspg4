#!/usr/bin/env python3
"""Visualize the MIRD-head run: phantom, beam, and emitter paths (Phase 1).

Renders, from a head run's emitters.csv:
  (A) the labelled phantom + beam: the three regions (scalp / skull / brain) in
      the beam-axis plane, filled, with the lateral pencil;
  (B) the emitter trail: production points coloured by isotope over the region
      outlines -- the activity the proton leaves through skull → brain;
  (C) the isotope mix per region (brain vs skull vs scalp) -- the payoff of the
      heterogeneous head: bone shifts the local O15/C11 balance;
  (D) the production profile along the beam.

Head geometry is mirrored from StageAConfig.hh (Phase 2 will write it into
run_meta.csv). Local frame: x=L-R, y=A-P, z=S-I, origin at the skull/scalp
centre; the brain sits +BRAIN_OFFSET toward the crown. The head is placed so the
brain centre is the world origin and the L-R axis is the beam (+z); world =
(z_loc - BRAIN_OFFSET, y_loc, -x_loc), so x_loc=-Z, y_loc=Y, z_loc=X+BRAIN_OFFSET.

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

# --- head geometry, mirrored from StageAConfig.hh (semi-axes / centres, mm) ---
BRAIN_OFFSET = 10.0  # brain centre vs skull/scalp centre (head-local z)
SCALP = dict(c=(0., 0., 0.), s=(72., 102., 87.))
CRAN_OUT = dict(c=(0., 0., 0.), s=(68., 98., 83.))
CRAN_IN = dict(c=(0., 0., BRAIN_OFFSET), s=(60., 90., 65.))  # = brain cavity
BRAIN = dict(c=(0., 0., BRAIN_OFFSET), s=(60., 90., 65.))

NAMES = {0: "O15", 1: "C11", 2: "N13", 3: "C10", 4: "O14"}
REG_COLOURS = {"scalp": "#d9b38c", "skull": "#e8e8e0", "brain": "#9aa8e0"}


def to_local(X, Y, Z):
    """World coords [mm] -> head-local coords [mm] (inverse of the placement)."""
    return -Z, Y, X + BRAIN_OFFSET


def _inside(xl, yl, zl, ell):
    cx, cy, cz = ell["c"]
    sx, sy, sz = ell["s"]
    return ((xl - cx) / sx) ** 2 + ((yl - cy) / sy) ** 2 + ((zl - cz) / sz) ** 2 <= 1.0


def classify(X, Y, Z):
    """Region of each world point: 'brain' | 'skull' | 'scalp'."""
    xl, yl, zl = to_local(X, Y, Z)
    in_brain = _inside(xl, yl, zl, BRAIN)
    in_shell = _inside(xl, yl, zl, CRAN_OUT) & ~_inside(xl, yl, zl, CRAN_IN)
    reg = np.full(len(X), "scalp", dtype=object)
    reg[in_shell] = "skull"
    reg[in_brain] = "brain"
    return reg


def _ellipse_zx(ell):
    """matplotlib Ellipse for a local ellipsoid in the z-x plane (z horizontal)."""
    cx, cy, cz = ell["c"]
    sx, sy, sz = ell["s"]
    # world centre (Z=-cx, X=cz-OFFSET); world semis: Z<-sx (local x), X<-sz (local z)
    return (-cx, cz - BRAIN_OFFSET), 2 * sx, 2 * sz


def _draw_regions(ax, fill):
    """Draw scalp/skull/brain in the z-x plane; filled (schematic) or outlines."""
    for ell, key in [(SCALP, "scalp"), (CRAN_OUT, "skull"), (BRAIN, "brain")]:
        (z_c, x_c), w, h = _ellipse_zx(ell)
        if fill:
            ax.add_patch(Ellipse((z_c, x_c), w, h, facecolor=REG_COLOURS[key],
                                 edgecolor="0.3", lw=1.2, zorder=1))
        else:
            ax.add_patch(Ellipse((z_c, x_c), w, h, fill=False, edgecolor="0.6",
                                 lw=1.0, zorder=1))


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
    reg = classify(e["anh_x_mm"].to_numpy(), e["anh_y_mm"].to_numpy(),
                   e["anh_z_mm"].to_numpy())
    regions, isos = ["brain", "skull", "scalp"], [0, 1, 2, 3, 4]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # (A) labelled phantom + beam (filled regions, no emitters).
    ax = axes[0, 0]
    _draw_regions(ax, fill=True)
    _beam_arrow(ax)
    ax.text(0, BRAIN_OFFSET, "brain", ha="center", va="center", fontweight="bold")
    ax.annotate("skull (bone)", xy=(0, -78), xytext=(48, -92),
                arrowprops=dict(arrowstyle="-", color="0.3"), fontsize=9)
    ax.annotate("scalp", xy=(0, 76), xytext=(48, 92),
                arrowprops=dict(arrowstyle="-", color="0.3"), fontsize=9)
    ax.set(xlabel="z [mm] (beam)", ylabel="x [mm] (L-R)",
           title="(A) phantom + lateral beam")
    ax.set_aspect("equal"); ax.set_xlim(-110, 95); ax.set_ylim(-105, 105)

    # (B) emitter trail by isotope over the region outlines.
    ax = axes[0, 1]
    _draw_regions(ax, fill=False)
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
        counts = [int(((reg == r) & (e["isotope_id"] == iid)).sum()) for r in regions]
        ax.bar(np.arange(len(regions)) + (k - 2) * width, counts, width, label=NAMES[iid])
    ax.set_xticks(range(len(regions))); ax.set_xticklabels(regions)
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

    nb, ns, nc = (int((reg == r).sum()) for r in regions)
    fig.suptitle(f"MIRD head — {len(e)} emitters (brain {nb}, skull {ns}, scalp {nc})",
                 fontsize=12)
    fig.tight_layout()
    out = args.out or os.path.join(args.data_dir, "mird_head.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")
    for r in regions:
        m = reg == r
        o15 = int(((reg == r) & (e["isotope_id"] == 0)).sum())
        c11 = int(((reg == r) & (e["isotope_id"] == 1)).sum())
        ratio = f"{o15/c11:.2f}" if c11 else "n/a"
        print(f"  {r:>5}: {int(m.sum()):5d} emitters   O15/C11 = {ratio}")


if __name__ == "__main__":
    main()
