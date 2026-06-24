#!/usr/bin/env python3
"""Visualize the MIRD-head run: phantom, beam, and emitter paths (Phase 1).

Renders, from a head run's emitters.csv:
  (A) the head in the beam-axis plane (z-x): the three regions (scalp / skull /
      brain) as cross-section outlines, the lateral beam, and the annihilation
      points coloured by region;
  (B) the same in the z-y plane;
  (C) the isotope mix per region (brain vs skull vs scalp) -- the payoff of the
      heterogeneous head: bone shifts the local O15/C11 balance;
  (D) the production profile along the beam (where activity is created).

The head geometry is mirrored from StageAConfig.hh (Phase 1 does not yet write it
into run_meta.csv -- that is Phase 2). Geometry frame: the head is built in a
local frame (x=L-R, y=A-P, z=S-I) and placed so the brain centre is at the world
origin with the L-R axis along the beam (+z); world = (z_loc - BRAIN_POS_Z,
y_loc, -x_loc), so the inverse is x_loc=-Z, y_loc=Y, z_loc=X+BRAIN_POS_Z.

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
# Head-local frame (x=L-R, y=A-P, z=S-I). Phase 2 will write these to run_meta.
BRAIN_POS_Z = 87.5
SKULL_POS_Z = 77.5
BRAIN = dict(c=(0., 0., BRAIN_POS_Z), s=(60., 90., 65.))
CRAN_OUT = dict(c=(0., 0., SKULL_POS_Z), s=(68., 98., 83.))
CRAN_IN = dict(c=(0., 0., SKULL_POS_Z + 10.), s=(60., 90., 65.))  # = brain cavity
HEAD_TUBE = dict(c=(0., 0., 0.), s=(70., 100., 77.5))       # elliptical tube
HEAD_CAP = dict(c=(0., 0., 77.5), s=(70., 100., 85.))        # ellipsoid cap


def to_local(X, Y, Z):
    """World coords [mm] -> head-local coords [mm] (inverse of the placement)."""
    return -Z, Y, X + BRAIN_POS_Z


def _inside(xl, yl, zl, ell):
    cx, cy, cz = ell["c"]
    ax, by, cz_ = ell["s"]
    return ((xl - cx) / ax) ** 2 + ((yl - cy) / by) ** 2 + ((zl - cz) / cz_) ** 2 <= 1.0


def classify(X, Y, Z):
    """Region of each world point: 'brain' | 'skull' | 'scalp'."""
    xl, yl, zl = to_local(X, Y, Z)
    in_brain = _inside(xl, yl, zl, BRAIN)
    in_shell = _inside(xl, yl, zl, CRAN_OUT) & ~_inside(xl, yl, zl, CRAN_IN)
    reg = np.full(len(X), "scalp", dtype=object)
    reg[in_shell] = "skull"
    reg[in_brain] = "brain"
    return reg


def _world_ellipse_zx(ell):
    """(centre_zx, width_z, height_x) for a local ellipsoid drawn in the z-x plane."""
    cx, cy, cz = ell["c"]
    ax, by, cz_ = ell["s"]
    # world centre = (cz - BRAIN_POS_Z, cy, -cx); world semis: x<-cz_, y<-by, z<-ax
    return (-cx, cz - BRAIN_POS_Z), 2 * ax, 2 * cz_  # (x_c, z_c), w=2*sz=2ax, h=2*sx=2cz_


def _world_ellipse_zy(ell):
    cx, cy, cz = ell["c"]
    ax, by, cz_ = ell["s"]
    return (cy, cz - BRAIN_POS_Z), 2 * ax, 2 * by  # (y_c, z_c), w=2*sz, h=2*sy


REG_COLOURS = {"brain": "C0", "skull": "C1", "scalp": "C7"}


def _draw_regions(ax, plane):
    """Outline head/skull/brain in the chosen plane ('zx' or 'zy'); z horizontal."""
    fn = _world_ellipse_zx if plane == "zx" else _world_ellipse_zy
    specs = [(HEAD_TUBE, "0.5", "head"), (HEAD_CAP, "0.5", None),
             (CRAN_OUT, "C1", "skull"), (BRAIN, "C0", "brain")]
    for ell, col, lab in specs:
        (vert_c, z_c), w, h = fn(ell)  # w spans z (horizontal), h spans vert
        ax.add_patch(Ellipse((z_c, vert_c), w, h, fill=False, edgecolor=col,
                             lw=1.6, label=lab))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", nargs="?", default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--out", default=None, help="output PNG (default data/mird_head.png)")
    args = ap.parse_args()

    e = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))
    names = {0: "O15", 1: "C11", 2: "N13", 3: "C10", 4: "O14"}
    reg = classify(e["anh_x_mm"].to_numpy(), e["anh_y_mm"].to_numpy(),
                   e["anh_z_mm"].to_numpy())

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # (A) beam-axis plane z-x, annihilation points coloured by region.
    ax = axes[0, 0]
    _draw_regions(ax, "zx")
    for r, col in REG_COLOURS.items():
        m = reg == r
        ax.scatter(e.loc[m, "anh_z_mm"], e.loc[m, "anh_x_mm"], s=4, alpha=0.4,
                   c=col, linewidths=0)
    ax.annotate("", xy=(-70, 0), xytext=(-110, 0),
                arrowprops=dict(arrowstyle="-|>", color="C2", lw=2))
    ax.annotate("beam +z", (-110, 0), textcoords="offset points", xytext=(0, 6), color="C2")
    ax.set(xlabel="z [mm] (beam)", ylabel="x [mm] (L-R)",
           title="(A) beam-axis plane: regions + annihilations")
    ax.set_aspect("equal"); ax.legend(loc="lower right", fontsize=8)

    # (B) z-y plane.
    ax = axes[0, 1]
    _draw_regions(ax, "zy")
    for r, col in REG_COLOURS.items():
        m = reg == r
        ax.scatter(e.loc[m, "anh_z_mm"], e.loc[m, "anh_y_mm"], s=4, alpha=0.4,
                   c=col, linewidths=0)
    ax.set(xlabel="z [mm] (beam)", ylabel="y [mm] (A-P)",
           title="(B) transverse plane: regions + annihilations")
    ax.set_aspect("equal")

    # (C) isotope mix per region -- the heterogeneity payoff.
    ax = axes[1, 0]
    regions = ["brain", "skull", "scalp"]
    isos = [0, 1, 2, 3, 4]
    width = 0.15
    for k, iid in enumerate(isos):
        counts = [int(((reg == r) & (e["isotope_id"] == iid)).sum()) for r in regions]
        ax.bar(np.arange(len(regions)) + (k - 2) * width, counts, width,
               label=names[iid])
    ax.set_xticks(range(len(regions))); ax.set_xticklabels(regions)
    ax.set(ylabel="emitters", title="(C) isotope mix per region")
    ax.legend(fontsize=8)

    # (D) production along the beam (where activity is created).
    ax = axes[1, 1]
    for iid in isos:
        m = e["isotope_id"] == iid
        ax.hist(e.loc[m, "prod_z_mm"], bins=40, range=(-70, 70), histtype="step",
                label=names[iid])
    ax.set(xlabel="production z [mm] (beam)", ylabel="emitters / bin",
           title="(D) production profile along the beam")
    ax.legend(fontsize=8)

    fig.suptitle(f"MIRD head — {len(e)} emitters "
                 f"(brain {int((reg=='brain').sum())}, skull {int((reg=='skull').sum())}, "
                 f"scalp {int((reg=='scalp').sum())})", fontsize=12)
    fig.tight_layout()
    out = args.out or os.path.join(args.data_dir, "mird_head.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")
    # Quick numeric check to stdout.
    for r in regions:
        m = reg == r
        n = int(m.sum())
        o15 = int(((reg == r) & (e["isotope_id"] == 0)).sum())
        c11 = int(((reg == r) & (e["isotope_id"] == 1)).sum())
        ratio = f"{o15/c11:.2f}" if c11 else "n/a"
        print(f"  {r:>5}: {n:5d} emitters   O15/C11 = {ratio}")


if __name__ == "__main__":
    main()
