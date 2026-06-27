#!/usr/bin/env python3
"""Draw the phantom medium + target box + beam (two orthogonal cross-sections).

From phantom_regions.csv + run_meta.csv: the medium regions as drawn solids (the
nested ellipses of the head, or the cylinder), the dose-normalization target box,
and the beam, in the z-x (beam + L/R) and z-y (beam + A/P) planes. Geometry-
agnostic. A target box poking outside the medium is then obvious by eye — the
visual companion to check_run.py. Writes <run_dir>/figures/phantom.png.

(Distinct from plot_geometry.py, which is the cylinder scenario-frame figure with
the source cloud, used in latex/04_source_reference.tex.)

Usage:
    python analysis_transport/plot_phantom.py <run_dir> [--out PNG]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse, Rectangle  # noqa: E402
import pandas as pd  # noqa: E402

REG_COLOURS = {
    "phantom": "#9db8d2", "head": "#9db8d2",
    "brain": "#7e9ed6", "skull": "#e8e4cf", "scalp": "#d8b89a",
}
_FALLBACK = ["#9db8d2", "#e8e4cf", "#d8b89a", "#a8c6a0", "#d2a0a0"]


def region_patch(r, centre_col, semi_col):
    """A patch for region row r in the (z, transverse) plane named by the columns."""
    colour = REG_COLOURS.get(r.region, _FALLBACK[r.priority % len(_FALLBACK)])
    cv = getattr(r, centre_col)       # transverse centre (cx or cy)
    semi_v = getattr(r, semi_col)     # transverse semi-axis / radius (a or b)
    if r.solid == "cylinder":         # (a,b,c) = (radius, radius, half-length)
        return Rectangle((r.cz_mm - r.c_mm, cv - semi_v), 2 * r.c_mm, 2 * semi_v,
                         facecolor=colour, edgecolor="k", alpha=0.55, lw=1.2,
                         label=r.region)
    return Ellipse((r.cz_mm, cv), 2 * r.c_mm, 2 * semi_v, facecolor=colour,
                   edgecolor="k", alpha=0.55, lw=1.2, label=r.region)


def draw_plane(ax, regions, meta, centre_col, semi_col, vlabel):
    """One cross-section: beam z (horizontal) vs a transverse axis (vertical)."""
    L = float(meta["phantom_length_mm"])
    entrance = -0.5 * L
    # Outer regions first (high priority) so inner ones draw on top.
    for r in regions.sort_values("priority", ascending=False).itertuples():
        ax.add_patch(region_patch(r, centre_col, semi_col))

    # Target box: z in [entrance+prox, entrance+dist], +-target_radius transverse.
    prox = float(meta["target_prox_depth_mm"])
    dist = float(meta["target_dist_depth_mm"])
    rad = float(meta["target_radius_mm"])
    ax.add_patch(Rectangle((entrance + prox, -rad), dist - prox, 2 * rad,
                           fill=False, edgecolor="crimson", lw=2.0, ls="--",
                           label="target box"))

    # Beam: enters the entrance face along +z, on the axis.
    ax.annotate("", xy=(entrance + 0.12 * L, 0), xytext=(entrance - 0.10 * L, 0),
                arrowprops=dict(arrowstyle="-|>", color="orangered", lw=2.2))
    ax.axhline(0, color="orangered", lw=0.8, ls=":", alpha=0.5)
    ax.axvline(entrance, color="gray", lw=0.8, ls=":", alpha=0.7)

    half = 0.62 * L
    ax.set(xlim=(entrance - 0.12 * L, -entrance + 0.05 * L), ylim=(-half, half),
           xlabel="z (beam) [mm]", ylabel=vlabel)
    ax.set_aspect("equal")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a Stage-A run directory")
    ap.add_argument("--out", default=None,
                    help="output PNG (default <run_dir>/figures/phantom.png)")
    args = ap.parse_args()

    meta = pd.read_csv(os.path.join(args.run_dir, "run_meta.csv")).iloc[0]
    regions = pd.read_csv(os.path.join(args.run_dir, "phantom_regions.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
    draw_plane(axes[0], regions, meta, "cx_mm", "a_mm", "x [mm]")
    draw_plane(axes[1], regions, meta, "cy_mm", "b_mm", "y [mm]")
    axes[0].set_title("z-x plane (beam + L/R)")
    axes[1].set_title("z-y plane (beam + A/P)")
    handles, labels = axes[0].get_legend_handles_labels()
    seen = dict(zip(labels, handles))  # dedupe
    axes[0].legend(seen.values(), seen.keys(), loc="upper right", fontsize=8)

    rad = float(meta["target_radius_mm"])
    prox, dist = float(meta["target_prox_depth_mm"]), float(meta["target_dist_depth_mm"])
    fig.suptitle(
        f"Phantom: {meta['geometry']} ({meta['phantom_material']})  -  "
        f"target box {2*rad/10:g} cm dia x {(dist-prox)/10:g} cm "
        f"at {prox/10:g}-{dist/10:g} cm depth",
        fontsize=12, weight="bold")
    fig.tight_layout()

    if args.out:
        out = args.out
    else:
        figdir = os.path.join(args.run_dir, "figures")
        os.makedirs(figdir, exist_ok=True)
        out = os.path.join(figdir, "phantom.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
