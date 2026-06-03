#!/usr/bin/env python3
"""Validate Stage-A (proton transport) output as a single dashboard.

Reads the CSVs written to data/ and renders one figure
(data/transport_validation.png) with:
  A positron range per isotope        E emitters-per-proton multiplicity
  B depth-activity profile            F Bragg depth-dose (total/primary/sec)
  C production map (z vs radius)       G activity vs Bragg overlay
  D radial production profile          + per-isotope summary table

The same per-isotope table is also printed to stdout. Uses the system Python.

Usage:
    python analysis_transport/validate_transport.py [data_dir]
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402

ISO_BY_ENDPOINT = sorted(ISOTOPES, key=lambda k: -ISOTOPES[k].endpoint_MeV)


def load(data_dir):
    meta = pd.read_csv(os.path.join(data_dir, "run_meta.csv")).iloc[0]
    emit = pd.read_csv(os.path.join(data_dir, "emitters.csv"))
    prod = emit[["prod_x_mm", "prod_y_mm", "prod_z_mm"]].to_numpy()
    anh = emit[["anh_x_mm", "anh_y_mm", "anh_z_mm"]].to_numpy()
    emit["range_mm"] = np.linalg.norm(anh - prod, axis=1)
    emit["prod_r_mm"] = np.hypot(emit["prod_x_mm"], emit["prod_y_mm"])
    depth = pd.read_csv(os.path.join(data_dir, "depth_dose.csv"))
    return meta, emit, depth


def _name(iid):
    return ISOTOPES[iid].name


def summary_rows(meta, emit):
    rows = []
    for iid in ISO_BY_ENDPOINT:
        iso = ISOTOPES[iid]
        sub = emit[emit.isotope_id == iid]
        n = len(sub)
        rows.append((iso.name, f"{iso.half_life_s:.1f}", f"{iso.endpoint_MeV:.2f}",
                     f"{n}", f"{n/int(meta['n_protons']):.2e}",
                     f"{sub.range_mm.mean():.3f}" if n else "-",
                     f"{sub.range_mm.median():.3f}" if n else "-"))
    return rows


def print_table(meta, emit):
    n = int(meta["n_protons"])
    print(f"\nprotons: {n}   emitters: {len(emit)} ({len(emit)/n:.3e}/p)   "
          f"dose: {meta['dose_total_Gy']:.3e} Gy   "
          f"{meta['physics_list']} G4 {meta['geant4_version']}\n")
    print(f"{'iso':>5} {'T1/2[s]':>8} {'endpt':>6} {'count':>7} {'yield/p':>10} "
          f"{'mean':>7} {'median':>7}")
    print("-" * 56)
    for r in summary_rows(meta, emit):
        print(f"{r[0]:>5} {r[1]:>8} {r[2]:>6} {r[3]:>7} {r[4]:>10} "
              f"{r[5]:>6}m {r[6]:>6}m")


# --- individual panels (each draws on a supplied Axes) ----------------------

def panel_range(ax, emit):  # A
    for iid in ISO_BY_ENDPOINT:
        sub = emit[emit.isotope_id == iid]
        if len(sub) < 2:
            continue
        ax.hist(sub.range_mm, bins=60, histtype="step", linewidth=1.5,
                label=f"{_name(iid)} (n={len(sub)})")
    ax.set(xlabel="positron range |anh-prod| [mm]", ylabel="emitters",
           title="A. Positron range per isotope")
    ax.legend(fontsize=8)


def panel_depth(ax, emit):  # B
    bins = np.linspace(emit.prod_z_mm.min(), emit.prod_z_mm.max(), 80)
    ax.hist(emit.prod_z_mm, bins=bins, histtype="stepfilled", alpha=0.2,
            color="k", label="all")
    for iid in ISO_BY_ENDPOINT:
        sub = emit[emit.isotope_id == iid]
        if len(sub) > 2:
            ax.hist(sub.prod_z_mm, bins=bins, histtype="step", linewidth=1.3,
                    label=_name(iid))
    ax.set(xlabel="production depth z [mm]", ylabel="emitters",
           title="B. Depth-activity profile")
    ax.legend(fontsize=8)


def panel_map(ax, fig, emit):  # C
    # Log color scale: production is concentrated near the entrance/axis, so a
    # linear scale washes everything else to zero. cmin=1 masks empty bins.
    h = ax.hist2d(emit.prod_z_mm, emit.prod_r_mm, bins=[80, 40], cmap="viridis",
                  norm=LogNorm(), cmin=1)
    fig.colorbar(h[3], ax=ax, label="emitters (log)")
    ax.set(xlabel="production depth z [mm]", ylabel="production radius [mm]",
           title="C. Production map")


def panel_radial(ax, emit):  # D
    ax.hist(emit.prod_r_mm, bins=60, histtype="step", linewidth=1.5, color="C3")
    ax.set(xlabel="production radius [mm]", ylabel="emitters",
           title="D. Radial production profile")


def panel_multiplicity(ax, emit, meta):  # E
    n = int(meta["n_protons"])
    per_evt = emit.groupby("event_id").size()
    counts = per_evt.value_counts().sort_index()
    k = [0] + list(counts.index)
    v = [n - len(per_evt)] + list(counts.values)
    ax.bar(k, v, color="C4")
    ax.set(xlabel="emitters per proton", ylabel="protons", yscale="log",
           title="E. Production multiplicity")


def panel_bragg(ax, depth):  # F
    ax.plot(depth.z_mm, depth.edep_total_MeV, label="total", color="k")
    ax.plot(depth.z_mm, depth.edep_primary_MeV, label="primary", color="C0",
            ls="--")
    ax.plot(depth.z_mm, depth.edep_total_MeV - depth.edep_primary_MeV,
            label="secondaries", color="C1", ls=":")
    ax.set(xlabel="depth z [mm]", ylabel="energy deposit [MeV/bin]",
           title="F. Bragg depth-dose")
    ax.legend(fontsize=8)


def panel_overlay(ax, emit, depth):  # G
    ax.plot(depth.z_mm, depth.edep_total_MeV, color="k", label="dose (Bragg)")
    ax.set(xlabel="depth z [mm]", ylabel="energy deposit [MeV/bin]")
    peak_z = depth.z_mm[depth.edep_total_MeV.idxmax()]
    ax.axvline(peak_z, color="k", ls=":", alpha=0.5)
    ax2 = ax.twinx()
    bins = np.linspace(depth.z_mm.min(), depth.z_mm.max(), 80)
    ax2.hist(emit.prod_z_mm, bins=bins, histtype="step", color="C2",
             linewidth=1.8, label="activity")
    ax2.set_ylabel("emitters", color="C2")
    ax.set_title(f"G. Activity vs Bragg (peak z={peak_z:.1f} mm)")
    lines = ax.get_lines()[:1] + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="upper left",
              fontsize=8)


def panel_table(ax, meta, emit):  # summary table
    ax.axis("off")
    cols = ["iso", "T1/2[s]", "endpt", "count", "yield/p", "mean[mm]", "med[mm]"]
    tbl = ax.table(cellText=summary_rows(meta, emit), colLabels=cols,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.4)
    for j in range(len(cols)):  # bold header
        tbl[0, j].set_text_props(weight="bold")
    ax.set_title("Per-isotope summary", fontsize=11)


def main(data_dir):
    meta, emit, depth = load(data_dir)
    print_table(meta, emit)

    fig = plt.figure(figsize=(18, 13), constrained_layout=True)
    gs = fig.add_gridspec(3, 3)
    n = int(meta["n_protons"])
    fig.suptitle(
        f"Stage A validation  -  {n:,} protons, {meta['beam_energy_MeV']:.0f} MeV "
        f"in {meta['phantom_material']}  -  {meta['physics_list']} "
        f"G4 {meta['geant4_version']}  -  dose {meta['dose_total_Gy']:.2e} Gy",
        fontsize=14, weight="bold")

    panel_range(fig.add_subplot(gs[0, 0]), emit)
    panel_bragg(fig.add_subplot(gs[0, 1]), depth)
    panel_overlay(fig.add_subplot(gs[0, 2]), emit, depth)
    panel_depth(fig.add_subplot(gs[1, 0]), emit)
    panel_map(fig.add_subplot(gs[1, 1]), fig, emit)
    panel_radial(fig.add_subplot(gs[1, 2]), emit)
    panel_multiplicity(fig.add_subplot(gs[2, 0]), emit, meta)
    panel_table(fig.add_subplot(gs[2, 1:3]), meta, emit)

    out = os.path.join(data_dir, "transport_validation.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"\nsaved dashboard -> {out}")


if __name__ == "__main__":
    ddir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "..", "data")
    main(ddir)
