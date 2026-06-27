#!/usr/bin/env python3
"""Design a Spread-Out Bragg Peak (SOBP): proton energy-layer weights.

Implements latex/02_beam_design.tex. Pristine-peak ranges span the target depth; energies
come from the Bortfeld power law R = alpha*E^p (water), scaled to the material by
density; flattening weights are the bin-integrated

    w_i ∝ (Rd - Ri + Δ/2)^(1/p) - (Rd - Ri - Δ/2)^(1/p)

of the continuous fluence Phi(R) ∝ (Rd - R)^(1/p - 1). Writes
data/sobp_layers.csv (energy_MeV, weight) for the Geant4 gun and a sanity plot
of the analytic SOBP (idealized peaks + range straggling; the G4 depth-dose is
the real check).

The field is phantom-specific (the medium enters via WEPL, the target window
too), so each design is written per field as data/field/<label>_sobp_layers.csv
(+ a provenance _meta.csv) and a run loads the table for its own geometry.

Two modes:
  * homogeneous (default): the cylinder/brain field, range scaled by --rho-rel.
  * WEPL (--from-run RUN_DIR): design through the run's heterogeneous phantom by
    ray-tracing water-equivalent path length (Sec. 5 of the doc), so the bone
    offset of the MIRD head is absorbed by construction.

Usage:
    python field_design/sobp.py [--d-prox CM] [--d-dist CM] [--n-layers N]
                                [--rho-rel R] [--label NAME]
    python field_design/sobp.py --from-run data/runs/mird_head_pencil_1e5
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

ALPHA = 0.0022  # cm / MeV^p, Bortfeld water
P = 1.77

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from phantom_material import MATERIALS, relative_stopping_power  # noqa: E402
from regions import material_at  # noqa: E402


def rsp_by_material(regions, energy_MeV=150.0):
    """RSP for each distinct material in the regions (air -> 0)."""
    return {m: (relative_stopping_power(MATERIALS[m], energy_MeV) if m in MATERIALS
                else 0.0)
            for m in regions["material"].unique()}


def wepl_curve(regions, z0_mm, d_max_mm, rsp_by, dz=0.05):
    """(depths, WEPL) [mm] along the central axis (x=y=0) from entrance z0:
    geometric depth reweighted by the RSP of the medium there (Eq. WEPL)."""
    depths = np.arange(0.0, d_max_mm + dz, dz)
    rsp = np.array([rsp_by.get(material_at(regions, 0.0, 0.0, z0_mm + d), 0.0)
                    for d in depths])
    wepl = np.concatenate([[0.0], np.cumsum(0.5 * (rsp[1:] + rsp[:-1]) * dz)])
    return depths, wepl


def energy_for_range(R_cm, rho_rel):
    """Proton energy [MeV] to reach range R_cm in a material of rel. density."""
    R_water = R_cm * rho_rel        # water-equivalent depth
    return (R_water / ALPHA) ** (1.0 / P)


def design(d_prox, d_dist, n_layers, rho_rel, mu):
    """Return per-layer ranges [cm], energies [MeV] and normalized weights."""
    R = np.linspace(d_prox, d_dist, n_layers)   # peak ranges, proximal->distal
    Rd = d_dist
    delta = (d_dist - d_prox) / (n_layers - 1)
    # Flat plateau requires fluence Phi(R) ∝ (Rd-R)^(-1/p) (Abel inversion of
    # D(z)=const); integrated over each range bin -> exponent (1 - 1/p). Clip at
    # the distal edge (no peaks beyond Rd) so the deepest bin stays finite.
    hi = np.clip(Rd - R + delta / 2, 0.0, None)
    lo = np.clip(Rd - R - delta / 2, 0.0, None)
    e = 1.0 - 1.0 / P
    w = hi ** e - lo ** e
    # Compensate the proximal->distal droop of the *delivered* field: protons
    # are lost (mainly to nuclear reactions) on the way in, so deeper layers
    # arrive depleted. Boost each layer by exp(mu*R) (mu [1/cm] tuned against the
    # simulated plateau).
    w *= np.exp(mu * R)
    w /= w.sum()
    E = energy_for_range(R, rho_rel)
    return R, E, w


def analytic_depth_dose(z, R, w, mu):
    """Idealized SOBP depth-dose on grid z [cm]: peaks + attenuation, smeared."""
    D = np.zeros_like(z)
    for Ri, wi in zip(R, w):
        m = z < Ri
        D[m] += wi * (Ri - z[m]) ** (1.0 / P - 1.0)
    D *= np.exp(-mu * z)   # beam attenuation reaching depth z (same model as above)
    # Range straggling ~1.2% of the distal range -> smooth the cusps.
    dz = z[1] - z[0]
    sigma = 0.012 * R[-1] / dz
    return gaussian_filter1d(D, sigma)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-run", default=None,
                    help="run dir: design in WEPL through its heterogeneous phantom")
    ap.add_argument("--d-prox", type=float, default=None,
                    help="proximal target depth [cm] (default: run's, or 5.5)")
    ap.add_argument("--d-dist", type=float, default=None,
                    help="distal target depth [cm] (default: run's, or 10.5)")
    ap.add_argument("--n-layers", type=int, default=20)
    ap.add_argument("--rho-rel", type=float, default=1.04,
                    help="homogeneous mode: density relative to water (brain=1.04)")
    ap.add_argument("--mu", type=float, default=0.025,
                    help="attenuation-correction coefficient [1/cm], tuned to flatten")
    ap.add_argument("--rsp-energy", type=float, default=150.0,
                    help="proton energy [MeV] at which RSP is evaluated")
    ap.add_argument("--label", default=None,
                    help="field label (default: the geometry, or 'cylinder')")
    ap.add_argument("--field-dir", default=os.path.join(_HERE, "..", "data", "field"))
    args = ap.parse_args()

    if args.from_run:  # ---- WEPL mode: design through the run's phantom --------
        regions = pd.read_csv(os.path.join(args.from_run, "phantom_regions.csv"))
        rmeta = pd.read_csv(os.path.join(args.from_run, "run_meta.csv")).iloc[0]
        geometry = str(rmeta["geometry"])
        d_prox = args.d_prox if args.d_prox is not None else float(rmeta["target_prox_depth_mm"]) / 10
        d_dist = args.d_dist if args.d_dist is not None else float(rmeta["target_dist_depth_mm"]) / 10
        z0 = -0.5 * float(rmeta["phantom_length_mm"])  # entrance face [mm]
        rsp_by = rsp_by_material(regions, args.rsp_energy)
        depths, wepl = wepl_curve(regions, z0, d_dist * 10.0, rsp_by)
        wp = float(np.interp(d_prox * 10.0, depths, wepl)) / 10.0  # WEPL window [cm]
        wd = float(np.interp(d_dist * 10.0, depths, wepl)) / 10.0
        R, E, w = design(wp, wd, args.n_layers, 1.0, args.mu)
        label = args.label or geometry
        design_lo, design_hi, x_is_wepl = wp, wd, True
        meta = dict(design_geometry=geometry, mode="wepl",
                    materials="|".join(dict.fromkeys(regions["material"])),
                    rsp_energy_MeV=args.rsp_energy,
                    wepl_prox_mm=wp * 10.0, wepl_dist_mm=wd * 10.0)
        summary = (f"{geometry}: target {d_prox}-{d_dist} cm geom -> "
                   f"{wp:.2f}-{wd:.2f} cm WEPL")
    else:  # ---- homogeneous mode (cylinder/brain), range scaled by rho_rel ----
        d_prox = args.d_prox if args.d_prox is not None else 5.5
        d_dist = args.d_dist if args.d_dist is not None else 10.5
        R, E, w = design(d_prox, d_dist, args.n_layers, args.rho_rel, args.mu)
        label = args.label or "cylinder"
        geometry = "cylinder"
        design_lo, design_hi, x_is_wepl = d_prox, d_dist, False
        meta = dict(design_geometry=geometry, mode="homogeneous", materials="",
                    rsp_energy_MeV=args.rsp_energy,
                    wepl_prox_mm=d_prox * 10.0 * args.rho_rel,
                    wepl_dist_mm=d_dist * 10.0 * args.rho_rel)
        summary = f"{geometry}: depth {d_prox}-{d_dist} cm (rho_rel={args.rho_rel})"

    # Common provenance (what the field was designed for -> check_run guard).
    meta.update(d_prox_mm=d_prox * 10.0, d_dist_mm=d_dist * 10.0, rho_rel=
                (1.0 if args.from_run else args.rho_rel), mu=args.mu,
                n_layers=args.n_layers, e_min_MeV=float(E.min()), e_max_MeV=float(E.max()))

    # --- write the per-field layer table + provenance meta --------------------
    os.makedirs(args.field_dir, exist_ok=True)
    out = os.path.join(args.field_dir, f"{label}_sobp_layers.csv")
    with open(out, "w") as f:
        f.write("energy_MeV,weight\n")
        for Ei, wi in zip(E, w):
            f.write(f"{Ei:.4f},{wi:.6e}\n")
    meta_path = os.path.join(args.field_dir, f"{label}_sobp_layers_meta.csv")
    pd.DataFrame([meta]).to_csv(meta_path, index=False)
    print(f"wrote {len(E)} layers -> {out}")
    print(f"  {summary}  ->  E {E[0]:.1f}-{E[-1]:.1f} MeV")

    # --- sanity plot: weights + analytic SOBP (in the design coordinate) -------
    z = np.linspace(0, design_hi * 1.15, 800)
    D = analytic_depth_dose(z, R, w, args.mu)
    falloff_cm = 0.5  # exclude the inherent distal falloff from the flatness
    plateau = D[(z >= design_lo) & (z <= design_hi - falloff_cm)]
    flatness = (plateau.max() - plateau.min()) / plateau.mean() * 100
    xlabel = "water-equivalent depth [cm]" if x_is_wepl else "depth [cm]"

    n_panels = 3 if x_is_wepl else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 4))
    a1, a2 = axes[0], axes[1]
    a1.stem(E, w)
    a1.set(xlabel="layer energy [MeV]", ylabel="weight",
           title=f"{label} SOBP layer weights")
    a2.plot(z, D / D[(z >= design_lo) & (z <= design_hi)].mean(), "k")
    a2.axvspan(design_lo, design_hi, color="C2", alpha=0.15, label="target")
    a2.set(xlabel=xlabel, ylabel="relative dose",
           title=f"analytic SOBP (plateau flatness {flatness:.1f}%)")
    a2.legend()
    if x_is_wepl:  # WEPL ray-trace: the bone kink that makes the field phantom-specific
        a3 = axes[2]
        dcm, wcm = depths / 10.0, wepl / 10.0
        a3.plot(dcm, wcm, "C0", label="WEPL(depth)")
        a3.plot([0, dcm.max()], [0, dcm.max()], "k--", lw=0.8, label="water (RSP=1)")
        a3.axvspan(d_prox, d_dist, color="C2", alpha=0.15, label="target")
        a3.set(xlabel="geometric depth [cm]", ylabel="WEPL [cm]",
               title="WEPL ray-trace through the phantom")
        a3.legend(fontsize=8)
    fig.tight_layout()
    out_png = os.path.join(args.field_dir, f"{label}_sobp_design.png")
    fig.savefig(out_png, dpi=120)
    print(f"  plateau flatness {flatness:.1f}%  ;  meta -> {meta_path}")
    print(f"  plot -> {out_png}")


if __name__ == "__main__":
    main()
