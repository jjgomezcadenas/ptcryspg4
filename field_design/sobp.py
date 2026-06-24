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

Usage:
    python field_design/sobp.py [--d-prox CM] [--d-dist CM] [--n-layers N]
                                [--rho-rel R] [--out PATH]
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

ALPHA = 0.0022  # cm / MeV^p, Bortfeld water
P = 1.77

_HERE = os.path.dirname(os.path.abspath(__file__))


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
    ap.add_argument("--d-prox", type=float, default=5.5, help="proximal depth [cm]")
    ap.add_argument("--d-dist", type=float, default=10.5, help="distal depth [cm]")
    ap.add_argument("--n-layers", type=int, default=20)
    ap.add_argument("--rho-rel", type=float, default=1.04,
                    help="material density relative to water (brain=1.04)")
    ap.add_argument("--mu", type=float, default=0.025,
                    help="attenuation-correction coefficient [1/cm], tuned to flatten")
    ap.add_argument("--out", default=os.path.join(_HERE, "..", "data", "sobp_layers.csv"))
    args = ap.parse_args()

    R, E, w = design(args.d_prox, args.d_dist, args.n_layers, args.rho_rel, args.mu)

    # --- write the layer table for the Geant4 gun -----------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        f.write("energy_MeV,weight\n")
        for Ei, wi in zip(E, w):
            f.write(f"{Ei:.4f},{wi:.6e}\n")
    print(f"wrote {len(E)} layers -> {args.out}")
    print(f"  depth {args.d_prox}-{args.d_dist} cm  ->  E "
          f"{E[0]:.1f}-{E[-1]:.1f} MeV  (rho_rel={args.rho_rel})")

    # --- sanity plot: weights + analytic SOBP ---------------------------------
    z = np.linspace(0, args.d_dist * 1.15, 800)
    D = analytic_depth_dose(z, R, w, args.mu)
    # Flatness over the plateau, excluding the inherent distal falloff (~last
    # 0.5 cm, where any SOBP must drop to zero).
    falloff_cm = 0.5
    mask = (z >= args.d_prox) & (z <= args.d_dist - falloff_cm)
    plateau = D[mask]
    flatness = (plateau.max() - plateau.min()) / plateau.mean() * 100

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.stem(E, w)
    a1.set(xlabel="layer energy [MeV]", ylabel="weight", title="SOBP layer weights")
    a2.plot(z, D / D[(z >= args.d_prox) & (z <= args.d_dist)].mean(), "k")
    a2.axvspan(args.d_prox, args.d_dist, color="C2", alpha=0.15, label="target")
    a2.set(xlabel="depth [cm]", ylabel="relative dose",
           title=f"analytic SOBP (plateau flatness {flatness:.1f}%)")
    a2.legend()
    fig.tight_layout()
    out_png = os.path.join(os.path.dirname(os.path.abspath(args.out)), "sobp_design.png")
    fig.savefig(out_png, dpi=120)
    print(f"  plateau flatness {flatness:.1f}%  ;  plot -> {out_png}")


if __name__ == "__main__":
    main()
