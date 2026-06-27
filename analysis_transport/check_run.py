#!/usr/bin/env python3
"""Independent sanity checks on a Stage-A run (catches geometry/normalization bugs).

Recomputes, in plain Python from run_meta.csv + phantom_regions.csv, what the
Geant4 RunAction claimed, and fails (exit 1) on any violation. Geometry-agnostic
(cylinder + head); needs no Geant4. Run it on any run dir — or a tiny smoke run —
as the gate before trusting a dose normalization.

Each check is independent; several would have caught the head depth-reference bug
(target box left in air, target mass using the wrong density):

  - target box centre sits inside a medium region (not air)
  - target_mass_g == pi r^2 L * rho(region at box centre)   [right density]
  - target_dose_Gy >= dose_total_Gy                         [focused box, not air]
  - Np_per_Gy == n_protons / target_dose_Gy                 [internal consistency]
  - target depths within the phantom; box radius within the transverse extent
  - regions inside the phantom bounding box; head carve order brain<skull<scalp

Usage:
    python analysis_transport/check_run.py <run_dir>
"""

import argparse
import math
import os
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from phantom_material import MATERIALS  # noqa: E402

# These are exact recomputes (geometry + density, no MC noise), so the tolerance
# is tight — loose enough only for float rounding, tight enough to catch a
# wrong-material density swap (e.g. brain 1.04 vs scalp 1.03, ~1% apart).
RTOL = 0.005


def contains(r, x, y, z):
    """True if world point (x,y,z) is inside region row r (mm)."""
    dx, dy, dz = x - r.cx_mm, y - r.cy_mm, z - r.cz_mm
    if r.solid == "cylinder":  # (a,b,c) = (radius, radius, half-length)
        return dx * dx + dy * dy <= r.a_mm ** 2 and abs(dz) <= r.c_mm
    return (dx / r.a_mm) ** 2 + (dy / r.b_mm) ** 2 + (dz / r.c_mm) ** 2 <= 1.0


def material_at(regions, x, y, z):
    """NIST material of the first priority-ordered region containing the point."""
    for r in regions.sort_values("priority").itertuples():
        if contains(r, x, y, z):
            return r.material
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a Stage-A run directory")
    args = ap.parse_args()

    meta = pd.read_csv(os.path.join(args.run_dir, "run_meta.csv")).iloc[0]
    regions = pd.read_csv(os.path.join(args.run_dir, "phantom_regions.csv"))

    L = float(meta["phantom_length_mm"])
    diam = float(meta["phantom_diameter_mm"])
    entrance = -0.5 * L  # beam enters at z = -L/2
    prox = float(meta["target_prox_depth_mm"])
    dist = float(meta["target_dist_depth_mm"])
    rad = float(meta["target_radius_mm"])
    z_c = entrance + 0.5 * (prox + dist)  # target-box centre along the beam
    n_p = float(meta["n_protons"])
    t_dose = float(meta["target_dose_Gy"])
    t_mass = float(meta["target_mass_g"])
    tot_dose = float(meta["dose_total_Gy"])
    np_gy = float(meta["Np_per_Gy"])
    geom = str(meta["geometry"])

    checks = []  # (name, ok, detail)

    def add(name, ok, detail):
        checks.append((name, bool(ok), detail))

    # 1. Target box centre is in a medium, not air.
    mat = material_at(regions, 0.0, 0.0, z_c)
    add("target box centre in medium", mat is not None,
        f"centre z={z_c:.1f} mm -> {mat or 'AIR'}")

    # 2. Reported target mass matches pi r^2 L * rho of the medium at the centre.
    if mat is not None and mat in MATERIALS:
        vol_cm3 = math.pi * rad ** 2 * (dist - prox) / 1000.0  # mm^3 -> cm^3
        exp_mass = vol_cm3 * MATERIALS[mat].density_g_cm3
        ok = abs(t_mass - exp_mass) <= RTOL * exp_mass
        add("target mass matches medium density", ok,
            f"reported {t_mass:.2f} g vs {exp_mass:.2f} g "
            f"(rho_{mat}={MATERIALS[mat].density_g_cm3:g})")
    else:
        add("target mass matches medium density", False,
            f"no density for material '{mat}' (not in registry)")

    # 3. The focused target box must see at least the whole-phantom average dose
    #    (it sits in the Bragg region); a box in air gives orders less.
    add("target dose >= whole-phantom dose", t_dose >= tot_dose,
        f"target {t_dose:.3e} Gy vs whole-phantom {tot_dose:.3e} Gy")

    # 4. Np(1 Gy) is internally consistent with n_protons / target_dose.
    exp_np = n_p / t_dose if t_dose > 0 else float("inf")
    ok = math.isfinite(exp_np) and abs(np_gy - exp_np) <= RTOL * exp_np
    add("Np_per_Gy consistent", ok, f"reported {np_gy:.3e} vs {exp_np:.3e}")

    # 5. Target box fits inside the phantom (depths and radius).
    add("target depths within phantom", 0.0 <= prox < dist <= L,
        f"[{prox:g}, {dist:g}] mm in (0, {L:g})")
    add("target radius within phantom", rad <= 0.5 * diam,
        f"r={rad:g} mm vs phantom half-width {0.5*diam:g} mm")

    # 6. Every region sits within the reported bounding box.
    bbox_ok, worst = True, ""
    for r in regions.itertuples():
        hz, hx, hy = abs(r.cz_mm) + r.c_mm, abs(r.cx_mm) + r.a_mm, abs(r.cy_mm) + r.b_mm
        if hz > 0.5 * L + 1e-6 or max(hx, hy) > 0.5 * diam + 1e-6:
            bbox_ok = False
            worst = r.region
    add("regions within bounding box", bbox_ok,
        "ok" if bbox_ok else f"region '{worst}' exceeds bbox")

    # 6b. Head carve order: brain before skull before scalp (priority ascending).
    if geom == "mird_head":
        order = list(regions.sort_values("priority")["region"])
        add("head carve order brain<skull<scalp", order == ["brain", "skull", "scalp"],
            " -> ".join(order))

    # --- report --------------------------------------------------------------
    width = max(len(n) for n, _, _ in checks)
    n_fail = 0
    print(f"check_run: {args.run_dir}  (geometry: {geom})")
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{width}}  {detail}")
        n_fail += not ok
    if n_fail:
        print(f"\n{n_fail} check(s) FAILED")
        sys.exit(1)
    print(f"\nall {len(checks)} checks passed")


if __name__ == "__main__":
    main()
