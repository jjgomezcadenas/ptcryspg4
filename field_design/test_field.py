#!/usr/bin/env python3
"""Unit checks for the WEPL field design (Step 2): RSP + the WEPL ray-trace.

Pure-function tests, no Geant4: relative stopping-power values, and the WEPL
ray-trace on synthetic homogeneous phantoms where the answer is exact (WEPL =
depth x RSP) plus the real head (bone makes WEPL > geometric). Exits non-zero on
any failure — the test gate for sobp.py's heterogeneous design.

Usage:
    python field_design/test_field.py
"""

import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from phantom_material import MATERIALS, relative_stopping_power as rsp  # noqa: E402
sys.path.insert(0, _HERE)
from sobp import rsp_by_material, wepl_curve  # noqa: E402

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))


def cyl_regions(material):
    """A single water-equivalent cylinder, z in [-80,80] mm on the axis."""
    return pd.DataFrame([dict(
        region="phantom", priority=0, material=material, solid="cylinder",
        a_mm=80.0, b_mm=80.0, c_mm=80.0, cx_mm=0.0, cy_mm=0.0, cz_mm=0.0)])


# --- RSP values ----------------------------------------------------------------
check("RSP(water) = 1", abs(rsp(MATERIALS["G4_WATER"]) - 1.0) < 1e-9)
check("RSP(brain) in [1.0,1.07]", 1.0 <= rsp(MATERIALS["G4_BRAIN_ICRP"]) <= 1.07,
      f"{rsp(MATERIALS['G4_BRAIN_ICRP']):.4f}")
check("RSP(soft) in [1.0,1.06]", 1.0 <= rsp(MATERIALS["G4_TISSUE_SOFT_ICRP"]) <= 1.06,
      f"{rsp(MATERIALS['G4_TISSUE_SOFT_ICRP']):.4f}")
check("RSP(bone) in [1.5,1.8]", 1.5 <= rsp(MATERIALS["G4_BONE_CORTICAL_ICRP"]) <= 1.8,
      f"{rsp(MATERIALS['G4_BONE_CORTICAL_ICRP']):.4f}")
_b100 = rsp(MATERIALS["G4_BONE_CORTICAL_ICRP"], 100.0)
_b200 = rsp(MATERIALS["G4_BONE_CORTICAL_ICRP"], 200.0)
check("RSP weakly energy-dependent", abs(_b100 - _b200) / _b100 < 0.01,
      f"bone {_b100:.4f}(100) vs {_b200:.4f}(200) MeV")

# --- WEPL on homogeneous phantoms is exact (= depth x RSP) ----------------------
for mat in ("G4_WATER", "G4_BRAIN_ICRP"):
    reg = cyl_regions(mat)
    depths, wepl = wepl_curve(reg, -80.0, 100.0, rsp_by_material(reg))
    expect = depths * rsp(MATERIALS[mat])
    err = float(np.max(np.abs(wepl - expect)))
    check(f"WEPL({mat}) = depth x RSP", err < 0.05, f"max err {err:.2e} mm")
    check(f"WEPL({mat}) monotonic", bool(np.all(np.diff(wepl) >= -1e-12)))

# --- WEPL on the real head: bone makes WEPL > geometric at the target -----------
_head_dir = os.path.join(_HERE, "..", "data", "runs", "mird_head_pencil_1e5")
if os.path.exists(os.path.join(_head_dir, "phantom_regions.csv")):
    head = pd.read_csv(os.path.join(_head_dir, "phantom_regions.csv"))
    depths, wepl = wepl_curve(head, -72.0, 105.0, rsp_by_material(head))
    wp = float(np.interp(55.0, depths, wepl))
    check("head WEPL(55mm) > 58 (bone adds)", wp > 58.0, f"WEPL(55mm) = {wp:.1f} mm")
    check("head WEPL monotonic", bool(np.all(np.diff(wepl) >= -1e-12)))
else:
    check("head WEPL (needs mird_head_pencil_1e5 run)", True, "skipped: no run")

# --- report --------------------------------------------------------------------
width = max(len(n) for n, _, _ in CHECKS)
n_fail = sum(not ok for _, ok, _ in CHECKS)
for name, ok, detail in CHECKS:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{width}}  {detail}")
if n_fail:
    print(f"\n{n_fail} check(s) FAILED")
    sys.exit(1)
print(f"\nall {len(CHECKS)} checks passed")
