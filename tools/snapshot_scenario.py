#!/usr/bin/env python3
"""Freeze a finished Stage-A run in data/ into a scenario snapshot.

Copies the run's source files and figures into <dest>/scenarios/<name>/, writes
isotopes.csv (from common/isotopes.py), the file-format note (SCHEMA.md), the
documentation PDFs (from latex/, into docs/), and a scenario README with the
run's numbers filled in. Each snapshot is standalone: it carries everything
needed to read it without the rest of this repo.

The doc PDFs are copied, not built: run `python latex/build_latex.py` first so
they exist in latex/. A missing PDF is warned about, not fatal.

Usage:
    python tools/snapshot_scenario.py <name> [--data-dir DIR] [--dest DIR]
"""

import argparse
import os
import shutil
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402
from phantom_material import MATERIALS, write_material_csv  # noqa: E402

# Annihilation energy [keV] -- the energy the downstream mu/attenuation is for.
ANNIHILATION_keV = 511.0

# Files copied verbatim from data/ into the snapshot.
DATA_FILES = [
    "emitters.csv", "run_meta.csv",
    "sampling_budget_inroom.csv", "sampling_budget_inroom_meta.csv",
    "sampling_budget_fast.csv", "sampling_budget_fast_meta.csv",
    "sampling_budget_offline.csv", "sampling_budget_offline_meta.csv",
    "sobp_layers.csv",
]
FIGURES = ["sobp_g4.png", "transport_validation.png", "activity.png"]

# Documentation PDFs copied from latex/ into the snapshot's docs/ (build first
# with latex/build_latex.py). The full doc set travels with the frozen scenario.
LATEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "latex")
DOC_PDFS = [
    "01_user_guide.pdf", "02_beam_design.pdf",
    "03_decay_kinetics.pdf", "04_source_reference.pdf",
]


def write_isotopes_csv(path):
    rows = ["isotope_id,name,half_life_s,endpoint_MeV,prompt_gamma"]
    for iid in sorted(ISOTOPES):
        iso = ISOTOPES[iid]
        rows.append(f"{iid},{iso.name},{iso.half_life_s:g},{iso.endpoint_MeV:g},"
                    f"{int(iso.prompt_gamma)}")
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")


def read_budget_legend(data_dir):
    """(scenario, t_irr, t_del, t_meas) for each sampling_budget_*_meta.csv."""
    legend = []
    for fn in DATA_FILES:
        if not fn.startswith("sampling_budget_") or not fn.endswith("_meta.csv"):
            continue
        bm = pd.read_csv(os.path.join(data_dir, fn)).iloc[0]
        legend.append((str(bm["scenario"]), float(bm["t_irr_s"]),
                       float(bm["t_del_s"]), float(bm["t_meas_s"])))
    return legend


def write_readme(path, name, meta, emit, legend):
    """Scenario README, numbers filled in from run_meta.csv and emitters.csv."""
    m = meta
    t_dose = float(m["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()
    # Production yield per Gy = count_j / target_dose (the absolute normalization).
    yields = {iid: counts.get(iid, 0) / t_dose for iid in sorted(ISOTOPES)}
    total = sum(yields.values())
    ystr = ", ".join(f"{ISOTOPES[i].name} {yields[i]:.2e}" for i in sorted(ISOTOPES))

    diam = float(m["phantom_diameter_mm"]) / 10.0
    length = float(m["phantom_length_mm"]) / 10.0
    half_z = float(m["phantom_length_mm"]) / 2.0
    box_diam = 2 * float(m["target_radius_mm"]) / 10.0
    box_len = (float(m["target_dist_depth_mm"]) - float(m["target_prox_depth_mm"])) / 10.0
    box_lo = float(m["target_prox_depth_mm"]) / 10.0
    box_hi = float(m["target_dist_depth_mm"]) / 10.0
    z_lo = float(m["target_prox_depth_mm"]) - half_z
    z_hi = float(m["target_dist_depth_mm"]) - half_z

    legrows = "\n".join(
        f"- `{s}`: t_irr={ti:g}, t_del={td:g}, t_meas={tm:g} s" for s, ti, td, tm in legend)

    text = f"""# {name}

Proton SOBP on {m['phantom_material']}, {int(m['n_protons']):g} protons, scaled to
1 Gy in the target box. Detector-independent positron-annihilation source for a
downstream PET simulation.

## Setup
- Phantom: {m['phantom_material']} cylinder, {diam:g} cm diameter x {length:g} cm
- Beam: proton SOBP (energy layers in sobp_layers.csv), uniform disk over the target
- Target box: {box_diam:g} cm diameter x {box_len:g} cm, at {box_lo:g}-{box_hi:g} cm depth
- Protons: {int(m['n_protons']):g}
- Geant4 {m['geant4_version']}, {m['physics_list']}, seed {int(m['random_seed'])}

## Coordinate frame
Phantom centred at the origin, axis along +z (the beam direction). It spans
z in [{-half_z:g}, {half_z:g}] mm and r <= {float(m['phantom_diameter_mm'])/2:g} mm; the beam enters at
z = {-half_z:g} mm. Depth from the entrance maps as z = depth - {half_z:g} mm, so the
target box sits at z in [{z_lo:g}, {z_hi:g}] mm. emitters.csv positions are in this frame.

## Normalization
- Target dose: {t_dose:.2e} Gy for {int(m['n_protons']):g} protons
- Np(1 Gy): {float(m['Np_per_Gy']):.3e} protons

## Yields per 1 Gy
{ystr} (total {total:.2e})
Parodi 2008 Table 2 comparison: about 2.2x overall.

## Timing scenarios (budgets)
The acquisition timing is the only difference between budgets; the spatial
source (emitters.csv) is shared. N_expected is for 1 Gy.
{legrows}

## How to use this source
1. Pick a timing budget; read N_expected per isotope from sampling_budget_<s>.csv.
2. For each isotope j, draw M_j ~ Poisson(N_expected_j) annihilation points by
   sampling (with replacement) the matching isotope_id rows of emitters.csv
   (anh_x/y/z_mm). Seed reproducibly (e.g. master_seed + realization) so every
   detector config sees the identical source.
3. Emit a back-to-back 511 keV photon pair (isotropic, with the chosen
   non-collinearity) from each sampled annihilation point.
4. Transport the photons through the phantom defined by phantom_material_*.csv +
   the geometry above; use the same medium for reconstruction attenuation
   correction.

## Files
emitters.csv, run_meta.csv, sampling_budget_{{inroom,fast,offline}}.csv (+ _meta),
sobp_layers.csv, phantom_material_*.csv (+ _meta: composition, density, mu at
511 keV for gamma transport + attenuation correction), figures/. Columns in
SCHEMA.md, isotope codes in isotopes.csv. The full documentation (user guide,
beam design, decay kinetics, source reference) is in docs/ as PDFs.
"""
    with open(path, "w") as f:
        f.write(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="scenario name, e.g. head_sobp_1e7")
    ap.add_argument("--data-dir", default=os.path.join(_HERE, "..", "data"))
    ap.add_argument("--dest", default=os.path.expanduser("~/Projects/ptcrysp-scenarios"))
    args = ap.parse_args()

    out = os.path.join(args.dest, "scenarios", args.name)
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)

    for fn in DATA_FILES:
        shutil.copy2(os.path.join(args.data_dir, fn), os.path.join(out, fn))
    for fn in FIGURES:
        shutil.copy2(os.path.join(args.data_dir, fn), os.path.join(figdir, fn))

    write_isotopes_csv(os.path.join(out, "isotopes.csv"))
    shutil.copy2(os.path.join(_HERE, "scenario_template", "SCHEMA.md"),
                 os.path.join(out, "SCHEMA.md"))

    meta = pd.read_csv(os.path.join(args.data_dir, "run_meta.csv")).iloc[0]
    emit = pd.read_csv(os.path.join(args.data_dir, "emitters.csv"))
    legend = read_budget_legend(args.data_dir)
    write_readme(os.path.join(out, "README.md"), args.name, meta, emit, legend)

    # Phantom medium for 511 keV gamma transport + reconstruction attenuation
    # correction: composition + mu, derived from the run's phantom_material.
    n_mat = 0
    mat_name = str(meta["phantom_material"])
    if mat_name in MATERIALS:
        write_material_csv(MATERIALS[mat_name], ANNIHILATION_keV, out)
        n_mat = 2
    else:
        print(f"WARNING: phantom_material '{mat_name}' not in the registry; "
              f"no phantom_material_*.csv written (add it to common/phantom_material.py).")

    # Documentation PDFs (built separately by latex/build_latex.py) -> docs/.
    docdir = os.path.join(out, "docs")
    os.makedirs(docdir, exist_ok=True)
    n_pdf = 0
    for fn in DOC_PDFS:
        src = os.path.join(LATEX_DIR, fn)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(docdir, fn))
            n_pdf += 1
        else:
            print(f"WARNING: {fn} not found in latex/; run "
                  f"`python latex/build_latex.py` first to include it.")

    n = len(DATA_FILES) + len(FIGURES) + n_pdf + 3 + n_mat
    print(f"snapshot '{args.name}': wrote {n} files -> {out}")


if __name__ == "__main__":
    main()
