#!/usr/bin/env python3
"""Freeze a finished Stage-A run directory into a scenario snapshot.

A snapshot is tied to a *specific run*: you name a run tag (the directory under
data/runs/, e.g. cylinder_sobp_1e7) and the run's identity, file set, and
figures all come from that directory — nothing is hand-typed or pulled from a
shared scratch dir. Copies the run's source files (whichever exist) and its
figures/ into <dest>/scenarios/<name>/, writes isotopes.csv (from
common/isotopes.py), the file-format note (SCHEMA.md), the documentation PDFs
(from latex/, into docs/), and a scenario README with the run's numbers filled
in. Each snapshot is standalone: it carries everything needed to read it without
the rest of this repo.

The doc PDFs are copied, not built: run `python latex/build_latex.py` first so
they exist in latex/. A missing PDF is warned about, not fatal. Figures are made
by `python analysis_transport/make_figures.py <run_dir>` — run that first.

Usage:
    python tools/snapshot_scenario.py <run_tag> [--name NAME] [--runs-dir DIR]
                                      [--dest DIR]
"""

import argparse
import glob
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

# Source files copied from the run dir if present. The exact set varies by run
# (a head pencil run has no sobp_layers; budgets exist only once budget.py ran),
# so each is optional -- we copy what is there, not a hardcoded must-have list.
CORE_FILES = [
    "emitters.csv", "run_meta.csv", "depth_dose.csv",
    "phantom_regions.csv", "sobp_layers.csv",
]
# Plus every sampling_budget_*.csv (+ _meta) found in the run dir (glob).

# Documentation PDFs copied from latex/ into the snapshot's docs/ (build first
# with latex/build_latex.py). The full doc set travels with the frozen scenario.
LATEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "latex")
DOC_PDFS = [
    "01_user_guide.pdf", "02_beam_design.pdf",
    "03_decay_kinetics.pdf", "04_source_reference.pdf",
]


def collect_data_files(run_dir):
    """Names of the source CSVs that actually exist in the run dir."""
    files = [fn for fn in CORE_FILES if os.path.exists(os.path.join(run_dir, fn))]
    for path in sorted(glob.glob(os.path.join(run_dir, "sampling_budget_*.csv"))):
        files.append(os.path.basename(path))
    return files


def write_isotopes_csv(path):
    rows = ["isotope_id,name,half_life_s,endpoint_MeV,prompt_gamma"]
    for iid in sorted(ISOTOPES):
        iso = ISOTOPES[iid]
        rows.append(f"{iid},{iso.name},{iso.half_life_s:g},{iso.endpoint_MeV:g},"
                    f"{int(iso.prompt_gamma)}")
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")


def read_budget_legend(run_dir):
    """(scenario, t_irr, t_del, t_meas) for each sampling_budget_*_meta.csv."""
    legend = []
    for path in sorted(glob.glob(os.path.join(run_dir, "sampling_budget_*_meta.csv"))):
        bm = pd.read_csv(path).iloc[0]
        legend.append((str(bm["scenario"]), float(bm["t_irr_s"]),
                       float(bm["t_del_s"]), float(bm["t_meas_s"])))
    return legend


def write_readme(path, name, meta, emit, legend, data_files):
    """Scenario README, numbers filled in from run_meta.csv and emitters.csv."""
    m = meta
    geometry = str(m["geometry"])
    has_sobp = "sobp_layers.csv" in data_files
    t_dose = float(m["target_dose_Gy"])
    counts = emit["isotope_id"].value_counts()
    n_emit = int(len(emit))  # raw produced rows = the source's spatial sample size
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

    # Phantom + beam descriptions adapt to the run (cylinder vs head; SOBP vs pencil).
    if geometry == "cylinder":
        phantom_desc = (f"{m['phantom_material']} cylinder, {diam:g} cm diameter "
                        f"x {length:g} cm")
    else:
        phantom_desc = (f"{geometry} ({m['phantom_material']}; regions in "
                        f"phantom_regions.csv), bounding box {diam:g} x {length:g} cm")
    # "multi" (a heterogeneous head) reads awkwardly in prose; name the geometry.
    medium = m['phantom_material'] if str(m['phantom_material']) != "multi" \
        else f"the {geometry} phantom"
    if has_sobp:
        beam_desc = ("proton SOBP (energy layers in sobp_layers.csv), uniform disk "
                     "over the target")
        beam_title = f"Proton SOBP on {medium}"
    else:
        beam_desc = (f"single proton pencil, {float(m['beam_energy_MeV']):g} MeV, "
                     f"sigma {float(m['beam_sigma_mm']):g} mm")
        beam_title = f"Proton pencil on {medium}"

    legrows = "\n".join(
        f"- `{s}`: t_irr={ti:g}, t_del={td:g}, t_meas={tm:g} s" for s, ti, td, tm in legend)
    if not legrows:
        legrows = "- (no timing budget frozen yet; run decay_sampling/budget.py)"

    files_line = ", ".join(data_files)
    # The Parodi Table 2 yield check is for the SOBP target-volume delivery; a
    # single pencil paints a thin track, so its integral yield is not comparable.
    parodi_note = ("Parodi 2008 Table 2 comparison: about 2.2x overall." if has_sobp
                   else "(Single-pencil run: thin-track yields, not comparable to "
                        "the Parodi Table 2 SOBP integrals.)")

    text = f"""# {name}

{beam_title}, {int(m['n_protons']):g} protons, scaled to 1 Gy in the target box.
Detector-independent positron-annihilation source for a downstream PET simulation.

## Setup
- Geometry: {geometry}
- Phantom: {phantom_desc}
- Beam: {beam_desc}
- Target box: {box_diam:g} cm diameter x {box_len:g} cm, at {box_lo:g}-{box_hi:g} cm depth
- Protons: {int(m['n_protons']):g}
- Geant4 {m['geant4_version']}, {m['physics_list']}, seed {int(m['random_seed'])}

## Coordinate frame
Phantom centred at the origin, axis along +z (the beam direction). It spans
z in [{-half_z:g}, {half_z:g}] mm and r <= {float(m['phantom_diameter_mm'])/2:g} mm; the beam enters at
z = {-half_z:g} mm. Depth from the entrance maps as z = depth - {half_z:g} mm, so the
target box sits at z in [{z_lo:g}, {z_hi:g}] mm. emitters.csv positions are in this frame.

## Normalization — where the dose comes from
The clinical dose is set at the handoff. Three quantities, three roles:

- **Protons simulated ({int(m['n_protons']):g})** — the Monte-Carlo statistics
  knob: more protons means less noise. The dose and the per-Gy yields are ratios
  in which the proton count cancels, so they hold at any statistics and Stage A
  runs once, reused for every detector.
- **Target dose ({t_dose:.2e} Gy)** — a measured output of the run: the dose these
  protons deposited in the target box.
- **Clinical dose D (1 Gy by default)** — the prescription, applied at the handoff
  by `budget.py --dose`. Yields are linear in dose, so each produced count is
  rescaled `P_j(D) = count_j · D / target_dose`. Delivering 1 Gy to the box takes
  Np(1 Gy) = {float(m['Np_per_Gy']):.3e} protons.

## Yields per 1 Gy
{ystr} (total {total:.2e})
{parodi_note}

## Timing scenarios (budgets)
The acquisition timing is the only difference between budgets; the spatial
source (emitters.csv) is shared. N_expected is for 1 Gy.
{legrows}

## How to use this source
emitters.csv holds {n_emit:g} annihilation points: the produced source's spatial
shape at this run's target dose. The measured count for the clinical dose is
N_expected (the budget below), realized by the Poisson draw in step 2.

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
{files_line},
phantom_material_*.csv (+ _meta: composition, density, mu at 511 keV per region,
for gamma transport + attenuation correction), figures/. Columns in SCHEMA.md,
isotope codes in isotopes.csv. The full documentation (user guide, beam design,
decay kinetics, source reference) is in docs/ as PDFs.
"""
    with open(path, "w") as f:
        f.write(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_tag", help="run directory under --runs-dir, e.g. cylinder_sobp_1e7")
    ap.add_argument("--name", default=None,
                    help="scenario name (default: the run tag)")
    ap.add_argument("--runs-dir", default=os.path.join(_HERE, "..", "data", "runs"))
    ap.add_argument("--dest", default=os.path.expanduser("~/Projects/ptcrysp-scenarios"))
    args = ap.parse_args()

    run_dir = os.path.join(args.runs_dir, args.run_tag)
    if not os.path.isdir(run_dir):
        sys.exit(f"run directory not found: {run_dir}")
    name = args.name or args.run_tag

    out = os.path.join(args.dest, "scenarios", name)
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)

    # Source CSVs: copy whatever this run actually produced.
    data_files = collect_data_files(run_dir)
    for fn in data_files:
        shutil.copy2(os.path.join(run_dir, fn), os.path.join(out, fn))

    # Figures: copy the run's figures/ verbatim (made by make_figures.py).
    src_figdir = os.path.join(run_dir, "figures")
    fig_names = []
    if os.path.isdir(src_figdir):
        for fn in sorted(os.listdir(src_figdir)):
            src = os.path.join(src_figdir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(figdir, fn))
                fig_names.append(fn)
    if not fig_names:
        print("WARNING: no figures in the run dir; run "
              "`python analysis_transport/make_figures.py <run_dir>` first.")

    write_isotopes_csv(os.path.join(out, "isotopes.csv"))
    shutil.copy2(os.path.join(_HERE, "scenario_template", "SCHEMA.md"),
                 os.path.join(out, "SCHEMA.md"))

    meta = pd.read_csv(os.path.join(run_dir, "run_meta.csv")).iloc[0]
    emit = pd.read_csv(os.path.join(run_dir, "emitters.csv"))
    legend = read_budget_legend(run_dir)
    write_readme(os.path.join(out, "README.md"), name, meta, emit, legend, data_files)

    # Phantom medium for 511 keV gamma transport + reconstruction attenuation
    # correction: composition + mu for each distinct material in
    # phantom_regions.csv (one homogeneous material for the cylinder; brain/bone/
    # soft tissue for the head).
    regions = pd.read_csv(os.path.join(run_dir, "phantom_regions.csv"))
    n_mat = 0
    for mat_name in dict.fromkeys(regions["material"].astype(str)):  # distinct, ordered
        if mat_name in MATERIALS:
            write_material_csv(MATERIALS[mat_name], ANNIHILATION_keV, out)
            n_mat += 2
        else:
            print(f"WARNING: material '{mat_name}' not in the registry; no "
                  f"phantom_material file (add it to common/phantom_material.py).")

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

    n = len(data_files) + len(fig_names) + n_pdf + 3 + n_mat
    print(f"snapshot '{name}' from {run_dir}: wrote {n} files -> {out}")


if __name__ == "__main__":
    main()
