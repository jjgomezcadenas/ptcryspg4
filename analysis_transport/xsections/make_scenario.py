"""Assemble a detector-consumable scenario from a data-driven sampling run.

The sampling run of stageA_xsection_ensemble writes emitters.csv,
depth_dose.csv and the exposure metadata.  This tool adds the remaining
scenario files in the stageA contract --- run_meta.csv, phantom_regions.csv,
isotopes.csv --- with additive provenance columns naming the production
model, and computes the per-isotope sampling budgets for the named handoff
scenarios with the existing decay_sampling/budget.py.  The result is read by
check_run.py, tools/snapshot_scenario.py and the PTCryspMC.jl scenario
reader without any change on their side.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from common.isotopes import ISOTOPES

from .exposure_folding import sha256

# World-frame medium regions per supported geometry (the stageA
# phantom_regions.csv rows, replicated from the geometry constants).
GEOMETRY_REGIONS = {
    "uniform_headep": [{
        "region": "head", "priority": 0, "material": "G4_BRAIN_ICRP",
        "solid": "ellipsoid", "a_mm": 72.0, "b_mm": 87.0, "c_mm": 102.0,
        "cx_mm": 0.0, "cy_mm": -30.0, "cz_mm": 0.0,
        "euler_x_deg": 0.0, "euler_y_deg": 0.0, "euler_z_deg": 0.0,
    }],
}
# Pencil-gun constants reported in run_meta (stageA convention; an SOBP
# run's beam record is the copied layer table).
BEAM_ENERGY_MEV = 100.0
BEAM_SIGMA_MM = 3.0
MEV_TO_J = 1.602176634e-13

RUN_META_FIELDS = (
    "n_protons", "beam_energy_MeV", "beam_sigma_mm", "geometry",
    "phantom_material", "phantom_diameter_mm", "phantom_length_mm",
    "phantom_mass_g", "edep_total_MeV", "dose_total_Gy", "target_dose_Gy",
    "target_mass_g", "target_radius_mm", "target_prox_depth_mm",
    "target_dist_depth_mm", "Np_per_Gy", "geant4_version", "physics_list",
    "random_seed",
    # additive provenance of the data-driven production
    "production_model", "sampling_curves_sha256", "fit_meta_sha256",
)


def assemble(run_dir: Path, repo: Path, scenarios=("d120s300", "d120s120", "d180s300", "d180s120", "d300s300")):
    meta = json.loads((run_dir / "run_meta_raw.json").read_text())
    geometry = meta["geometry"]
    if geometry not in GEOMETRY_REGIONS:
        raise ValueError(f"no region table for geometry '{geometry}'")
    regions = GEOMETRY_REGIONS[geometry]

    # Bounding box from the regions (stageA convention: transverse diameter,
    # beam-axis length).
    hx = max(abs(r["cx_mm"]) + r["a_mm"] for r in regions)
    hy = max(abs(r["cy_mm"]) + r["b_mm"] for r in regions)
    hz = max(abs(r["cz_mm"]) + r["c_mm"] for r in regions)

    mass_g = float(meta["phantom_mass_g"])
    edep_MeV = float(meta["edep_total_MeV"])
    dose_total_Gy = edep_MeV * MEV_TO_J / (mass_g * 1.0e-3)

    sampling = json.loads((run_dir / "sampling_meta.json").read_text())
    folding = json.loads((run_dir / "folding/folding_meta.json").read_text())

    row = {
        "n_protons": int(meta["n_protons"]),
        "beam_energy_MeV": BEAM_ENERGY_MEV,
        "beam_sigma_mm": BEAM_SIGMA_MM,
        "geometry": geometry,
        "phantom_material": regions[0]["material"] if len(regions) == 1 else "multi",
        "phantom_diameter_mm": 2.0 * max(hx, hy),
        "phantom_length_mm": 2.0 * hz,
        "phantom_mass_g": mass_g,
        "edep_total_MeV": edep_MeV,
        "dose_total_Gy": dose_total_Gy,
        "target_dose_Gy": float(meta["target_dose_Gy"]),
        "target_mass_g": float(meta["target_mass_g"]),
        "target_radius_mm": float(meta["target_radius_mm"]),
        "target_prox_depth_mm": float(meta["target_prox_mm"]),
        "target_dist_depth_mm": float(meta["target_dist_mm"]),
        "Np_per_Gy": float(meta["n_protons"]) / float(meta["target_dose_Gy"]),
        "geant4_version": meta["geant4_version"],
        "physics_list": meta["physics_list"],
        "random_seed": int(meta["random_seed"]),
        "production_model": "data-driven-nominal",
        "sampling_curves_sha256": sha256(Path(sampling["curves_file"])),
        "fit_meta_sha256": folding["fit_meta_sha256"],
    }
    with (run_dir / "run_meta.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RUN_META_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)

    with (run_dir / "phantom_regions.csv").open("w", newline="",
                                                encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(regions[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(regions)

    with (run_dir / "isotopes.csv").open("w", encoding="utf-8") as f:
        f.write("isotope_id,name,half_life_s,endpoint_MeV,prompt_gamma\n")
        for iid in sorted(ISOTOPES):
            iso = ISOTOPES[iid]
            f.write(f"{iid},{iso.name},{iso.half_life_s:g},"
                    f"{iso.endpoint_MeV:g},{int(iso.prompt_gamma)}\n")

    for scenario in scenarios:
        subprocess.run(
            [sys.executable, str(repo / "decay_sampling/budget.py"),
             str(run_dir), "--scenario", scenario],
            check=True, capture_output=True, text=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    row = assemble(args.run_dir.resolve(), args.repo.resolve())
    print(f"scenario files written for {row['geometry']} "
          f"({row['production_model']}): run_meta, phantom_regions, isotopes, "
          f"budgets")


if __name__ == "__main__":
    main()
