#!/usr/bin/env python3
"""Build the documented analytic reference for the folding procedure."""

from __future__ import annotations

import argparse
import json
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from decay_sampling.scenarios import resolve_scenario

from .exposure_convergence import (
    coarsen_exposure,
    compare_folding_results,
    regular_coarse_edges,
)
from .exposure_folding import (
    CrossSectionEnsemble,
    ExposureMetadata,
    fold_exposure,
    sha256,
    validate_exposure_table,
    write_result,
)
from .make_folding_plots import generate as generate_plots


REFERENCE_PARAMETERS = {
    "description": (
        "Analytic longitudinal proton-exposure field used to demonstrate "
        "cross-section folding; it is not a Geant4 treatment result."),
    "energy_min_MeV": 5.0,
    "energy_max_MeV": 120.0,
    "fine_energy_width_MeV": 0.5,
    "depth_min_mm": 0.0,
    "depth_max_mm": 122.0,
    "depth_width_mm": 2.0,
    "reference_range_mm": 112.0,
    "entrance_mean_energy_MeV": 112.0,
    "energy_spread_MeV": 3.0,
    "range_energy_exponent": 0.56,
    "target_exposure_scales": {"C12": 0.12, "N14": 0.03, "O16": 0.55},
    "exposure_scale_cm2_inv": 1.0e27,
    "n_protons": 1.0e6,
    "target_dose_Gy": 0.01,
    "scenario": "inroom",
    "candidate_energy_widths_MeV": [0.5, 1.0, 2.0, 5.0],
    "relative_tolerance_fraction": 0.10,
    "r50_absolute_tolerance_mm": 0.10,
}


def analytic_exposure() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Return a deterministic target-resolved exposure with a distal fall-off."""

    parameters = REFERENCE_PARAMETERS
    energy_edges = np.arange(
        parameters["energy_min_MeV"],
        parameters["energy_max_MeV"]
        + 0.5 * parameters["fine_energy_width_MeV"],
        parameters["fine_energy_width_MeV"],
    )
    depth_edges = np.arange(
        parameters["depth_min_mm"],
        parameters["depth_max_mm"] + 0.5 * parameters["depth_width_mm"],
        parameters["depth_width_mm"],
    )
    energy_centres = 0.5 * (energy_edges[:-1] + energy_edges[1:])
    depth_centres = 0.5 * (depth_edges[:-1] + depth_edges[1:])
    rows = []
    for depth_index, depth in enumerate(depth_centres):
        residual_fraction = max(
            (parameters["reference_range_mm"] - depth)
            / parameters["reference_range_mm"],
            0.0,
        )
        if residual_fraction > 0:
            mean_energy = (
                parameters["energy_min_MeV"]
                + (parameters["entrance_mean_energy_MeV"]
                   - parameters["energy_min_MeV"])
                * residual_fraction ** parameters["range_energy_exponent"]
            )
            energy_weights = np.exp(
                -0.5
                * ((energy_centres - mean_energy)
                   / parameters["energy_spread_MeV"]) ** 2)
            energy_weights[energy_weights < 1.0e-8 * energy_weights.max()] = 0.0
            energy_weights /= energy_weights.sum()
            fluence = 1.0 - 0.12 * depth / parameters["reference_range_mm"]
        else:
            energy_weights = np.zeros_like(energy_centres)
            energy_weights[0] = 1.0
            fluence = 0.20
        occupied = np.flatnonzero(energy_weights > 0)
        for target, target_scale in parameters["target_exposure_scales"].items():
            for energy_index in occupied:
                rows.append({
                    "target": target,
                    "energy_low_MeV": energy_edges[energy_index],
                    "energy_high_MeV": energy_edges[energy_index + 1],
                    "energy_mean_MeV": energy_centres[energy_index],
                    "depth_low_mm": depth_edges[depth_index],
                    "depth_high_mm": depth_edges[depth_index + 1],
                    "depth_mean_mm": depth,
                    "target_exposure_cm2_inv": (
                        parameters["exposure_scale_cm2_inv"]
                        * target_scale * fluence * energy_weights[energy_index]),
                })
    return validate_exposure_table(pd.DataFrame(rows)), energy_edges, depth_edges


def _metadata(energy_edges: np.ndarray, depth_edges: np.ndarray) -> ExposureMetadata:
    parameters = REFERENCE_PARAMETERS
    n_protons = float(parameters["n_protons"])
    target_dose = float(parameters["target_dose_Gy"])
    return ExposureMetadata(
        schema_version=1,
        run_id="analytic_folding_reference_v1",
        exposure_file="analytic_reference_exposure.csv",
        exposure_sha256="generated-in-memory",
        n_protons=n_protons,
        target_dose_Gy=target_dose,
        Np_per_Gy=n_protons / target_dose,
        physics_list="analytic-reference",
        geant4_version="not-applicable",
        random_seed=0,
        software_revision="tracked-generator",
        beam_axis="z",
        depth_origin="analytic reference entrance",
        depth_unit="mm",
        energy_edges_MeV=energy_edges,
        depth_edges_mm=depth_edges,
    )


def generate(repo: Path) -> list[Path]:
    """Run the real folding and convergence code on the analytic reference."""

    repo = Path(repo).resolve()
    exposure, energy_edges, depth_edges = analytic_exposure()
    metadata = _metadata(energy_edges, depth_edges)
    ensemble = CrossSectionEnsemble.from_fit_directory(
        repo / "data/xsections/fits")
    scenario = resolve_scenario(
        REFERENCE_PARAMETERS["scenario"],
        repo / "config/handoff_scenarios.toml")
    fine_result = fold_exposure(exposure, ensemble, metadata, scenario)

    convergence_rows = []
    for width in REFERENCE_PARAMETERS["candidate_energy_widths_MeV"]:
        coarse_edges = regular_coarse_edges(energy_edges, float(width))
        coarse_exposure = coarsen_exposure(exposure, coarse_edges)
        coarse_result = fold_exposure(
            coarse_exposure,
            ensemble,
            replace(metadata, energy_edges_MeV=coarse_edges),
            scenario,
        )
        comparison = compare_folding_results(
            fine_result,
            coarse_result,
            relative_tolerance_fraction=float(
                REFERENCE_PARAMETERS["relative_tolerance_fraction"]),
            r50_absolute_tolerance_mm=float(
                REFERENCE_PARAMETERS["r50_absolute_tolerance_mm"]),
        )
        comparison.insert(0, "candidate_width_MeV", float(width))
        convergence_rows.append(comparison)
    convergence = pd.concat(convergence_rows, ignore_index=True)

    figure_directory = repo / "docs/figures/xsection_folding/reference"
    generated_directory = repo / "docs/generated/xsection_folding/reference"
    figure_directory.mkdir(parents=True, exist_ok=True)
    generated_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ptcrysp-folding-reference-") as temporary:
        folding_directory = Path(temporary) / "folding"
        folding_directory.mkdir()
        write_result(
            fine_result,
            folding_directory,
            exposure_metadata=metadata,
            cross_sections=ensemble,
            scenario=scenario,
            write_replica_profiles=False,
        )
        convergence_path = Path(temporary) / "exposure_convergence.csv"
        convergence.to_csv(convergence_path, index=False, lineterminator="\n")
        figures = generate_plots(
            folding_directory,
            figure_directory,
            convergence_csv=convergence_path,
            generated_directory=generated_directory,
            context_label="Analytic exposure reference (not a transport result)",
            source_label="analytic_folding_reference_v1",
        )

    fine_result.nominal_isotope_profiles.to_csv(
        generated_directory / "nominal_profiles.csv", index=False,
        lineterminator="\n")
    fine_result.profile_bands.to_csv(
        generated_directory / "profile_bands.csv", index=False,
        lineterminator="\n")
    fine_result.production_summary.to_csv(
        generated_directory / "production_summary.csv", index=False,
        lineterminator="\n")
    fine_result.uncertainty_summary.to_csv(
        generated_directory / "uncertainty_summary.csv", index=False,
        lineterminator="\n")
    convergence.to_csv(
        generated_directory / "exposure_convergence.csv", index=False,
        lineterminator="\n")
    reference_metadata = {
        "schema_version": 1,
        "parameters": REFERENCE_PARAMETERS,
        "fit_meta_sha256": sha256(repo / "data/xsections/fits/fit_meta.json"),
        "handoff_scenario_sha256": scenario.config_sha256,
        "exposure_rows": len(exposure),
        "replicas": len(ensemble.replica_ids),
        "all_candidates_pass": bool(convergence["pass"].all()),
    }
    (generated_directory / "reference_definition.json").write_text(
        json.dumps(reference_metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return figures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo", type=Path, default=Path(__file__).resolve().parents[2])
    arguments = parser.parse_args()
    figures = generate(arguments.repo)
    print(f"generated {len(figures)} analytic folding-reference figures")


if __name__ == "__main__":
    main()
