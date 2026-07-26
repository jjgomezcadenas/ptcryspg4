#!/usr/bin/env python3
"""Coarsen a fine exposure table exactly and test folding-grid convergence."""

from __future__ import annotations

import argparse
import json
import math
import tomllib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from decay_sampling.scenarios import DEFAULT_SCENARIO_CONFIG, resolve_scenario

from .exposure_folding import (
    CrossSectionEnsemble,
    ExposureMetadata,
    FoldingResult,
    distal_r50,
    fold_exposure,
    load_exposure_metadata,
    load_exposure_table,
    validate_exposure_table,
)


def _containing_bins(low, high, edges, label: str) -> np.ndarray:
    grid = np.asarray(edges, dtype=float)
    if grid.ndim != 1 or len(grid) < 2 or np.any(np.diff(grid) <= 0):
        raise ValueError(f"coarse {label} edges must be strictly increasing")
    result = np.empty(len(low), dtype=int)
    for row, (lower, upper) in enumerate(zip(low, high, strict=True)):
        matches = np.flatnonzero(
            (lower >= grid[:-1] - 1.0e-10)
            & (upper <= grid[1:] + 1.0e-10)
        )
        if len(matches) != 1:
            raise ValueError(
                f"fine {label} bin [{lower:g}, {upper:g}] does not lie wholly "
                "inside one coarse bin")
        result[row] = int(matches[0])
    return result


def coarsen_exposure(
    fine_exposure: pd.DataFrame,
    coarse_energy_edges_MeV,
    coarse_depth_edges_mm=None,
) -> pd.DataFrame:
    """Recombine exposure and its first moments without information loss."""

    fine = validate_exposure_table(fine_exposure).copy()
    energy_edges = np.asarray(coarse_energy_edges_MeV, dtype=float)
    if coarse_depth_edges_mm is None:
        depth_edges = np.unique(np.concatenate([
            fine["depth_low_mm"].to_numpy(dtype=float),
            fine["depth_high_mm"].to_numpy(dtype=float),
        ]))
    else:
        depth_edges = np.asarray(coarse_depth_edges_mm, dtype=float)
    fine["_energy_index"] = _containing_bins(
        fine["energy_low_MeV"].to_numpy(dtype=float),
        fine["energy_high_MeV"].to_numpy(dtype=float),
        energy_edges,
        "energy",
    )
    fine["_depth_index"] = _containing_bins(
        fine["depth_low_mm"].to_numpy(dtype=float),
        fine["depth_high_mm"].to_numpy(dtype=float),
        depth_edges,
        "depth",
    )
    fine["_energy_moment"] = (
        fine["target_exposure_cm2_inv"] * fine["energy_mean_MeV"])
    fine["_depth_moment"] = (
        fine["target_exposure_cm2_inv"] * fine["depth_mean_mm"])
    grouped = (
        fine.groupby(["target", "_energy_index", "_depth_index"], sort=True)
        .agg(
            target_exposure_cm2_inv=("target_exposure_cm2_inv", "sum"),
            energy_moment=("_energy_moment", "sum"),
            depth_moment=("_depth_moment", "sum"),
        )
        .reset_index()
    )
    grouped = grouped.loc[grouped["target_exposure_cm2_inv"] > 0].copy()
    exposure = grouped["target_exposure_cm2_inv"].to_numpy(dtype=float)
    energy_index = grouped["_energy_index"].to_numpy(dtype=int)
    depth_index = grouped["_depth_index"].to_numpy(dtype=int)
    result = pd.DataFrame({
        "target": grouped["target"],
        "energy_low_MeV": energy_edges[energy_index],
        "energy_high_MeV": energy_edges[energy_index + 1],
        "energy_mean_MeV": grouped["energy_moment"] / exposure,
        "depth_low_mm": depth_edges[depth_index],
        "depth_high_mm": depth_edges[depth_index + 1],
        "depth_mean_mm": grouped["depth_moment"] / exposure,
        "target_exposure_cm2_inv": exposure,
    })
    return validate_exposure_table(result)


def regular_coarse_edges(fine_edges, width: float) -> np.ndarray:
    """Build a regular candidate grid while retaining both fine endpoints."""

    fine = np.asarray(fine_edges, dtype=float)
    if width <= 0:
        raise ValueError("candidate energy width must be positive")
    start, stop = float(fine[0]), float(fine[-1])
    edges = start + width * np.arange(int(np.floor((stop - start) / width)) + 1)
    if not np.isclose(edges[-1], stop, atol=1.0e-10, rtol=0):
        edges = np.append(edges, stop)
    else:
        edges[-1] = stop
    return edges


def _half_width(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return float("nan")
    q16, q84 = np.quantile(finite, [0.16, 0.84])
    return float(0.5 * (q84 - q16))


def compare_folding_results(
    fine: FoldingResult,
    coarse: FoldingResult,
    *,
    relative_tolerance_fraction: float = 0.10,
    r50_absolute_tolerance_mm: float = 0.10,
) -> pd.DataFrame:
    """Apply paired replica-relative yield and distal-edge criteria."""

    if fine.profile_labels != coarse.profile_labels:
        raise ValueError("fine and coarse foldings have different profile labels")
    if fine.replica_profile_values.shape[:2] != coarse.replica_profile_values.shape[:2]:
        raise ValueError("fine and coarse foldings have different replica ensembles")
    if relative_tolerance_fraction <= 0 or r50_absolute_tolerance_mm <= 0:
        raise ValueError("convergence tolerances must be positive")

    rows = []
    for position, label in enumerate(fine.profile_labels):
        fine_yields = fine.replica_profile_values[:, position, :].sum(axis=1)
        coarse_yields = coarse.replica_profile_values[:, position, :].sum(axis=1)
        fine_nominal_yield = float(fine.nominal_profile_values[position].sum())
        coarse_nominal_yield = float(coarse.nominal_profile_values[position].sum())
        yield_half_width = _half_width(fine_yields)
        yield_tolerance = relative_tolerance_fraction * yield_half_width
        max_yield_change = float(max(
            abs(coarse_nominal_yield - fine_nominal_yield),
            np.max(np.abs(coarse_yields - fine_yields)),
        ))
        coarse_yield_half_width = _half_width(coarse_yields)
        yield_band_change = abs(coarse_yield_half_width - yield_half_width)

        fine_r50 = np.asarray([
            distal_r50(fine.depth_mm, values)
            for values in fine.replica_profile_values[:, position, :]
        ])
        coarse_r50 = np.asarray([
            distal_r50(coarse.depth_mm, values)
            for values in coarse.replica_profile_values[:, position, :]
        ])
        fine_nominal_r50 = distal_r50(
            fine.depth_mm, fine.nominal_profile_values[position])
        coarse_nominal_r50 = distal_r50(
            coarse.depth_mm, coarse.nominal_profile_values[position])
        fine_r50_all = np.concatenate([[fine_nominal_r50], fine_r50])
        coarse_r50_all = np.concatenate([[coarse_nominal_r50], coarse_r50])
        jointly_finite = np.isfinite(fine_r50_all) & np.isfinite(coarse_r50_all)
        finiteness_matches = bool(np.array_equal(
            np.isfinite(fine_r50_all), np.isfinite(coarse_r50_all)))
        if np.any(jointly_finite):
            max_r50_change = float(np.max(np.abs(
                coarse_r50_all[jointly_finite] - fine_r50_all[jointly_finite])))
        else:
            max_r50_change = 0.0 if finiteness_matches else float("inf")
        fine_shifts = fine_r50 - fine_nominal_r50
        r50_half_width = _half_width(fine_shifts)
        if not math.isfinite(r50_half_width) and finiteness_matches:
            r50_half_width = 0.0
        r50_relative_limit = relative_tolerance_fraction * r50_half_width
        r50_tolerance = min(r50_absolute_tolerance_mm, r50_relative_limit)
        coarse_half_width = _half_width(coarse_r50 - coarse_nominal_r50)
        if not math.isfinite(coarse_half_width) and finiteness_matches:
            coarse_half_width = 0.0
        r50_band_change = abs(coarse_half_width - r50_half_width)

        numerical_floor = 1.0e-12
        yield_pass = max_yield_change <= max(yield_tolerance, numerical_floor)
        r50_pass = (
            finiteness_matches
            and max_r50_change <= max(r50_tolerance, numerical_floor)
        )
        yield_band_pass = yield_band_change <= max(yield_tolerance, numerical_floor)
        r50_band_pass = r50_band_change <= max(r50_tolerance, numerical_floor)
        rows.append({
            "profile_label": label,
            "max_paired_yield_change_run": max_yield_change,
            "yield_replica_half_width_run": yield_half_width,
            "yield_tolerance_run": yield_tolerance,
            "yield_half_width_change_run": yield_band_change,
            "max_paired_R50_change_mm": max_r50_change,
            "R50_replica_half_width_mm": r50_half_width,
            "R50_tolerance_mm": r50_tolerance,
            "R50_half_width_change_mm": r50_band_change,
            "yield_pass": yield_pass,
            "R50_pass": r50_pass,
            "yield_band_pass": yield_band_pass,
            "R50_band_pass": r50_band_pass,
            "pass": yield_pass and r50_pass and yield_band_pass and r50_band_pass,
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exposure_csv", type=Path)
    parser.add_argument("--exposure-meta", type=Path, required=True)
    parser.add_argument("--fit-dir", type=Path, default=Path("data/xsections/fits"))
    parser.add_argument(
        "--config", type=Path,
        default=Path("config/xsection_exposure_convergence.toml"))
    parser.add_argument("--scenario", default="d120s300")
    parser.add_argument("--scenario-config", type=Path, default=DEFAULT_SCENARIO_CONFIG)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()

    settings = tomllib.loads(arguments.config.read_text(encoding="utf-8"))
    if int(settings.get("schema_version", 0)) != 1:
        raise ValueError("unsupported exposure-convergence schema version")
    fine_exposure = load_exposure_table(arguments.exposure_csv)
    metadata = load_exposure_metadata(arguments.exposure_meta, arguments.exposure_csv)
    ensemble = CrossSectionEnsemble.from_fit_directory(arguments.fit_dir)
    scenario = resolve_scenario(arguments.scenario, arguments.scenario_config)
    fine_result = fold_exposure(fine_exposure, ensemble, metadata, scenario)
    rows = []
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    candidate_records = []
    for width in settings["candidate_energy_widths_MeV"]:
        edges = regular_coarse_edges(metadata.energy_edges_MeV, float(width))
        coarse_exposure = coarsen_exposure(fine_exposure, edges)
        width_tag = f"{float(width):g}".replace(".", "p")
        coarse_name = f"exposure_coarse_{width_tag}MeV.csv"
        coarse_exposure.to_csv(arguments.output_dir / coarse_name, index=False)
        coarse_metadata = replace(metadata, energy_edges_MeV=edges)
        coarse_result = fold_exposure(
            coarse_exposure, ensemble, coarse_metadata, scenario)
        comparison = compare_folding_results(
            fine_result,
            coarse_result,
            relative_tolerance_fraction=float(
                settings["relative_tolerance_fraction"]),
            r50_absolute_tolerance_mm=float(
                settings["r50_absolute_tolerance_mm"]),
        )
        comparison.insert(0, "candidate_width_MeV", float(width))
        rows.append(comparison)
        candidate_records.append({
            "candidate_width_MeV": float(width),
            "energy_edges_MeV": edges.tolist(),
            "coarsened_exposure_file": coarse_name,
            "pass": bool(comparison["pass"].all()),
        })
    result = pd.concat(rows, ignore_index=True)
    result.to_csv(arguments.output_dir / "exposure_convergence.csv", index=False)
    summary = {
        "schema_version": 1,
        "fine_exposure_sha256": metadata.exposure_sha256,
        "fine_energy_edges_MeV": metadata.energy_edges_MeV.tolist(),
        "candidates": candidate_records,
        "relative_tolerance_fraction": settings["relative_tolerance_fraction"],
        "r50_absolute_tolerance_mm": settings["r50_absolute_tolerance_mm"],
        "all_candidates_pass": bool(result["pass"].all()),
    }
    (arguments.output_dir / "exposure_convergence_meta.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
