#!/usr/bin/env python3
"""Fold fitted production cross sections with proton exposure from transport.

The physical exposure table will be written by a dedicated Geant4 run.  This
module is independent of Geant4 so its units, interpolation and aggregation can
be validated with analytic synthetic tables first.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .channels import CHANNELS


MB_TO_CM2 = 1.0e-27
SCHEMA_VERSION = 1
EXPOSURE_FIELDS = (
    "target",
    "energy_low_MeV",
    "energy_high_MeV",
    "energy_mean_MeV",
    "depth_low_mm",
    "depth_high_mm",
    "depth_mean_mm",
    "target_exposure_cm2_inv",
)
NUMERIC_EXPOSURE_FIELDS = EXPOSURE_FIELDS[1:]
STEP_FIELDS = (
    "target",
    "proton_weight",
    "target_number_density_cm3",
    "step_length_cm",
    "energy_MeV",
    "depth_mm",
)
ISOTOPES = tuple(dict.fromkeys(channel.residual for channel in CHANNELS))
REPLICA_COLUMN = re.compile(r"^sigma_([0-9]+(?:\.[0-9]+)?)_MeV_mb$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validated_grid(energy, values, label: str) -> tuple[np.ndarray, np.ndarray]:
    grid = np.asarray(energy, dtype=float)
    array = np.asarray(values, dtype=float)
    if grid.ndim != 1 or len(grid) < 2:
        raise ValueError(f"{label}: energy grid must contain at least two points")
    if not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0):
        raise ValueError(f"{label}: energy grid must be finite and increasing")
    if array.shape[-1] != len(grid):
        raise ValueError(f"{label}: cross-section array does not match energy grid")
    if not np.isfinite(array).all() or np.any(array < 0):
        raise ValueError(f"{label}: cross sections must be finite and non-negative")
    return grid, array


def interpolate_curves(
    grid,
    values,
    query,
    *,
    threshold_MeV: float,
    label: str,
) -> np.ndarray:
    """Linearly interpolate one or more curves along their last axis.

    Energies at and below the reaction threshold return zero.  Energies below
    the stored grid are also zero; energies above the fitted range are rejected.
    """

    grid, array = _validated_grid(grid, values, label)
    energies = np.asarray(query, dtype=float)
    if not np.isfinite(energies).all() or np.any(energies < 0):
        raise ValueError(f"{label}: query energies must be finite and non-negative")
    if np.any(energies > grid[-1] + 1.0e-12):
        highest = float(np.max(energies))
        raise ValueError(
            f"{label}: energy {highest:g} MeV exceeds fitted limit {grid[-1]:g} MeV")

    original_shape = energies.shape
    flat = energies.reshape(-1)
    leading_shape = array.shape[:-1]
    result = np.zeros(leading_shape + (len(flat),), dtype=float)
    active = (flat > threshold_MeV) & (flat >= grid[0])
    if np.any(active):
        selected = flat[active]
        upper = np.searchsorted(grid, selected, side="right")
        upper = np.clip(upper, 1, len(grid) - 1)
        lower = upper - 1
        fraction = (selected - grid[lower]) / (grid[upper] - grid[lower])
        low_values = np.take(array, lower, axis=-1)
        high_values = np.take(array, upper, axis=-1)
        interpolated = low_values + (high_values - low_values) * fraction
        result[..., active] = interpolated
    return result.reshape(leading_shape + original_shape)


@dataclass(frozen=True)
class ChannelCurves:
    threshold_MeV: float
    nominal_energy_MeV: np.ndarray
    nominal_sigma_mb: np.ndarray
    replica_energy_MeV: np.ndarray
    replica_sigma_mb: np.ndarray


@dataclass(frozen=True)
class CrossSectionEnsemble:
    channels: dict[str, ChannelCurves]
    replica_ids: np.ndarray
    fit_directory: Path | None = None

    @classmethod
    def from_fit_directory(cls, fit_directory: Path | str) -> "CrossSectionEnsemble":
        fit_directory = Path(fit_directory)
        metadata = json.loads((fit_directory / "fit_meta.json").read_text())
        thresholds = metadata["threshold_MeV"]
        curves: dict[str, ChannelCurves] = {}
        common_replica_ids: np.ndarray | None = None

        for channel in CHANNELS:
            curve = pd.read_csv(fit_directory / f"{channel.channel_id}_curve.csv")
            replicas = pd.read_csv(
                fit_directory / f"{channel.channel_id}_replicas.csv")
            columns = []
            energies = []
            for column in replicas.columns:
                match = REPLICA_COLUMN.match(column)
                if match:
                    columns.append(column)
                    energies.append(float(match.group(1)))
            if not columns:
                raise ValueError(f"{channel.channel_id}: no replica cross-section columns")

            order = np.argsort(energies)
            columns = [columns[index] for index in order]
            replica_energy = np.asarray(energies, dtype=float)[order]
            replica_ids = replicas["replica_id"].to_numpy(dtype=int)
            if len(np.unique(replica_ids)) != len(replica_ids):
                raise ValueError(f"{channel.channel_id}: duplicate replica identifiers")
            if common_replica_ids is None:
                common_replica_ids = replica_ids
            elif not np.array_equal(common_replica_ids, replica_ids):
                raise ValueError("replica identifiers differ between channels")

            nominal_energy, nominal_sigma = _validated_grid(
                curve["energy_MeV"], curve["sigma_nominal_mb"],
                f"{channel.channel_id} nominal")
            replica_energy, replica_sigma = _validated_grid(
                replica_energy, replicas[columns].to_numpy(dtype=float),
                f"{channel.channel_id} replicas")
            curves[channel.channel_id] = ChannelCurves(
                threshold_MeV=float(thresholds[channel.channel_id]),
                nominal_energy_MeV=nominal_energy,
                nominal_sigma_mb=nominal_sigma,
                replica_energy_MeV=replica_energy,
                replica_sigma_mb=replica_sigma,
            )

        if common_replica_ids is None:
            raise ValueError("fit directory contains no channel replicas")
        expected = int(metadata["n_replicas"])
        if len(common_replica_ids) != expected:
            raise ValueError(
                f"fit metadata specifies {expected} replicas; found {len(common_replica_ids)}")
        return cls(curves, common_replica_ids, fit_directory.resolve())

    def nominal(self, channel_id: str, energies) -> np.ndarray:
        curve = self.channels[channel_id]
        return interpolate_curves(
            curve.nominal_energy_MeV,
            curve.nominal_sigma_mb,
            energies,
            threshold_MeV=curve.threshold_MeV,
            label=f"{channel_id} nominal",
        )

    def replicas(self, channel_id: str, energies) -> np.ndarray:
        curve = self.channels[channel_id]
        return interpolate_curves(
            curve.replica_energy_MeV,
            curve.replica_sigma_mb,
            energies,
            threshold_MeV=curve.threshold_MeV,
            label=f"{channel_id} replicas",
        )


def validate_exposure_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize a target-resolved proton exposure table."""

    missing = [field for field in EXPOSURE_FIELDS if field not in frame.columns]
    if missing:
        raise ValueError(f"exposure table missing columns: {', '.join(missing)}")
    table = frame.loc[:, EXPOSURE_FIELDS].copy()
    for field in NUMERIC_EXPOSURE_FIELDS:
        table[field] = pd.to_numeric(table[field], errors="raise")
    numeric = table.loc[:, NUMERIC_EXPOSURE_FIELDS].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("exposure table contains non-finite values")

    known_targets = {channel.target for channel in CHANNELS}
    unknown = sorted(set(table["target"]) - known_targets)
    if unknown:
        raise ValueError(f"exposure table contains unknown targets: {', '.join(unknown)}")
    if np.any(table["energy_low_MeV"] >= table["energy_high_MeV"]):
        raise ValueError("energy-bin lower edges must be below upper edges")
    if np.any(table["depth_low_mm"] >= table["depth_high_mm"]):
        raise ValueError("depth-bin lower edges must be below upper edges")
    if np.any(
        (table["energy_mean_MeV"] < table["energy_low_MeV"])
        | (table["energy_mean_MeV"] > table["energy_high_MeV"])
    ):
        raise ValueError("mean proton energy must lie inside its bin")
    if np.any(
        (table["depth_mean_mm"] < table["depth_low_mm"])
        | (table["depth_mean_mm"] > table["depth_high_mm"])
    ):
        raise ValueError("mean depth must lie inside its bin")
    if np.any(table["target_exposure_cm2_inv"] < 0):
        raise ValueError("target exposure must be non-negative")

    key = [
        "target", "energy_low_MeV", "energy_high_MeV",
        "depth_low_mm", "depth_high_mm",
    ]
    if table.duplicated(key).any():
        raise ValueError("exposure table contains duplicate target-energy-depth bins")

    intervals = (
        table[["depth_low_mm", "depth_high_mm"]]
        .drop_duplicates()
        .sort_values(["depth_low_mm", "depth_high_mm"])
        .to_numpy(dtype=float)
    )
    if len(intervals) > 1 and np.any(intervals[1:, 0] < intervals[:-1, 1] - 1.0e-12):
        raise ValueError("depth bins overlap")
    return table.sort_values(
        ["depth_low_mm", "target", "energy_low_MeV"]
    ).reset_index(drop=True)


def load_exposure_table(path: Path | str) -> pd.DataFrame:
    return validate_exposure_table(pd.read_csv(path))


def _bin_indices(values: np.ndarray, edges: np.ndarray, label: str) -> np.ndarray:
    if edges.ndim != 1 or len(edges) < 2:
        raise ValueError(f"{label} edges must contain at least two values")
    if not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError(f"{label} edges must be finite and increasing")
    if np.any(values < edges[0]) or np.any(values > edges[-1]):
        raise ValueError(f"one or more {label} values lie outside the bin edges")
    indices = np.searchsorted(edges, values, side="right") - 1
    indices[values == edges[-1]] = len(edges) - 2
    return indices


def accumulate_step_exposure(
    steps: pd.DataFrame,
    energy_edges_MeV,
    depth_edges_mm,
) -> pd.DataFrame:
    """Reference accumulation of synthetic proton steps into exposure bins.

    The production Geant4 scorer will perform this accumulation during
    transport.  This Python implementation supplies an independently testable
    definition of the output table.
    """

    missing = [field for field in STEP_FIELDS if field not in steps.columns]
    if missing:
        raise ValueError(f"step table missing columns: {', '.join(missing)}")
    frame = steps.loc[:, STEP_FIELDS].copy()
    numeric_fields = STEP_FIELDS[1:]
    for field in numeric_fields:
        frame[field] = pd.to_numeric(frame[field], errors="raise")
    numeric = frame.loc[:, numeric_fields].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("step table contains non-finite values")
    if np.any(frame[[
        "proton_weight", "target_number_density_cm3", "step_length_cm",
    ]].to_numpy(dtype=float) < 0):
        raise ValueError("step weights, number densities and lengths must be non-negative")
    known_targets = {channel.target for channel in CHANNELS}
    unknown = sorted(set(frame["target"]) - known_targets)
    if unknown:
        raise ValueError(f"step table contains unknown targets: {', '.join(unknown)}")

    energy_edges = np.asarray(energy_edges_MeV, dtype=float)
    depth_edges = np.asarray(depth_edges_mm, dtype=float)
    frame["energy_index"] = _bin_indices(
        frame["energy_MeV"].to_numpy(dtype=float), energy_edges, "energy")
    frame["depth_index"] = _bin_indices(
        frame["depth_mm"].to_numpy(dtype=float), depth_edges, "depth")
    frame["_exposure"] = (
        frame["proton_weight"]
        * frame["target_number_density_cm3"]
        * frame["step_length_cm"]
    )
    frame = frame.loc[frame["_exposure"] > 0].copy()
    if frame.empty:
        return pd.DataFrame(columns=EXPOSURE_FIELDS)
    frame["_energy_moment"] = frame["_exposure"] * frame["energy_MeV"]
    frame["_depth_moment"] = frame["_exposure"] * frame["depth_mm"]
    grouped = (
        frame.groupby(["target", "energy_index", "depth_index"], sort=True)
        .agg(
            target_exposure_cm2_inv=("_exposure", "sum"),
            energy_moment=("_energy_moment", "sum"),
            depth_moment=("_depth_moment", "sum"),
        )
        .reset_index()
    )
    energy_index = grouped["energy_index"].to_numpy(dtype=int)
    depth_index = grouped["depth_index"].to_numpy(dtype=int)
    exposure = grouped["target_exposure_cm2_inv"].to_numpy(dtype=float)
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


def distal_r50(depth_mm, production) -> float:
    """Return the distal half-maximum crossing using linear interpolation."""

    depth = np.asarray(depth_mm, dtype=float)
    values = np.asarray(production, dtype=float)
    if depth.ndim != 1 or values.shape != depth.shape:
        raise ValueError("R50 requires equal one-dimensional arrays")
    if len(depth) < 2:
        return float("nan")
    if not np.isfinite(depth).all() or np.any(np.diff(depth) <= 0):
        raise ValueError("R50 depths must be finite and increasing")
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("R50 production must be finite and non-negative")
    maximum = float(np.max(values))
    if maximum <= 0:
        return float("nan")
    half = 0.5 * maximum
    above = np.flatnonzero(values >= half)
    last = int(above[-1])
    if last == len(values) - 1:
        return float("nan")
    x0, x1 = depth[last], depth[last + 1]
    y0, y1 = values[last], values[last + 1]
    if y1 == y0:
        return float(x0)
    return float(x0 + (half - y0) * (x1 - x0) / (y1 - y0))


@dataclass(frozen=True)
class FoldingResult:
    nominal_channel_contributions: pd.DataFrame
    nominal_isotope_profiles: pd.DataFrame
    replica_isotope_profiles: pd.DataFrame
    production_summary: pd.DataFrame


def _depth_bins(exposure: pd.DataFrame) -> pd.DataFrame:
    bins = (
        exposure[["depth_low_mm", "depth_high_mm"]]
        .drop_duplicates()
        .sort_values(["depth_low_mm", "depth_high_mm"])
        .reset_index(drop=True)
    )
    bins["depth_mm"] = 0.5 * (bins["depth_low_mm"] + bins["depth_high_mm"])
    bins["depth_index"] = np.arange(len(bins), dtype=int)
    return bins


def fold_exposure(
    exposure: pd.DataFrame,
    cross_sections: CrossSectionEnsemble,
) -> FoldingResult:
    """Fold nominal and replica cross sections with a validated exposure table."""

    table = validate_exposure_table(exposure)
    depth_bins = _depth_bins(table)
    depth_key = {
        (row.depth_low_mm, row.depth_high_mm): int(row.depth_index)
        for row in depth_bins.itertuples(index=False)
    }
    n_depth = len(depth_bins)
    n_replicas = len(cross_sections.replica_ids)
    isotope_index = {isotope: index for index, isotope in enumerate(ISOTOPES)}

    nominal_array = np.zeros((len(ISOTOPES), n_depth), dtype=float)
    replica_array = np.zeros((n_replicas, len(ISOTOPES), n_depth), dtype=float)
    channel_frames = []

    for channel in CHANNELS:
        selected = table.loc[table["target"] == channel.target].copy()
        if selected.empty:
            continue
        energy = selected["energy_mean_MeV"].to_numpy(dtype=float)
        exposure_values = selected["target_exposure_cm2_inv"].to_numpy(dtype=float)
        nominal_sigma = cross_sections.nominal(channel.channel_id, energy)
        replica_sigma = cross_sections.replicas(channel.channel_id, energy)
        nominal_yield = exposure_values * nominal_sigma * MB_TO_CM2
        replica_yield = replica_sigma * exposure_values[np.newaxis, :] * MB_TO_CM2

        selected["channel_id"] = channel.channel_id
        selected["residual"] = channel.residual
        selected["sigma_nominal_mb"] = nominal_sigma
        selected["expected_nuclei"] = nominal_yield
        channel_frames.append(selected)

        isotope = isotope_index[channel.residual]
        for column, row in enumerate(selected.itertuples(index=False)):
            depth = depth_key[(row.depth_low_mm, row.depth_high_mm)]
            nominal_array[isotope, depth] += nominal_yield[column]
            replica_array[:, isotope, depth] += replica_yield[:, column]

    if not channel_frames:
        raise ValueError("exposure table has no rows for the configured channels")
    nominal_channels = pd.concat(channel_frames, ignore_index=True)

    nominal_profile_rows = []
    replica_profile_rows = []
    depth_records = list(depth_bins.itertuples(index=False))
    profile_labels = ISOTOPES + ("all",)
    nominal_profiles = np.concatenate(
        [nominal_array, nominal_array.sum(axis=0, keepdims=True)], axis=0)
    replica_profiles = np.concatenate(
        [replica_array, replica_array.sum(axis=1, keepdims=True)], axis=1)

    for isotope_position, isotope in enumerate(profile_labels):
        for depth_position, depth in enumerate(depth_records):
            nominal_profile_rows.append({
                "isotope": isotope,
                "depth_low_mm": depth.depth_low_mm,
                "depth_high_mm": depth.depth_high_mm,
                "depth_mm": depth.depth_mm,
                "expected_nuclei": nominal_profiles[isotope_position, depth_position],
            })
        for replica_position, replica_id in enumerate(cross_sections.replica_ids):
            for depth_position, depth in enumerate(depth_records):
                replica_profile_rows.append({
                    "replica_id": int(replica_id),
                    "isotope": isotope,
                    "depth_low_mm": depth.depth_low_mm,
                    "depth_high_mm": depth.depth_high_mm,
                    "depth_mm": depth.depth_mm,
                    "expected_nuclei": replica_profiles[
                        replica_position, isotope_position, depth_position],
                })

    nominal_profiles_frame = pd.DataFrame(nominal_profile_rows)
    replica_profiles_frame = pd.DataFrame(replica_profile_rows)
    depth_centres = depth_bins["depth_mm"].to_numpy(dtype=float)

    summary_rows = []
    for isotope_position, isotope in enumerate(profile_labels):
        nominal_values = nominal_profiles[isotope_position]
        nominal_r50 = distal_r50(depth_centres, nominal_values)
        summary_rows.append({
            "model": "nominal",
            "replica_id": "",
            "isotope": isotope,
            "expected_nuclei": float(np.sum(nominal_values)),
            "R50_prod_mm": nominal_r50,
            "R50_shift_mm": 0.0,
        })
        for replica_position, replica_id in enumerate(cross_sections.replica_ids):
            values = replica_profiles[replica_position, isotope_position]
            r50 = distal_r50(depth_centres, values)
            summary_rows.append({
                "model": "replica",
                "replica_id": int(replica_id),
                "isotope": isotope,
                "expected_nuclei": float(np.sum(values)),
                "R50_prod_mm": r50,
                "R50_shift_mm": r50 - nominal_r50,
            })

    return FoldingResult(
        nominal_channel_contributions=nominal_channels,
        nominal_isotope_profiles=nominal_profiles_frame,
        replica_isotope_profiles=replica_profiles_frame,
        production_summary=pd.DataFrame(summary_rows),
    )


def write_result(
    result: FoldingResult,
    output_directory: Path | str,
    *,
    exposure_path: Path | None = None,
    cross_sections: CrossSectionEnsemble | None = None,
) -> None:
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    products = {
        "nominal_channel_contributions.csv": result.nominal_channel_contributions,
        "nominal_isotope_profiles.csv": result.nominal_isotope_profiles,
        "replica_isotope_profiles.csv": result.replica_isotope_profiles,
        "production_summary.csv": result.production_summary,
    }
    for name, frame in products.items():
        frame.to_csv(output_directory / name, index=False)

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "mb_to_cm2": MB_TO_CM2,
        "n_replicas": (
            len(cross_sections.replica_ids) if cross_sections is not None else None),
        "exposure_file": str(exposure_path.resolve()) if exposure_path else None,
        "exposure_sha256": sha256(exposure_path) if exposure_path else None,
        "fit_directory": (
            str(cross_sections.fit_directory)
            if cross_sections is not None and cross_sections.fit_directory else None),
    }
    if cross_sections is not None and cross_sections.fit_directory is not None:
        fit_meta = cross_sections.fit_directory / "fit_meta.json"
        metadata["fit_meta_sha256"] = sha256(fit_meta)
        metadata["fit_files_sha256"] = {
            path.name: sha256(path)
            for channel in CHANNELS
            for path in (
                cross_sections.fit_directory / f"{channel.channel_id}_curve.csv",
                cross_sections.fit_directory / f"{channel.channel_id}_replicas.csv",
            )
        }
    (output_directory / "folding_meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exposure_csv", type=Path)
    parser.add_argument(
        "--fit-dir", type=Path, default=Path("data/xsections/fits"),
        help="directory containing fit_meta.json and channel fit products")
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()

    exposure = load_exposure_table(arguments.exposure_csv)
    cross_sections = CrossSectionEnsemble.from_fit_directory(arguments.fit_dir)
    result = fold_exposure(exposure, cross_sections)
    write_result(
        result,
        arguments.output_dir,
        exposure_path=arguments.exposure_csv,
        cross_sections=cross_sections,
    )
    print(
        f"folded {len(exposure)} exposure bins with "
        f"{len(cross_sections.replica_ids)} cross-section replicas")


if __name__ == "__main__":
    main()
