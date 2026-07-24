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
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .channels import CHANNELS
from common.isotopes import ISOTOPES as ISOTOPE_DATA, NAME_TO_ID
from decay_sampling.scenarios import (
    DEFAULT_SCENARIO_CONFIG,
    HandoffScenario,
    resolve_scenario,
)


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
EXPOSURE_META_FIELDS = (
    "schema_version",
    "run_id",
    "exposure_file",
    "exposure_sha256",
    "n_protons",
    "target_dose_Gy",
    "Np_per_Gy",
    "physics_list",
    "geant4_version",
    "random_seed",
    "software_revision",
    "beam_axis",
    "depth_origin",
    "depth_unit",
    "energy_edges_MeV",
    "depth_edges_mm",
)


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


def validate_threshold_coverage(grid, threshold_MeV: float, label: str) -> None:
    """Require a fitted grid to bracket its physical production threshold."""

    energy = np.asarray(grid, dtype=float)
    threshold = float(threshold_MeV)
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError(f"{label}: reaction threshold must be finite and non-negative")
    if energy[0] > threshold + 1.0e-12:
        raise ValueError(
            f"{label}: fitted grid starts at {energy[0]:g} MeV, above the "
            f"{threshold:g} MeV reaction threshold")
    if energy[-1] <= threshold:
        raise ValueError(
            f"{label}: fitted grid ends at or below the reaction threshold")


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
            threshold = float(thresholds[channel.channel_id])
            validate_threshold_coverage(
                nominal_energy, threshold, f"{channel.channel_id} nominal")
            validate_threshold_coverage(
                replica_energy, threshold, f"{channel.channel_id} replicas")
            curves[channel.channel_id] = ChannelCurves(
                threshold_MeV=threshold,
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


@dataclass(frozen=True)
class ExposureMetadata:
    schema_version: int
    run_id: str
    exposure_file: str
    exposure_sha256: str
    n_protons: float
    target_dose_Gy: float
    Np_per_Gy: float
    physics_list: str
    geant4_version: str
    random_seed: int
    software_revision: str
    beam_axis: str
    depth_origin: str
    depth_unit: str
    energy_edges_MeV: np.ndarray
    depth_edges_mm: np.ndarray
    source_path: Path | None = None

    @classmethod
    def synthetic(cls) -> "ExposureMetadata":
        """Unit normalization for analytic tests with no physical run metadata."""

        return cls(
            schema_version=SCHEMA_VERSION,
            run_id="synthetic",
            exposure_file="synthetic.csv",
            exposure_sha256="synthetic",
            n_protons=1.0,
            target_dose_Gy=1.0,
            Np_per_Gy=1.0,
            physics_list="synthetic",
            geant4_version="synthetic",
            random_seed=0,
            software_revision="synthetic",
            beam_axis="z",
            depth_origin="synthetic",
            depth_unit="mm",
            energy_edges_MeV=np.asarray([], dtype=float),
            depth_edges_mm=np.asarray([], dtype=float),
        )


def _validated_edges(values, label: str, *, allow_empty: bool = False) -> np.ndarray:
    edges = np.asarray(values, dtype=float)
    if allow_empty and edges.size == 0:
        return edges
    if edges.ndim != 1 or len(edges) < 2:
        raise ValueError(f"{label} must contain at least two edges")
    if not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError(f"{label} must be finite and strictly increasing")
    return edges


def load_exposure_metadata(
    path: Path | str,
    exposure_path: Path | str | None = None,
) -> ExposureMetadata:
    """Read and validate the companion provenance and normalization record."""

    source = Path(path).resolve()
    document = json.loads(source.read_text(encoding="utf-8"))
    missing = [field for field in EXPOSURE_META_FIELDS if field not in document]
    if missing:
        raise ValueError(f"exposure metadata missing fields: {', '.join(missing)}")
    if int(document["schema_version"]) != SCHEMA_VERSION:
        raise ValueError("unsupported exposure-metadata schema version")
    for field in (
        "run_id", "exposure_file", "exposure_sha256", "physics_list",
        "geant4_version", "software_revision", "beam_axis", "depth_origin",
        "depth_unit",
    ):
        if not isinstance(document[field], str) or not document[field].strip():
            raise ValueError(f"exposure metadata field '{field}' must be non-empty text")
    if document["depth_unit"] != "mm":
        raise ValueError("exposure metadata depth_unit must be 'mm'")
    n_protons = float(document["n_protons"])
    target_dose = float(document["target_dose_Gy"])
    np_per_gy = float(document["Np_per_Gy"])
    if not all(math.isfinite(value) and value > 0 for value in (
        n_protons, target_dose, np_per_gy,
    )):
        raise ValueError("n_protons, target_dose_Gy and Np_per_Gy must be positive")
    expected_np_per_gy = n_protons / target_dose
    if not math.isclose(np_per_gy, expected_np_per_gy, rel_tol=1.0e-10):
        raise ValueError("Np_per_Gy is inconsistent with n_protons/target_dose_Gy")
    try:
        random_seed = int(document["random_seed"])
    except (TypeError, ValueError) as error:
        raise ValueError("random_seed must be an integer") from error

    energy_edges = _validated_edges(document["energy_edges_MeV"], "energy edges")
    depth_edges = _validated_edges(document["depth_edges_mm"], "depth edges")
    if exposure_path is not None:
        exposure = Path(exposure_path).resolve()
        if Path(document["exposure_file"]).name != exposure.name:
            raise ValueError("exposure metadata names a different exposure file")
        if document["exposure_sha256"].lower() != sha256(exposure):
            raise ValueError("exposure file SHA-256 does not match its metadata")

    return ExposureMetadata(
        schema_version=SCHEMA_VERSION,
        run_id=document["run_id"],
        exposure_file=document["exposure_file"],
        exposure_sha256=document["exposure_sha256"].lower(),
        n_protons=n_protons,
        target_dose_Gy=target_dose,
        Np_per_Gy=np_per_gy,
        physics_list=document["physics_list"],
        geant4_version=document["geant4_version"],
        random_seed=random_seed,
        software_revision=document["software_revision"],
        beam_axis=document["beam_axis"],
        depth_origin=document["depth_origin"],
        depth_unit=document["depth_unit"],
        energy_edges_MeV=energy_edges,
        depth_edges_mm=depth_edges,
        source_path=source,
    )


def validate_exposure_against_metadata(
    exposure: pd.DataFrame,
    metadata: ExposureMetadata,
) -> None:
    """Check that every recorded bin is drawn from the declared run grids."""

    if metadata.energy_edges_MeV.size == 0 and metadata.depth_edges_mm.size == 0:
        return
    for low_field, high_field, edges, label in (
        ("energy_low_MeV", "energy_high_MeV", metadata.energy_edges_MeV, "energy"),
        ("depth_low_mm", "depth_high_mm", metadata.depth_edges_mm, "depth"),
    ):
        for field in (low_field, high_field):
            values = exposure[field].to_numpy(dtype=float)
            represented = np.any(
                np.isclose(values[:, np.newaxis], edges[np.newaxis, :], rtol=0, atol=1e-10),
                axis=1,
            )
            if not np.all(represented):
                raise ValueError(f"exposure {label} bins do not belong to metadata grid")


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

    for group_key, group in table.groupby(
        ["target", "depth_low_mm", "depth_high_mm"], sort=False
    ):
        energy_intervals = (
            group[["energy_low_MeV", "energy_high_MeV"]]
            .drop_duplicates()
            .sort_values(["energy_low_MeV", "energy_high_MeV"])
            .to_numpy(dtype=float)
        )
        if (
            len(energy_intervals) > 1
            and np.any(energy_intervals[1:, 0] < energy_intervals[:-1, 1] - 1.0e-12)
        ):
            target, depth_low, depth_high = group_key
            raise ValueError(
                "energy bins overlap for "
                f"{target} at depth [{depth_low:g}, {depth_high:g}] mm")

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
    profile_bands: pd.DataFrame
    production_summary: pd.DataFrame
    uncertainty_summary: pd.DataFrame
    profile_labels: tuple[str, ...]
    depth_mm: np.ndarray
    nominal_profile_values: np.ndarray
    replica_profile_values: np.ndarray


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
    metadata: ExposureMetadata | None = None,
    scenario: HandoffScenario | None = None,
) -> FoldingResult:
    """Fold nominal and replica cross sections with a validated exposure table."""

    table = validate_exposure_table(exposure)
    metadata = metadata or ExposureMetadata.synthetic()
    scenario = scenario or resolve_scenario("inroom")
    validate_exposure_against_metadata(table, metadata)
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
        selected["expected_nuclei_run"] = nominal_yield
        selected["expected_nuclei_per_proton"] = nominal_yield / metadata.n_protons
        selected["expected_nuclei_per_Gy"] = nominal_yield / metadata.target_dose_Gy
        channel_frames.append(selected)

        isotope = isotope_index[channel.residual]
        depth_indices = np.fromiter(
            (
                depth_key[(row.depth_low_mm, row.depth_high_mm)]
                for row in selected.itertuples(index=False)
            ),
            dtype=int,
            count=len(selected),
        )
        assignment = np.zeros((len(selected), n_depth), dtype=float)
        assignment[np.arange(len(selected)), depth_indices] = 1.0
        nominal_array[isotope] += nominal_yield @ assignment
        replica_array[:, isotope, :] += replica_yield @ assignment

    if not channel_frames:
        raise ValueError("exposure table has no rows for the configured channels")
    nominal_channels = pd.concat(channel_frames, ignore_index=True)

    activity_factors = np.asarray([
        scenario.measured_fraction(ISOTOPE_DATA[NAME_TO_ID[name]].lam)
        for name in ISOTOPES
    ])
    profile_labels = ISOTOPES + ("all_production", f"all_{scenario.name}")
    nominal_profiles = np.concatenate([
        nominal_array,
        nominal_array.sum(axis=0, keepdims=True),
        (nominal_array * activity_factors[:, np.newaxis]).sum(
            axis=0, keepdims=True),
    ], axis=0)
    replica_profiles = np.concatenate([
        replica_array,
        replica_array.sum(axis=1, keepdims=True),
        (replica_array * activity_factors[np.newaxis, :, np.newaxis]).sum(
            axis=1, keepdims=True),
    ], axis=1)

    label_count = len(profile_labels)
    label_column = np.repeat(np.asarray(profile_labels), n_depth)
    depth_low = np.tile(depth_bins["depth_low_mm"].to_numpy(dtype=float), label_count)
    depth_high = np.tile(depth_bins["depth_high_mm"].to_numpy(dtype=float), label_count)
    depth_centres = depth_bins["depth_mm"].to_numpy(dtype=float)
    depth_column = np.tile(depth_centres, label_count)
    nominal_run = nominal_profiles.reshape(-1)
    nominal_profiles_frame = pd.DataFrame({
        "profile_label": label_column,
        "quantity": np.where(
            np.char.startswith(label_column.astype(str), "all_"),
            np.where(label_column == "all_production", "production_nuclei", "measured_decays"),
            "production_nuclei",
        ),
        "depth_low_mm": depth_low,
        "depth_high_mm": depth_high,
        "depth_mm": depth_column,
        "expected_count_run": nominal_run,
        "expected_count_per_proton": nominal_run / metadata.n_protons,
        "expected_count_per_Gy": nominal_run / metadata.target_dose_Gy,
    })

    replica_count = n_replicas * label_count * n_depth
    replica_run = replica_profiles.reshape(-1)
    replica_profiles_frame = pd.DataFrame({
        "replica_id": np.repeat(
            cross_sections.replica_ids, label_count * n_depth),
        "profile_label": np.tile(label_column, n_replicas),
        "quantity": np.tile(
            nominal_profiles_frame["quantity"].to_numpy(), n_replicas),
        "depth_low_mm": np.tile(depth_low, n_replicas),
        "depth_high_mm": np.tile(depth_high, n_replicas),
        "depth_mm": np.tile(depth_column, n_replicas),
        "expected_count_run": replica_run,
        "expected_count_per_proton": replica_run / metadata.n_protons,
        "expected_count_per_Gy": replica_run / metadata.target_dose_Gy,
    })
    if len(replica_profiles_frame) != replica_count:
        raise RuntimeError("internal replica-profile shape mismatch")

    quantiles = np.quantile(replica_profiles, [0.16, 0.50, 0.84], axis=0)
    band_rows = []
    for label_position, label in enumerate(profile_labels):
        for depth_position in range(n_depth):
            row = {
                "profile_label": label,
                "quantity": (
                    "measured_decays"
                    if label.startswith("all_") and label != "all_production"
                    else "production_nuclei"),
                "depth_low_mm": depth_bins.loc[depth_position, "depth_low_mm"],
                "depth_high_mm": depth_bins.loc[depth_position, "depth_high_mm"],
                "depth_mm": depth_centres[depth_position],
                "nominal_run": nominal_profiles[label_position, depth_position],
                "nominal_per_proton": (
                    nominal_profiles[label_position, depth_position]
                    / metadata.n_protons),
                "nominal_per_Gy": (
                    nominal_profiles[label_position, depth_position]
                    / metadata.target_dose_Gy),
            }
            for quantile_position, quantile_name in enumerate(("q16", "q50", "q84")):
                value = quantiles[quantile_position, label_position, depth_position]
                row[f"{quantile_name}_run"] = value
                row[f"{quantile_name}_per_proton"] = value / metadata.n_protons
                row[f"{quantile_name}_per_Gy"] = value / metadata.target_dose_Gy
            band_rows.append(row)
    profile_bands = pd.DataFrame(band_rows)

    summary_rows = []
    uncertainty_rows = []
    for label_position, label in enumerate(profile_labels):
        nominal_values = nominal_profiles[label_position]
        nominal_r50 = distal_r50(depth_centres, nominal_values)
        nominal_yield = float(np.sum(nominal_values))
        summary_rows.append({
            "model": "nominal",
            "replica_id": "",
            "profile_label": label,
            "expected_count_run": nominal_yield,
            "expected_count_per_proton": nominal_yield / metadata.n_protons,
            "expected_count_per_Gy": nominal_yield / metadata.target_dose_Gy,
            "R50_prod_mm": nominal_r50,
            "R50_shift_mm": 0.0,
        })
        replica_yields = np.sum(replica_profiles[:, label_position, :], axis=1)
        replica_r50 = np.asarray([
            distal_r50(depth_centres, values)
            for values in replica_profiles[:, label_position, :]
        ])
        replica_shifts = replica_r50 - nominal_r50
        for replica_position, replica_id in enumerate(cross_sections.replica_ids):
            replica_yield = float(replica_yields[replica_position])
            summary_rows.append({
                "model": "replica",
                "replica_id": int(replica_id),
                "profile_label": label,
                "expected_count_run": replica_yield,
                "expected_count_per_proton": replica_yield / metadata.n_protons,
                "expected_count_per_Gy": replica_yield / metadata.target_dose_Gy,
                "R50_prod_mm": replica_r50[replica_position],
                "R50_shift_mm": replica_shifts[replica_position],
            })
        yield_q = np.quantile(replica_yields, [0.16, 0.50, 0.84])
        finite_shifts = replica_shifts[np.isfinite(replica_shifts)]
        shift_q = (
            np.quantile(finite_shifts, [0.16, 0.50, 0.84])
            if len(finite_shifts) else np.full(3, np.nan)
        )
        uncertainty_rows.append({
            "profile_label": label,
            "nominal_yield_run": nominal_yield,
            "nominal_yield_per_proton": nominal_yield / metadata.n_protons,
            "nominal_yield_per_Gy": nominal_yield / metadata.target_dose_Gy,
            "yield_q16_run": yield_q[0],
            "yield_q50_run": yield_q[1],
            "yield_q84_run": yield_q[2],
            "yield_half_width_run": 0.5 * (yield_q[2] - yield_q[0]),
            "yield_q16_per_proton": yield_q[0] / metadata.n_protons,
            "yield_q50_per_proton": yield_q[1] / metadata.n_protons,
            "yield_q84_per_proton": yield_q[2] / metadata.n_protons,
            "yield_half_width_per_proton": (
                0.5 * (yield_q[2] - yield_q[0]) / metadata.n_protons),
            "yield_q16_per_Gy": yield_q[0] / metadata.target_dose_Gy,
            "yield_q50_per_Gy": yield_q[1] / metadata.target_dose_Gy,
            "yield_q84_per_Gy": yield_q[2] / metadata.target_dose_Gy,
            "yield_half_width_per_Gy": (
                0.5 * (yield_q[2] - yield_q[0]) / metadata.target_dose_Gy),
            "R50_nominal_mm": nominal_r50,
            "R50_shift_q16_mm": shift_q[0],
            "R50_shift_q50_mm": shift_q[1],
            "R50_shift_q84_mm": shift_q[2],
            "R50_shift_half_width_mm": 0.5 * (shift_q[2] - shift_q[0]),
        })

    return FoldingResult(
        nominal_channel_contributions=nominal_channels,
        nominal_isotope_profiles=nominal_profiles_frame,
        replica_isotope_profiles=replica_profiles_frame,
        profile_bands=profile_bands,
        production_summary=pd.DataFrame(summary_rows),
        uncertainty_summary=pd.DataFrame(uncertainty_rows),
        profile_labels=profile_labels,
        depth_mm=depth_centres,
        nominal_profile_values=nominal_profiles,
        replica_profile_values=replica_profiles,
    )


def write_result(
    result: FoldingResult,
    output_directory: Path | str,
    *,
    exposure_path: Path | None = None,
    exposure_metadata: ExposureMetadata | None = None,
    cross_sections: CrossSectionEnsemble | None = None,
    scenario: HandoffScenario | None = None,
    native_route_summary_path: Path | None = None,
    write_replica_profiles: bool = True,
) -> None:
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    products = {
        "nominal_channel_contributions.csv": result.nominal_channel_contributions,
        "nominal_isotope_profiles.csv": result.nominal_isotope_profiles,
        "profile_bands.csv": result.profile_bands,
        "production_summary.csv": result.production_summary,
        "uncertainty_summary.csv": result.uncertainty_summary,
    }
    for name, frame in products.items():
        frame.to_csv(output_directory / name, index=False)
    if write_replica_profiles:
        result.replica_isotope_profiles.to_csv(
            output_directory / "replica_isotope_profiles.csv.gz",
            index=False,
            compression="gzip",
        )

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
        "output_normalization": ["per_run", "per_proton", "per_Gy"],
        "replica_profiles_written": write_replica_profiles,
    }
    if exposure_metadata is not None:
        metadata["exposure_metadata"] = {
            field: (
                getattr(exposure_metadata, field).tolist()
                if isinstance(getattr(exposure_metadata, field), np.ndarray)
                else getattr(exposure_metadata, field)
            )
            for field in EXPOSURE_META_FIELDS
        }
        metadata["exposure_meta_file"] = (
            str(exposure_metadata.source_path)
            if exposure_metadata.source_path is not None else None)
        metadata["exposure_meta_sha256"] = (
            sha256(exposure_metadata.source_path)
            if exposure_metadata.source_path is not None else None)
    if scenario is not None:
        metadata["handoff_scenario"] = {
            "name": scenario.name,
            "description": scenario.description,
            "t_irr_s": scenario.t_irr_s,
            "t_del_s": scenario.t_del_s,
            "t_meas_s": scenario.t_meas_s,
            "config_path": str(scenario.config_path),
            "config_sha256": scenario.config_sha256,
            "measured_fraction_by_isotope": {
                name: scenario.measured_fraction(ISOTOPE_DATA[NAME_TO_ID[name]].lam)
                for name in ISOTOPES
            },
        }
    if native_route_summary_path is not None:
        route_path = Path(native_route_summary_path).resolve()
        route_summary = pd.read_csv(route_path)
        required = {
            "profile_label", "modeled_count", "unmodeled_count",
            "total_count", "unmodeled_fraction",
        }
        missing = sorted(required - set(route_summary.columns))
        if missing:
            raise ValueError(
                "native-route summary missing columns: " + ", ".join(missing))
        required_labels = {"all_production"}
        if scenario is not None:
            required_labels.add(f"all_{scenario.name}")
        if not required_labels.issubset(set(route_summary["profile_label"])):
            raise ValueError("native-route summary lacks required profile labels")
        metadata["native_route_diagnostic"] = {
            "source_file": str(route_path),
            "source_sha256": sha256(route_path),
            "fractions": route_summary.loc[
                route_summary["profile_label"].isin(required_labels),
                ["profile_label", "unmodeled_fraction"],
            ].to_dict(orient="records"),
            "affects_folded_source": False,
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
    parser.add_argument("--exposure-meta", type=Path, required=True)
    parser.add_argument(
        "--fit-dir", type=Path, default=Path("data/xsections/fits"),
        help="directory containing fit_meta.json and channel fit products")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--scenario", default="inroom")
    parser.add_argument(
        "--scenario-config", type=Path, default=DEFAULT_SCENARIO_CONFIG)
    parser.add_argument(
        "--no-replica-profiles", action="store_true",
        help="omit the large compressed per-replica profile table")
    parser.add_argument(
        "--native-route-summary", type=Path,
        help="optional diagnostic summary to record in folding provenance")
    arguments = parser.parse_args()

    exposure = load_exposure_table(arguments.exposure_csv)
    exposure_metadata = load_exposure_metadata(
        arguments.exposure_meta, arguments.exposure_csv)
    scenario = resolve_scenario(arguments.scenario, arguments.scenario_config)
    cross_sections = CrossSectionEnsemble.from_fit_directory(arguments.fit_dir)
    result = fold_exposure(exposure, cross_sections, exposure_metadata, scenario)
    write_result(
        result,
        arguments.output_dir,
        exposure_path=arguments.exposure_csv,
        exposure_metadata=exposure_metadata,
        cross_sections=cross_sections,
        scenario=scenario,
        native_route_summary_path=arguments.native_route_summary,
        write_replica_profiles=not arguments.no_replica_profiles,
    )
    print(
        f"folded {len(exposure)} exposure bins with "
        f"{len(cross_sections.replica_ids)} cross-section replicas")


if __name__ == "__main__":
    main()
