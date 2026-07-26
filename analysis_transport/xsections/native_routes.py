#!/usr/bin/env python3
"""Classify native Geant4 emitters outside the five-channel folded source."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common.isotopes import ISOTOPES, NAME_TO_ID
from decay_sampling.scenarios import DEFAULT_SCENARIO_CONFIG, resolve_scenario

from .channels import CHANNELS


NATIVE_ROUTE_FIELDS = (
    "projectile",
    "target",
    "residual",
    "depth_low_mm",
    "depth_high_mm",
    "production_count",
)
SELECTED_ROUTES = {(channel.target, channel.residual) for channel in CHANNELS}


def validate_native_routes(frame: pd.DataFrame) -> pd.DataFrame:
    missing = [field for field in NATIVE_ROUTE_FIELDS if field not in frame.columns]
    if missing:
        raise ValueError(f"native-route table missing columns: {', '.join(missing)}")
    table = frame.loc[:, NATIVE_ROUTE_FIELDS].copy()
    for field in ("depth_low_mm", "depth_high_mm", "production_count"):
        table[field] = pd.to_numeric(table[field], errors="raise")
    numeric = table[[
        "depth_low_mm", "depth_high_mm", "production_count",
    ]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("native-route table contains non-finite values")
    if np.any(table["depth_low_mm"] >= table["depth_high_mm"]):
        raise ValueError("native-route depth-bin lower edges must be below upper edges")
    if np.any(table["production_count"] < 0):
        raise ValueError("native-route production counts must be non-negative")
    unknown = sorted(set(table["residual"]) - set(NAME_TO_ID))
    if unknown:
        raise ValueError(f"unknown beta-plus residuals: {', '.join(unknown)}")
    for field in ("projectile", "target", "residual"):
        if table[field].astype(str).str.strip().eq("").any():
            raise ValueError(f"native-route field '{field}' must be non-empty")
    return table.reset_index(drop=True)


def classify_native_routes(
    frame: pd.DataFrame,
    *,
    scenario_names=("d120s300",),
    scenario_config: Path | str = DEFAULT_SCENARIO_CONFIG,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return route-level classifications and modeled/unmodeled fractions."""

    table = validate_native_routes(frame)
    represented = np.asarray([
        projectile.lower() == "proton" and (target, residual) in SELECTED_ROUTES
        for projectile, target, residual in table[
            ["projectile", "target", "residual"]
        ].itertuples(index=False, name=None)
    ])
    table["represented_by_fold"] = represented
    table["route_class"] = np.where(represented, "modeled", "unmodeled")

    weight_columns = {"all_production": np.ones(len(table), dtype=float)}
    for name in scenario_names:
        scenario = resolve_scenario(name, scenario_config)
        weight_columns[f"all_{name}"] = np.asarray([
            scenario.measured_fraction(ISOTOPES[NAME_TO_ID[residual]].lam)
            for residual in table["residual"]
        ])

    summary_rows = []
    for profile_label, factors in weight_columns.items():
        table[f"weight_{profile_label}"] = factors
        table[f"weighted_count_{profile_label}"] = table["production_count"] * factors
        values = table[f"weighted_count_{profile_label}"].to_numpy(dtype=float)
        modeled = float(np.sum(values[represented]))
        unmodeled = float(np.sum(values[~represented]))
        total = modeled + unmodeled
        summary_rows.append({
            "profile_label": profile_label,
            "modeled_count": modeled,
            "unmodeled_count": unmodeled,
            "total_count": total,
            "unmodeled_fraction": unmodeled / total if total > 0 else float("nan"),
        })
    return table, pd.DataFrame(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_route_csv", type=Path)
    parser.add_argument(
        "--scenario", action="append", default=None,
        help="named beam-off scenario; may be supplied more than once")
    parser.add_argument("--scenario-config", type=Path, default=DEFAULT_SCENARIO_CONFIG)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()
    scenarios = arguments.scenario or ["d120s300"]
    detailed, summary = classify_native_routes(
        pd.read_csv(arguments.native_route_csv),
        scenario_names=scenarios,
        scenario_config=arguments.scenario_config,
    )
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    detailed.to_csv(arguments.output_dir / "native_route_classification.csv", index=False)
    summary.to_csv(arguments.output_dir / "native_route_summary.csv", index=False)
    source = arguments.native_route_csv.resolve()
    scenario_records = []
    for name in scenarios:
        scenario = resolve_scenario(name, arguments.scenario_config)
        scenario_records.append({
            "name": scenario.name,
            "t_irr_s": scenario.t_irr_s,
            "t_del_s": scenario.t_del_s,
            "t_meas_s": scenario.t_meas_s,
            "config_sha256": scenario.config_sha256,
        })
    metadata = {
        "schema_version": 1,
        "source_file": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "scenarios": scenario_records,
        "selected_channel_count": len(SELECTED_ROUTES),
        "affects_folded_source": False,
    }
    (arguments.output_dir / "native_route_meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
