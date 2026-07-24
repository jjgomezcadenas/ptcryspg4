"""Authoritative named acquisition scenarios shared by folding and handoff."""

from __future__ import annotations

import hashlib
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path


DEFAULT_SCENARIO_CONFIG = (
    Path(__file__).resolve().parents[1] / "config/handoff_scenarios.toml")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def survival(lam: float, t_irr: float, t_del: float, t_meas: float):
    """Return build-up, delay-survival, and acquisition-window factors."""

    if lam <= 0 or t_irr <= 0 or t_del < 0 or t_meas <= 0:
        raise ValueError("decay constant and acquisition times are outside their domains")
    build = -math.expm1(-lam * t_irr) / (lam * t_irr)
    transport = math.exp(-lam * t_del)
    window = -math.expm1(-lam * t_meas)
    return build, transport, window


@dataclass(frozen=True)
class HandoffScenario:
    name: str
    description: str
    t_irr_s: float
    t_del_s: float
    t_meas_s: float
    config_path: Path
    config_sha256: str

    def factors(self, lam: float) -> tuple[float, float, float]:
        return survival(lam, self.t_irr_s, self.t_del_s, self.t_meas_s)

    def measured_fraction(self, lam: float) -> float:
        build, transport, window = self.factors(lam)
        return build * transport * window


def resolve_scenario(
    name: str,
    config_path: Path | str = DEFAULT_SCENARIO_CONFIG,
) -> HandoffScenario:
    path = Path(config_path).resolve()
    document = tomllib.loads(path.read_text(encoding="utf-8"))
    if int(document.get("schema_version", 0)) != 1:
        raise ValueError("unsupported handoff-scenario schema version")
    scenarios = document.get("scenario", {})
    if name not in scenarios:
        available = ", ".join(sorted(scenarios))
        raise ValueError(f"unknown handoff scenario '{name}'; available: {available}")
    row = scenarios[name]
    times = {field: float(row[field]) for field in ("t_irr_s", "t_del_s", "t_meas_s")}
    if times["t_irr_s"] <= 0 or times["t_del_s"] < 0 or times["t_meas_s"] <= 0:
        raise ValueError(f"handoff scenario '{name}' contains invalid times")
    return HandoffScenario(
        name=name,
        description=str(row.get("description", "")),
        config_path=path,
        config_sha256=file_sha256(path),
        **times,
    )
