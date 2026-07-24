"""Reader for JENDL residual-production text tables."""

import re
from pathlib import Path

from .common import point


BLOCK_RE = re.compile(r"^# Production of ([A-Za-z]+)-(\d+)\s*$")


def read(path: Path, residual: str):
    wanted = re.sub(r"(\D+)(\d+)", r"\1-\2", residual)
    values = []
    active = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        match = BLOCK_RE.match(raw_line.strip())
        if match:
            active = f"{match.group(1)}-{match.group(2)}" == wanted
            continue
        if active and raw_line.startswith("# Production"):
            break
        if active and raw_line.strip() and not raw_line.startswith("#"):
            fields = raw_line.split()
            if len(fields) >= 2:
                values.append((float(fields[0]), float(fields[1])))
    if not values:
        raise ValueError(f"Residual {residual} is absent from {path}")
    target = path.stem.replace("0", "", 1) if path.stem[1:2] == "0" else path.stem
    dataset_id = f"jendl40he_p_{target}_x_{residual}"
    return {
        "dataset_id": dataset_id,
        "library": "JENDL-4.0/HE",
        "target": target,
        "residual": residual,
        "label": "JENDL-4.0/HE",
        "original_energy_unit": "eV",
        "original_cross_section_unit": "b",
        "transformation": "JAEA precomputed MF=3 MT=5 times MF=6 MT=5 residual production",
        "source_file": str(path),
        "points": [point(i, energy / 1.0e6, sigma * 1.0e3)
                   for i, (energy, sigma) in enumerate(values)],
    }
