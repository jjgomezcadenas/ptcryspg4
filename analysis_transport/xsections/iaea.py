"""Reader for the IAEA recommended O-16(p,alpha)N-13 table."""

from pathlib import Path

from .common import point


def read(path: Path):
    values = []
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        fields = raw_line.split()
        if len(fields) < 3:
            continue
        try:
            values.append(tuple(float(value) for value in fields[:3]))
        except ValueError:
            continue
    if not values:
        raise ValueError(f"No numerical IAEA points in {path}")
    return {
        "dataset_id": "iaea2021_p_O16_a_N13",
        "library": "IAEA recommended",
        "target": "O16",
        "residual": "N13",
        "label": "IAEA recommended O-16(p,alpha)N-13",
        "original_energy_unit": "MeV",
        "original_cross_section_unit": "mb",
        "role": "recommended_evaluation",
        "source_file": str(path),
        "points": [point(i, energy, sigma, sigma_unc=uncertainty)
                   for i, (energy, sigma, uncertainty) in enumerate(values)],
    }
