"""Focused ENDF-6 reader for residual production in MF=3/6, MT=5."""

import math
from pathlib import Path

import numpy as np

from .common import point


def endf_float(field: str) -> float:
    text = field.strip()
    if not text:
        return 0.0
    if "e" not in text.lower():
        for index in range(len(text) - 1, 0, -1):
            if text[index] in "+-" and text[index - 1].isdigit():
                text = text[:index] + "e" + text[index:]
                break
    return float(text)


def _record(line):
    fields = [line[index:index + 11] for index in range(0, 66, 11)]
    return (
        endf_float(fields[0]), endf_float(fields[1]),
        int(endf_float(fields[2])), int(endf_float(fields[3])),
        int(endf_float(fields[4])), int(endf_float(fields[5])),
    )


def _section(lines, mf, mt):
    return [line for line in lines if int(line[70:72] or 0) == mf and int(line[72:75] or 0) == mt]


def _tab1(section, header_index):
    _, _, _, _, nr, npairs = _record(section[header_index])
    n_interpolation_lines = math.ceil(2 * nr / 6)
    first_data = header_index + 1 + n_interpolation_lines
    numbers = []
    for line in section[first_data:first_data + math.ceil(2 * npairs / 6)]:
        numbers.extend(endf_float(line[index:index + 11]) for index in range(0, 66, 11))
    return np.asarray(numbers[0:2 * npairs:2]), np.asarray(numbers[1:2 * npairs:2])


def read(path: Path, target: str, residual: str, residual_z: int, residual_a: int):
    lines = path.read_text(encoding="ascii").splitlines()
    mf3 = _section(lines, 3, 5)
    if len(mf3) < 2:
        raise ValueError(f"MF=3 MT=5 is absent from {path}")
    total_energy, total_sigma = _tab1(mf3, 1)

    zap = residual_z * 1000 + residual_a
    mf6 = _section(lines, 6, 5)
    product_index = None
    for index, line in enumerate(mf6[1:], 1):
        c1, _, _, law, nr, npairs = _record(line)
        if round(c1) == zap and law >= 0 and nr > 0 and npairs > 0:
            product_index = index
            break
    if product_index is None:
        raise ValueError(f"Residual ZAP={zap} is absent from MF=6 MT=5 in {path}")
    energy, yield_value = _tab1(mf6, product_index)
    sigma = np.interp(energy, total_energy, total_sigma) * yield_value
    dataset_id = f"lanl_endfb71_p_{target}_x_{residual}"
    return {
        "dataset_id": dataset_id,
        "library": "LANL ENDF/B-VII.1",
        "target": target,
        "residual": residual,
        "label": "LANL ENDF/B-VII.1",
        "original_energy_unit": "eV",
        "original_cross_section_unit": "b",
        "transformation": "MF=3 MT=5 cross section times MF=6 MT=5 residual yield",
        "source_file": str(path),
        "points": [point(i, e / 1.0e6, s * 1.0e3)
                   for i, (e, s) in enumerate(zip(energy, sigma))],
    }
