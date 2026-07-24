"""Reader for the pinned EXFOR Plot File CX4 representation."""

import re
from pathlib import Path

from .common import point


REACTION_RE = re.compile(r"(\d+)-([A-Z][A-Z]?)-(\d+)\(P,[^)]+\)(\d+)-([A-Z][A-Z]?)-(\d+)", re.I)


def read(path: Path):
    metadata = {}
    values = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith("#") and ":" in line:
            key, value = line[1:].split(":", 1)
            metadata[key.strip()] = value.strip()
        elif line and not line.startswith("#"):
            fields = line.split()
            if len(fields) >= 4:
                values.append(tuple(float(value) for value in fields[:4]))
    reaction = metadata.get("REACTION", "")
    match = REACTION_RE.search(reaction)
    if not match:
        raise ValueError(f"Cannot identify reaction in {path}: {reaction}")
    target = f"{match.group(2).upper().title()}{int(match.group(3))}"
    residual = f"{match.group(5).upper().title()}{int(match.group(6))}"
    accession = metadata.get("EXFOR #", path.stem).strip()
    dataset_id = "exfor_" + accession.replace(".", "_") + f"_p_{target}_x_{residual}"
    points = [
        point(index, energy / 1.0e6, sigma * 1.0e3,
              energy_unc / 1.0e6, sigma_unc * 1.0e3)
        for index, (energy, energy_unc, sigma, sigma_unc) in enumerate(values)
    ]
    normalized = {
        "dataset_id": dataset_id,
        "library": "EXFOR",
        "target": target,
        "residual": residual,
        "label": f"{metadata.get('Author', accession)} ({accession})",
        "accession": accession,
        "reaction": reaction,
        "author": metadata.get("Author", ""),
        "reference": metadata.get("Reference", ""),
        "doi": metadata.get("DOI", ""),
        "quantity": metadata.get("Quantity", ""),
        "original_energy_unit": "eV",
        "original_cross_section_unit": "b",
        "role": "experimental",
        "source_file": str(path),
        "points": points,
    }
    if accession.upper().startswith("E2568"):
        normalized["role"] = "shape_validation"
        normalized["normalization_note"] = (
            "Masuda oxygen curves share an external normalization through "
            "O-16(p,x)O-15 = 76.8 mb at 35 MeV")
    return normalized
