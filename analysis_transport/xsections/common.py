"""Common normalized-table helpers."""

import csv
from pathlib import Path


POINT_FIELDS = (
    "point_id",
    "energy_MeV",
    "energy_unc_minus_MeV",
    "energy_unc_plus_MeV",
    "sigma_mb",
    "sigma_unc_stat_mb",
    "sigma_unc_sys_mb",
    "sigma_unc_minus_mb",
    "sigma_unc_plus_mb",
)


def point(point_id, energy, sigma, energy_unc="", sigma_unc=""):
    return {
        "point_id": point_id,
        "energy_MeV": energy,
        "energy_unc_minus_MeV": energy_unc,
        "energy_unc_plus_MeV": energy_unc,
        "sigma_mb": sigma,
        "sigma_unc_stat_mb": "",
        "sigma_unc_sys_mb": "",
        "sigma_unc_minus_mb": sigma_unc,
        "sigma_unc_plus_mb": sigma_unc,
    }


def write_points(path: Path, points):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=POINT_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(points)
