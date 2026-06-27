"""Shared phantom-medium geometry: point-in-region test and material lookup.

The medium is `phantom_regions.csv`: priority-ordered, world-frame, axis-aligned
solids (ellipsoid or cylinder). A point takes the material of the first
(lowest-priority) region that contains it, else `None`/`"air"`. This is the one
implementation of that rule, shared by `check_run.py`, `plot_mird_head.py`, and
the WEPL ray-trace in `field_design/sobp.py`.

A "region" is any object with attributes `solid, a_mm, b_mm, c_mm, cx_mm, cy_mm,
cz_mm` (e.g. a pandas `itertuples` row); for ellipsoids (a,b,c) are the semi-axes,
for a cylinder (a,b,c) = (radius, radius, half-length).
"""

import numpy as np


def contains(region, x, y, z):
    """True if world point (x, y, z) [mm] is inside `region`."""
    dx, dy, dz = x - region.cx_mm, y - region.cy_mm, z - region.cz_mm
    if region.solid == "cylinder":
        return dx * dx + dy * dy <= region.a_mm ** 2 and abs(dz) <= region.c_mm
    return (dx / region.a_mm) ** 2 + (dy / region.b_mm) ** 2 + (dz / region.c_mm) ** 2 <= 1.0


def material_at(regions, x, y, z):
    """NIST material of the first priority-ordered region containing the point,
    else None (air). `regions` is the `phantom_regions.csv` DataFrame."""
    for r in regions.sort_values("priority").itertuples():
        if contains(r, x, y, z):
            return r.material
    return None


def classify(regions, X, Y, Z):
    """Region name of each point in the arrays (X, Y, Z) [mm], 'air' if in none.
    Priority-ordered: the lowest-priority region containing a point wins."""
    X, Y, Z = np.asarray(X, float), np.asarray(Y, float), np.asarray(Z, float)
    out = np.full(X.shape, "air", dtype=object)
    assigned = np.zeros(X.shape, dtype=bool)
    for _, r in regions.sort_values("priority").iterrows():
        dx, dy, dz = X - r.cx_mm, Y - r.cy_mm, Z - r.cz_mm
        if r.solid == "cylinder":
            inside = (dx * dx + dy * dy <= r.a_mm ** 2) & (np.abs(dz) <= r.c_mm)
        else:
            inside = (dx / r.a_mm) ** 2 + (dy / r.b_mm) ** 2 + (dz / r.c_mm) ** 2 <= 1.0
        take = inside & ~assigned
        out[take] = r.region
        assigned |= inside
    return out
