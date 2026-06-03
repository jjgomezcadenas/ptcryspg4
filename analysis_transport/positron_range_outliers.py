#!/usr/bin/env python3
"""Diagnose the positron-range tail in Stage-A output.

The mean |anh-prod| runs well above the median because a small fraction of
positrons escape the phantom into the surrounding air (where their range is
huge) and annihilate hundreds of mm away. This script quantifies that tail:
how many annihilations fall OUTSIDE the phantom, where those emitters were
produced, and what the per-isotope range statistics look like once the escaped
positrons are removed (which should leave clean, endpoint-ordered values).

Phantom dimensions are read from run_meta.csv, so this stays correct if the
geometry changes.

Usage:
    python analysis_transport/positron_range_outliers.py [data_dir]
"""

import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
from isotopes import ISOTOPES  # noqa: E402

ISO_BY_ENDPOINT = sorted(ISOTOPES, key=lambda k: -ISOTOPES[k].endpoint_MeV)


def main(data_dir: str) -> None:
    meta = pd.read_csv(os.path.join(data_dir, "run_meta.csv")).iloc[0]
    d = pd.read_csv(os.path.join(data_dir, "emitters.csv"))

    radius = 0.5 * float(meta["phantom_diameter_mm"])
    half_z = 0.5 * float(meta["phantom_length_mm"])

    prod = d[["prod_x_mm", "prod_y_mm", "prod_z_mm"]].to_numpy()
    anh = d[["anh_x_mm", "anh_y_mm", "anh_z_mm"]].to_numpy()
    d["range_mm"] = np.linalg.norm(anh - prod, axis=1)
    d["prod_r_mm"] = np.hypot(d.prod_x_mm, d.prod_y_mm)
    d["anh_r_mm"] = np.hypot(d.anh_x_mm, d.anh_y_mm)

    # An annihilation is "escaped" if it lands outside the phantom volume.
    inside = (d.anh_r_mm < radius) & (d.prod_z_mm.between(-half_z, half_z)) \
        & (d.anh_z_mm.between(-half_z, half_z))
    escaped = ~inside

    n = len(d)
    print(f"\nphantom: radius {radius:.0f} mm, half-length {half_z:.0f} mm")
    print(f"emitters total      : {n}")
    print(f"escaped (anh outside): {escaped.sum()}  "
          f"({100 * escaped.mean():.2f}%)")
    print(f"range > 10 mm        : {(d.range_mm > 10).sum()}")
    print(f"  of those, escaped  : {(escaped & (d.range_mm > 10)).sum()}"
          f"   (escape explains the tail if these match)")

    esc = d[escaped]
    if len(esc):
        print(f"\nescaped emitters were produced at:")
        print(f"  prod radius: mean {esc.prod_r_mm.mean():.1f} mm "
              f"(bulk {d[inside].prod_r_mm.mean():.1f} mm)")
        print(f"  prod z     : min {esc.prod_z_mm.min():.1f}, "
              f"max {esc.prod_z_mm.max():.1f} mm  "
              f"(near a face means the e+ exits into air)")

    print(f"\nper-isotope positron range, escaped positrons removed:")
    print(f"{'iso':>5} {'endpt':>6} {'count':>7} {'mean':>8} {'median':>8} "
          f"{'mean(all)':>10}")
    print("-" * 52)
    kept = d[inside]
    for iid in ISO_BY_ENDPOINT:
        iso = ISOTOPES[iid]
        s = kept[kept.isotope_id == iid].range_mm
        s_all = d[d.isotope_id == iid].range_mm
        if len(s) == 0:
            continue
        print(f"{iso.name:>5} {iso.endpoint_MeV:>6.2f} {len(s):>7} "
              f"{s.mean():>6.3f}mm {s.median():>6.3f}mm {s_all.mean():>8.3f}mm")
    print("\n(mean now tracks the endpoint ordering once escapes are removed;\n"
          " median was already robust to the tail.)")


if __name__ == "__main__":
    ddir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "..", "data")
    main(ddir)
