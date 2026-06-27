#!/usr/bin/env python3
"""Geometry-aware control-figure dispatcher for one Stage-A run directory.

Reads <run_dir>/run_meta.csv, picks the plotters that make sense for the run's
geometry, and writes their PNGs into <run_dir>/figures/ (each plotter creates
that subdir). The snapshot tool then copies <run_dir>/figures/ verbatim, so a
frozen scenario always ships the figures appropriate to what it actually is.

Dispatch:
  any geometry          transport_validation.png  (Stage-A dashboard)
                        activity.png              (decay/activity curves)
  cylinder              + sobp_g4.png             (realized depth-dose plateau)
  uniform_head/mird_head + mird_head.png          (phantom + beam + emitter trail)

Plotters that need a budget (none here) or other optional inputs are skipped
silently if their inputs are absent.

Usage:
    python analysis_transport/make_figures.py <run_dir>
"""

import argparse
import os
import subprocess
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")


def run_plotter(script, run_dir):
    """Invoke a plotter as a subprocess pointed at the run dir."""
    path = os.path.join(_ROOT, script)
    print(f"  -> {script}")
    subprocess.run([sys.executable, path, run_dir], check=True)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a run directory, e.g. data/runs/cylinder_sobp_1e7")
    args = ap.parse_args()

    meta_path = os.path.join(args.run_dir, "run_meta.csv")
    if not os.path.exists(meta_path):
        sys.exit(f"no run_meta.csv in {args.run_dir} — not a Stage-A run directory")
    geometry = str(pd.read_csv(meta_path).iloc[0]["geometry"])

    print(f"figures for {args.run_dir}  (geometry: {geometry})")
    # Geometry-independent: phantom+target drawing, central-axis depth dose,
    # dashboard, activity curves.
    run_plotter("analysis_transport/plot_phantom.py", args.run_dir)
    run_plotter("analysis_transport/plot_depth_dose.py", args.run_dir)
    run_plotter("analysis_transport/validate_transport.py", args.run_dir)
    run_plotter("decay_sampling/activity_plot.py", args.run_dir)
    # Geometry-specific control plot.
    if geometry == "cylinder":
        run_plotter("field_design/plot_sobp.py", args.run_dir)
    else:  # uniform_head | mird_head
        run_plotter("analysis_transport/plot_mird_head.py", args.run_dir)

    print(f"figures -> {os.path.join(args.run_dir, 'figures')}")


if __name__ == "__main__":
    main()
