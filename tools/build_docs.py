#!/usr/bin/env python3
"""Build the HTML code reference from the docstrings and comments.

Two generators, one output tree (docs_api/, gitignored):
  docs_api/python/  pdoc renders every Python module's docstrings as-is
                    (plain prose, the house style) into one linked site.
  docs_api/stageA/  Doxygen renders the Geant4 C++ app from its /// comments
                    (config: tools/Doxyfile).

test_field.py is a test gate that runs its checks at import, so it is not part
of the reference. Requires pdoc (pip install pdoc) and doxygen (brew install
doxygen); a missing generator is reported and the other still runs.

Usage:
    python3 tools/build_docs.py            # build both
    python3 tools/build_docs.py python     # only the Python reference
    python3 tools/build_docs.py stageA     # only the C++ reference
"""

import argparse
import os
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)
OUT = os.path.join(ROOT, "docs_api")

# Every documented module, by pipeline stage. pdoc imports each file as a
# top-level module, so the basenames must stay distinct.
PY_MODULES = [
    "common/isotopes.py",
    "common/regions.py",
    "common/phantom_material.py",
    "field_design/sobp.py",
    "field_design/plot_sobp.py",
    "decay_sampling/budget.py",
    "decay_sampling/budget_gen.py",
    "decay_sampling/activity_plot.py",
    "analysis_transport/check_run.py",
    "analysis_transport/validate_transport.py",
    "analysis_transport/sobp_metrics.py",
    "analysis_transport/distal_pool.py",
    "analysis_transport/parodi_cross_check.py",
    "analysis_transport/make_figures.py",
    "analysis_transport/compare_beams.py",
    "analysis_transport/bragg_profile.py",
    "analysis_transport/plot_depth_dose.py",
    "analysis_transport/plot_dose_activity.py",
    "analysis_transport/plot_geometry.py",
    "analysis_transport/plot_phantom.py",
    "analysis_transport/plot_mird_head.py",
    "tools/snapshot_scenario.py",
    "latex/build_latex.py",
]


def build_python():
    """pdoc over PY_MODULES -> docs_api/python (index + one page per module)."""
    out = os.path.join(OUT, "python")
    cmd = [sys.executable, "-m", "pdoc", "-o", out] + PY_MODULES
    subprocess.run(cmd, cwd=ROOT, check=True)
    print(f"  python reference -> {out}/index.html")


def build_stageA():
    """doxygen with tools/Doxyfile -> docs_api/stageA."""
    if shutil.which("doxygen") is None:
        print("  doxygen not found (brew install doxygen) — skipped stageA")
        return
    os.makedirs(os.path.join(OUT, "stageA"), exist_ok=True)
    subprocess.run(["doxygen", os.path.join(_HERE, "Doxyfile")],
                   cwd=ROOT, check=True)
    print(f"  stageA reference -> {os.path.join(OUT, 'stageA')}/index.html")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", default="all",
                    choices=["all", "python", "stageA"])
    args = ap.parse_args()

    print(f"building code reference in {OUT}")
    if args.target in ("all", "python"):
        build_python()
    if args.target in ("all", "stageA"):
        build_stageA()


if __name__ == "__main__":
    main()
