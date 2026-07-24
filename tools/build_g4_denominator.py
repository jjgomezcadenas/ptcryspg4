#!/usr/bin/env python3
"""Build, run, plot, validate, and document the Geant4 denominator pilot."""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def run(command, repo, environment=None):
    print(" ".join(str(item) for item in command))
    subprocess.run([str(item) for item in command], cwd=repo, check=True,
                   env=environment)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--skip-pilot", action="store_true")
    parser.add_argument("--skip-latex", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    python = Path(sys.executable)
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", str(repo / ".cache/matplotlib"))
    environment.setdefault("PYTHONPYCACHEPREFIX", "/tmp/ptcrysp-pycache")
    run(["cmake", "-S", "xsections_g4", "-B", "xsections_g4/build"], repo)
    run(["cmake", "--build", "xsections_g4/build", "-j4"], repo)
    if not args.skip_pilot:
        run([python, "-m", "analysis_transport.xsections.run_g4_denominator",
             "pilot"], repo, environment)
    run([python, "-m", "analysis_transport.xsections.build_g4_denominator"],
        repo, environment)
    run([python, "-m", "analysis_transport.xsections.make_g4_denominator_plots"],
        repo, environment)
    run([python, "-m", "analysis_transport.xsections.validate_g4_denominator"],
        repo, environment)
    if not args.skip_latex:
        run([python, "docs/build_latex.py", "xsections_plan"], repo, environment)


if __name__ == "__main__":
    main()
