#!/usr/bin/env python3
"""Regenerate and compile the standalone EXFOR cross-section fit report."""

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
    parser.add_argument("--skip-normalize", action="store_true")
    parser.add_argument("--skip-fit", action="store_true")
    parser.add_argument("--skip-latex", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    python = Path(sys.executable)
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", str(repo / ".cache/matplotlib"))
    environment.setdefault("PYTHONPYCACHEPREFIX", "/tmp/ptcrysp-pycache")
    if not args.skip_normalize:
        run([python, "-m", "analysis_transport.xsections.normalize"],
            repo, environment)
    run([python, "-m", "analysis_transport.xsections.reaction_thresholds"],
        repo, environment)
    run([python, "-m", "analysis_transport.xsections.curate_exfor"],
        repo, environment)
    if not args.skip_fit:
        run([python, "-m", "analysis_transport.xsections.fit_exfor"],
            repo, environment)
    run([python, "-m", "analysis_transport.xsections.make_fit_plots"],
        repo, environment)
    run([python, "-m", "analysis_transport.xsections.validate_fit"],
        repo, environment)
    if not args.skip_latex:
        run([python, "docs/build_latex.py", "xsection_fit"], repo)


if __name__ == "__main__":
    main()
