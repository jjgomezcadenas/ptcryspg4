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
    parser.add_argument(
        "--folding-dir", type=Path,
        help="complete folding output directory used for propagated figures")
    parser.add_argument(
        "--convergence-csv", type=Path,
        help="energy-grid convergence table used for the fourth propagated figure")
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
    run([python, "-m", "analysis_transport.xsections.make_folding_reference"],
        repo, environment)
    if args.convergence_csv is not None and args.folding_dir is None:
        parser.error("--convergence-csv requires --folding-dir")
    if args.folding_dir is not None:
        plot_command = [
            python, "-m", "analysis_transport.xsections.make_folding_plots",
            args.folding_dir.resolve(),
            "--output-dir", repo / "docs/figures/xsection_folding",
            "--generated-dir", repo / "docs/generated/xsection_folding",
        ]
        if args.convergence_csv is not None:
            plot_command.extend([
                "--convergence-csv", args.convergence_csv.resolve()])
        run(plot_command, repo, environment)
    if not args.skip_latex:
        run([python, "docs/build_latex.py", "xsection_fit"], repo)


if __name__ == "__main__":
    main()
