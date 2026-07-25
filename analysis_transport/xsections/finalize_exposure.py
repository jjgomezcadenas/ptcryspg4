"""Complete and validate the metadata of an exposure run.

The Geant4 exposure application writes ``run_meta_raw.json`` next to its
``proton_exposure.csv``.  This tool adds the exposure-file digest, the
software revision and the per-Gy normalization, writes the validated
``exposure_meta.json``, and runs the exposure-table, metadata and
native-route validators on the finished run directory.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from .exposure_folding import (
    load_exposure_metadata,
    load_exposure_table,
    sha256,
)
from .native_routes import validate_native_routes


def software_revision(repo: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True)
    return result.stdout.strip()


def finalize(run_dir: Path, repo: Path) -> Path:
    raw_path = run_dir / "run_meta_raw.json"
    document = json.loads(raw_path.read_text(encoding="utf-8"))
    exposure_path = run_dir / document["exposure_file"]

    n_protons = float(document["n_protons"])
    target_dose = float(document["target_dose_Gy"])
    if target_dose <= 0:
        raise ValueError(
            "target_dose_Gy is not positive; the run's target box does not "
            "intersect the phantom or the run deposited no target dose")
    document["exposure_sha256"] = sha256(exposure_path)
    document["software_revision"] = software_revision(repo)
    document["Np_per_Gy"] = n_protons / target_dose

    meta_path = run_dir / "exposure_meta.json"
    meta_path.write_text(json.dumps(document, indent=2) + "\n",
                         encoding="utf-8")

    metadata = load_exposure_metadata(meta_path, exposure_path)
    table = load_exposure_table(exposure_path)
    native = validate_native_routes(
        pd.read_csv(run_dir / "native_route_counts.csv"))
    print(f"finalized {meta_path}")
    print(f"  run {metadata.run_id}: {metadata.n_protons:.0f} protons, "
          f"target dose {metadata.target_dose_Gy:.4e} Gy, "
          f"Np_per_Gy {metadata.Np_per_Gy:.6e}")
    print(f"  exposure rows {len(table)}, native-route rows {len(native)}")
    return meta_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    finalize(args.run_dir.resolve(), args.repo.resolve())


if __name__ == "__main__":
    main()
