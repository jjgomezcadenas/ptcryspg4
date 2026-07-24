"""Run direct Geant4 inelastic final-state sampling for the denominator."""

import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from .config import load


COUNT_FIELDS = ("n_interactions", "n_c11", "n_n13", "n_o15",
                "n_secondaries", "n_nuclei")


def invoke(executable: Path, target, energy, interactions, seed,
           physics_list, output: Path, log: Path):
    command = [
        str(executable), "--target", target, "--energy-mev", str(energy),
        "--interactions", str(interactions), "--seed", str(seed),
        "--physics-list", physics_list, "--output", str(output),
    ]
    with log.open("w", encoding="utf-8") as stream:
        subprocess.run(command, check=True, stdout=stream,
                       stderr=subprocess.STDOUT)
    with output.open(newline="", encoding="utf-8") as stream:
        return next(csv.DictReader(stream))


def merge_rows(rows):
    if not rows:
        raise ValueError("Cannot merge an empty denominator batch list")
    merged = dict(rows[0])
    for field in COUNT_FIELDS:
        merged[field] = sum(int(row[field]) for row in rows)
    for field in ("sigma_inelastic_mb", "target_z", "target_a"):
        if len({row[field] for row in rows}) != 1:
            raise ValueError(f"Inconsistent denominator batch field: {field}")
    for field in ("cross_section_data_sets", "physics_list", "geant4_version"):
        values = {row[field] for row in rows}
        if len(values) != 1:
            raise ValueError(f"Inconsistent denominator provenance: {field}")
    merged["batch_count"] = len(rows)
    merged["seeds"] = ";".join(row["seed"] for row in rows)
    return merged


def _write_rows(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def active_residuals(config, target, energy):
    target_config = config["target"][target]
    return [residual for residual in target_config["active_residuals"]
            if energy >= target_config["threshold_MeV"][residual]]


def run(config, executable: Path, destination: Path, mode, targets, energies):
    pilot = mode == "pilot"
    output_name = ("denominator_pilot_counts.csv" if pilot
                   else "denominator_counts.csv")
    logs = destination / "logs" / f"denominator_{mode}"
    logs.mkdir(parents=True, exist_ok=True)
    rows = []
    with tempfile.TemporaryDirectory(
            prefix="ptcrysp-g4-denominator-") as temporary:
        for target in targets:
            for energy in energies:
                batches = []
                total = 0
                batch_index = 0
                residuals = active_residuals(config, target, energy)
                while True:
                    if pilot:
                        interactions = int(config["denominator_pilot_interactions"])
                    else:
                        interactions = min(
                            int(config["denominator_batch_interactions"]),
                            int(config["denominator_maximum_interactions"]) - total)
                    output = Path(temporary) / "sample.csv"
                    if output.exists():
                        output.unlink()
                    seed = (3_000_000 + sum(ord(character) for character in target) * 10_000
                            + int(energy * 10) * 100 + batch_index)
                    row = invoke(
                        executable, target, energy, interactions, seed,
                        config["physics_list"], output,
                        logs / f"{target}_{energy:g}_{batch_index:04d}.log")
                    batches.append(row)
                    total += interactions
                    merged = merge_rows(batches)
                    if pilot or not residuals:
                        break
                    counts = [int(merged[f"n_{residual.lower()}"])
                              for residual in residuals]
                    if min(counts) >= int(config["denominator_minimum_residual_count"]):
                        break
                    if total >= int(config["denominator_maximum_interactions"]):
                        break
                    batch_index += 1
                if not residuals:
                    merged["stopping_condition"] = "below_threshold"
                elif pilot:
                    merged["stopping_condition"] = "pilot_interactions"
                elif min(int(merged[f"n_{residual.lower()}"])
                         for residual in residuals) >= int(
                             config["denominator_minimum_residual_count"]):
                    merged["stopping_condition"] = "residual_count"
                else:
                    merged["stopping_condition"] = "interaction_cap"
                rows.append(merged)
    _write_rows(destination / output_name, rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("pilot", "production"))
    parser.add_argument("--config", type=Path,
                        default=Path("config/xsections_denominator.toml"))
    parser.add_argument("--executable", type=Path,
                        default=Path("xsections_g4/build/denominator_sampler"))
    parser.add_argument("--output-root", type=Path,
                        default=Path("data/xsections/g4"))
    parser.add_argument("--targets", nargs="*")
    parser.add_argument("--energies", nargs="*", type=float)
    args = parser.parse_args()
    config_path = args.config.resolve()
    executable = args.executable.resolve()
    config = load(config_path)
    targets = args.targets or config["targets"]
    energies = args.energies or config["energies_MeV"]
    destination = (args.output_root.resolve()
                   / f"geant4-{config['geant4_version']}"
                   / config["physics_list"])
    destination.mkdir(parents=True, exist_ok=True)
    rows = run(config, executable, destination, args.mode, targets, energies)
    metadata = {
        "schema_version": 1,
        "mode": args.mode,
        "config": str(config_path),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "executable": str(executable),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "targets": targets,
        "energies_MeV": energies,
        "row_count": len(rows),
    }
    (destination / f"denominator_{args.mode}_meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    print(f"Wrote {len(rows)} direct Geant4 denominator {args.mode} rows "
          f"under {destination}")


if __name__ == "__main__":
    main()
