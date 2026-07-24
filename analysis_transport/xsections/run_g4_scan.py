"""Run the configured Geant4 thin-target pilot and production scan."""

import argparse
import csv
import json
import math
import subprocess
import tempfile
from pathlib import Path

from .config import load


COUNT_FIELDS = ("n_protons", "n_inelastic", "n_c11", "n_n13", "n_o15")
ENERGY_SUM_FIELDS = (
    ("mean_continuous_loss_MeV", "n_protons"),
    ("mean_inelastic_energy_MeV", "n_inelastic"),
    ("mean_c11_energy_MeV", "n_c11"),
    ("mean_n13_energy_MeV", "n_n13"),
    ("mean_o15_energy_MeV", "n_o15"),
)


def invoke(executable: Path, target, energy, thickness, protons, seed, threads,
           physics_list, output: Path, log: Path):
    command = [
        str(executable), "--target", target, "--energy-mev", str(energy),
        "--areal-mg-cm2", str(thickness), "--protons", str(protons),
        "--seed", str(seed), "--threads", str(threads),
        "--physics-list", physics_list, "--output", str(output),
    ]
    with log.open("w", encoding="utf-8") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
    with output.open(newline="", encoding="utf-8") as stream:
        return next(csv.DictReader(stream))


def select_thickness(rows, config):
    max_loss = config["maximum_mean_energy_loss_MeV"]
    max_probability = config["maximum_inelastic_probability"]
    for row in rows:
        n = int(row["n_protons"])
        # Three events is the 95% upper mean for an observed zero.
        conservative_probability = (int(row["n_inelastic"]) + 3.0) / n
        if (float(row["mean_continuous_loss_MeV"]) < max_loss
                and conservative_probability < max_probability):
            return float(row["areal_mg_cm2"])
    raise RuntimeError("No configured thickness passes the pilot criteria")


def merge_rows(rows):
    if not rows:
        raise ValueError("Cannot merge an empty run list")
    merged = dict(rows[0])
    for output_field, count_field in ENERGY_SUM_FIELDS:
        denominator = sum(int(row[count_field]) for row in rows)
        numerator = sum(float(row[output_field]) * int(row[count_field]) for row in rows)
        merged[output_field] = numerator / denominator if denominator else 0.0
    for field in COUNT_FIELDS:
        merged[field] = sum(int(row[field]) for row in rows)
    merged["batch_count"] = len(rows)
    return merged


def _write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_pilot(config, executable, destination, threads, targets, energies):
    rows = []
    logs = destination / "logs/pilot"
    logs.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ptcrysp-xsections-") as temporary:
        for target in targets:
            for energy in energies:
                candidate_rows = []
                for candidate_index, thickness in enumerate(config["thickness_candidates_mg_cm2"]):
                    output = Path(temporary) / "counts.csv"
                    if output.exists(): output.unlink()
                    seed = 1000000 + sum(ord(c) for c in target) * 1000 + int(energy * 10) * 10 + candidate_index
                    row = invoke(executable, target, energy, thickness,
                                 config["pilot_protons"], seed, threads,
                                 config["physics_list"], output,
                                 logs / f"{target}_{energy:g}_{thickness:g}.log")
                    candidate_rows.append(row)
                    try:
                        selected = select_thickness(candidate_rows, config)
                        row["selected"] = 1
                        rows.extend(candidate_rows)
                        break
                    except RuntimeError:
                        row["selected"] = 0
                else:
                    raise RuntimeError(f"Pilot failed for {target} at {energy:g} MeV")
    _write_rows(destination / "pilot_counts.csv", rows)
    return rows


def _selected_map(pilot_rows):
    return {(row["target"], float(row["energy_MeV"])): float(row["areal_mg_cm2"])
            for row in pilot_rows if int(row.get("selected", 0)) == 1}


def run_production(config, executable, destination, threads, targets, energies,
                   maximum_protons=None):
    pilot_path = destination / "pilot_counts.csv"
    with pilot_path.open(newline="", encoding="utf-8") as stream:
        selected = _selected_map(list(csv.DictReader(stream)))
    maximum_protons = maximum_protons or config["maximum_protons"]
    counts_path = destination / "thin_target_counts.csv"
    if counts_path.exists():
        with counts_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
    else:
        rows = []
    completed = {(row["target"], float(row["energy_MeV"])) for row in rows}
    logs = destination / "logs/production"
    logs.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ptcrysp-xsections-") as temporary:
        for target in targets:
            for energy in energies:
                if (target, energy) in completed:
                    continue
                residual_fields = [
                    f"n_{name.lower()}"
                    for name in config["target"][target]["active_residuals"]
                    if energy >= config["target"][target]["threshold_MeV"][name]
                ]
                thickness = selected[(target, energy)]
                batches = []
                total = 0
                batch_index = 0
                while total < maximum_protons:
                    n_batch = min(config["batch_protons"], maximum_protons - total)
                    output = Path(temporary) / "counts.csv"
                    if output.exists(): output.unlink()
                    seed = 2000000 + sum(ord(c) for c in target) * 10000 + int(energy * 10) * 100 + batch_index
                    row = invoke(executable, target, energy, thickness, n_batch, seed, threads,
                                 config["physics_list"], output,
                                 logs / f"{target}_{energy:g}_{batch_index:04d}.log")
                    batches.append(row)
                    total += n_batch
                    merged = merge_rows(batches)
                    if not residual_fields:
                        break
                    active_counts = [int(merged[field]) for field in residual_fields]
                    if active_counts and min(active_counts) >= config["minimum_residual_count"]:
                        break
                    batch_index += 1
                merged = merge_rows(batches)
                if not residual_fields:
                    merged["stopping_condition"] = "below_threshold"
                elif total < maximum_protons:
                    merged["stopping_condition"] = "residual_count"
                else:
                    merged["stopping_condition"] = "proton_cap"
                rows.append(merged)
                _write_rows(counts_path, rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("pilot", "production"))
    parser.add_argument("--config", type=Path, default=Path("config/xsections_scan.toml"))
    parser.add_argument("--executable", type=Path, default=Path("xsections_g4/build/thin_target"))
    parser.add_argument("--output-root", type=Path, default=Path("data/xsections/g4"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--targets", nargs="*")
    parser.add_argument("--energies", nargs="*", type=float)
    parser.add_argument("--maximum-protons", type=int)
    parser.add_argument("--force", action="store_true",
                        help="replace existing outputs with different runtime metadata")
    args = parser.parse_args()
    config = load(args.config)
    targets = args.targets or config["targets"]
    energies = args.energies or config["energies_MeV"]
    destination = (args.output_root / f"geant4-{config['geant4_version']}"
                   / config["physics_list"])
    destination.mkdir(parents=True, exist_ok=True)
    metadata = {key: value for key, value in config.items() if key != "target"}
    metadata["executable"] = str(args.executable)
    metadata["mode"] = args.mode
    metadata["runtime_targets"] = targets
    metadata["runtime_energies_MeV"] = energies
    metadata["runtime_threads"] = args.threads
    metadata["runtime_maximum_protons"] = args.maximum_protons or config["maximum_protons"]
    metadata_path = destination / "run_meta.json"
    counts_path = destination / "thin_target_counts.csv"
    if args.mode == "production" and args.force and counts_path.exists():
        counts_path.unlink()
    if metadata_path.exists() and args.mode == "production":
        previous = json.loads(metadata_path.read_text())
        comparable = ("config_sha256", "runtime_targets", "runtime_energies_MeV",
                      "runtime_maximum_protons")
        if (previous.get("mode") == "production"
                and any(previous.get(key) != metadata.get(key) for key in comparable)):
            if not args.force:
                raise RuntimeError("Existing scan metadata differ; use --force to replace the production output")
            if counts_path.exists(): counts_path.unlink()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    if args.mode == "pilot":
        rows = run_pilot(config, args.executable, destination, args.threads, targets, energies)
    else:
        rows = run_production(config, args.executable, destination, args.threads,
                              targets, energies, args.maximum_protons)
    print(f"Wrote {len(rows)} {args.mode} rows under {destination}")


if __name__ == "__main__":
    main()
