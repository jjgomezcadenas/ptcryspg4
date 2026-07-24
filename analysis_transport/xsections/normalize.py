"""Normalize all frozen external cross-section inputs."""

import argparse
import csv
import json
import shutil
from pathlib import Path

from . import endf6, exfor, iaea, jendl
from .channels import BY_PAIR, CHANNELS
from .common import write_points


def collect(repo: Path):
    raw = repo / "data/xsections/raw"
    datasets = []
    for path in sorted((raw / "exfor/2026-06-29/cx4").rglob("*.cx4")):
        dataset = exfor.read(path)
        if (dataset["target"], dataset["residual"]) in BY_PAIR:
            datasets.append(dataset)

    for channel in CHANNELS:
        target_file = raw / "jendl/4.0he" / f"{channel.target[0]}0{channel.target[1:]}.txt"
        datasets.append(jendl.read(target_file, channel.residual))
        endf_file = raw / "tendl/2023" / f"p-{channel.target[0]}0{channel.target[1:]}.tendl"
        datasets.append(endf6.read(endf_file, channel.target, channel.residual,
                                   channel.residual_z, channel.residual_a))
    datasets.append(iaea.read(raw / "iaea/medical/o6p13nt.txt"))
    return datasets


def normalize(repo: Path, clean=True):
    destination = repo / "data/xsections/normalized"
    if clean and destination.exists():
        for path in destination.iterdir():
            if path.name != ".gitkeep":
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
    destination.mkdir(parents=True, exist_ok=True)
    datasets = collect(repo)
    seen = set()
    catalog = []
    for dataset in datasets:
        dataset_id = dataset["dataset_id"]
        if dataset_id in seen:
            raise ValueError(f"Duplicate normalized dataset id: {dataset_id}")
        seen.add(dataset_id)
        csv_path = destination / f"{dataset_id}.csv"
        write_points(csv_path, dataset["points"])
        metadata = {key: value for key, value in dataset.items() if key != "points"}
        metadata["source_file"] = str(Path(metadata["source_file"]).relative_to(repo))
        metadata["point_count"] = len(dataset["points"])
        metadata_path = destination / f"{dataset_id}.json"
        metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        catalog.append({
            "dataset_id": dataset_id,
            "library": dataset["library"],
            "target": dataset["target"],
            "residual": dataset["residual"],
            "label": dataset["label"],
            "point_count": len(dataset["points"]),
            "point_file": str(csv_path.relative_to(repo)),
            "metadata_file": str(metadata_path.relative_to(repo)),
        })
    with (destination / "datasets.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=catalog[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(catalog)
    return catalog


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    catalog = normalize(args.repo.resolve())
    print(f"Normalized {len(catalog)} datasets into data/xsections/normalized")


if __name__ == "__main__":
    main()
