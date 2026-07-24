"""Convert aggregate Geant4 thin-target counters into production curves."""

import argparse
import csv
import math
from pathlib import Path

from .channels import CHANNELS


ZERO_COUNT_95_UPPER = -math.log(0.05)


def estimate(count, protons, nuclear_areal_density):
    scale = 1.0e27 / (protons * nuclear_areal_density)
    sigma = count * scale
    uncertainty = math.sqrt(count) * scale if count else ""
    upper = ZERO_COUNT_95_UPPER * scale if count == 0 else ""
    return sigma, uncertainty, upper


def estimate_factorized(count, interactions, sigma_inelastic_mb):
    probability = count / interactions
    sigma = sigma_inelastic_mb * probability
    uncertainty = (sigma_inelastic_mb * math.sqrt(
        probability * (1.0 - probability) / interactions)
        if count else "")
    upper = (sigma_inelastic_mb * ZERO_COUNT_95_UPPER / interactions
             if count == 0 else "")
    return sigma, uncertainty, upper


def convert(counts_path: Path, destination: Path):
    with counts_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    destination.mkdir(parents=True, exist_ok=True)
    outputs = []
    for channel in CHANNELS:
        points = []
        for row in rows:
            if row["target"] != channel.target:
                continue
            count = int(row[f"n_{channel.residual.lower()}"])
            sigma, uncertainty, upper = estimate(
                count, int(row["n_protons"]), float(row["nuclear_areal_density_cm2"]))
            points.append({
                "point_id": len(points),
                "energy_MeV": float(row["energy_MeV"]),
                "sigma_mb": sigma,
                "sigma_unc_stat_mb": uncertainty,
                "sigma_upper_95_mb": upper,
                "residual_count": count,
                "n_protons": int(row["n_protons"]),
                "areal_mg_cm2": float(row["areal_mg_cm2"]),
                "mean_interaction_energy_MeV": float(row[f"mean_{channel.residual.lower()}_energy_MeV"]),
                "mean_continuous_loss_MeV": float(row["mean_continuous_loss_MeV"]),
                "inelastic_probability": int(row["n_inelastic"]) / int(row["n_protons"]),
                "stopping_condition": row.get("stopping_condition", ""),
            })
        if not points:
            continue
        points.sort(key=lambda point: point["energy_MeV"])
        path = destination / f"g4_qgsp_bic_hp_{channel.channel_id}.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=points[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(points)
        outputs.append(path)
    for target in sorted({row["target"] for row in rows}):
        points = []
        for row in rows:
            if row["target"] != target:
                continue
            count = int(row["n_inelastic"])
            sigma, uncertainty, upper = estimate(
                count, int(row["n_protons"]), float(row["nuclear_areal_density_cm2"]))
            points.append({
                "point_id": len(points), "energy_MeV": float(row["energy_MeV"]),
                "sigma_mb": sigma, "sigma_unc_stat_mb": uncertainty,
                "sigma_upper_95_mb": upper, "interaction_count": count,
                "n_protons": int(row["n_protons"]),
            })
        points.sort(key=lambda point: point["energy_MeV"])
        path = destination / f"g4_qgsp_bic_hp_p_{target}_inelastic.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=points[0].keys(), lineterminator="\n")
            writer.writeheader(); writer.writerows(points)
        outputs.append(path)
    catalog_path = destination / "models.csv"
    fieldnames = ("model_id", "channel_id", "curve_file", "parent_dataset_ids",
                  "construction", "energy_min_MeV", "energy_max_MeV",
                  "interpolation", "extrapolation", "version")
    existing = []
    if catalog_path.exists():
        with catalog_path.open(newline="", encoding="utf-8") as stream:
            existing = [row for row in csv.DictReader(stream)
                        if not row["model_id"].startswith("g4_qgsp_bic_hp_")]
    for channel in CHANNELS:
        path = destination / f"g4_qgsp_bic_hp_{channel.channel_id}.csv"
        if not path.exists(): continue
        with path.open(newline="", encoding="utf-8") as stream:
            points = list(csv.DictReader(stream))
        existing.append({
            "model_id": f"g4_qgsp_bic_hp_{channel.channel_id}",
            "channel_id": channel.channel_id,
            "curve_file": str(path), "parent_dataset_ids": "",
            "construction": "Geant4 11.4.1 isotopically pure thin-target residual count",
            "energy_min_MeV": min(float(row["energy_MeV"]) for row in points),
            "energy_max_MeV": max(float(row["energy_MeV"]) for row in points),
            "interpolation": "linear", "extrapolation": "none", "version": "validation-1",
        })
    with catalog_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader(); writer.writerows(existing)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path,
                        default=Path("data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP/thin_target_counts.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/xsections/models"))
    args = parser.parse_args()
    outputs = convert(args.counts, args.output)
    print(f"Wrote {len(outputs)} Geant4 cross-section curves")


if __name__ == "__main__":
    main()
