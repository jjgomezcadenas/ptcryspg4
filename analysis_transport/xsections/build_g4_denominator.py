"""Build effective Geant4 residual-production denominators and support audit."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .channels import CHANNELS
from .g4_cross_sections import estimate_factorized


def build(counts_path: Path, fit_dir: Path, destination: Path):
    counts = pd.read_csv(counts_path)
    destination.mkdir(parents=True, exist_ok=True)
    summaries = []
    for channel in CHANNELS:
        selected = counts[counts.target == channel.target].copy()
        residual_field = f"n_{channel.residual.lower()}"
        rows = []
        for point_id, row in enumerate(selected.itertuples(index=False)):
            count = int(getattr(row, residual_field))
            interactions = int(row.n_interactions)
            sigma, uncertainty, upper = estimate_factorized(
                count, interactions, float(row.sigma_inelastic_mb))
            rows.append({
                "point_id": point_id,
                "energy_MeV": row.energy_MeV,
                "sigma_inelastic_mb": row.sigma_inelastic_mb,
                "residual_probability": count / interactions,
                "sigma_mb": sigma,
                "sigma_unc_stat_mb": uncertainty,
                "sigma_upper_95_mb": upper,
                "residual_count": count,
                "n_interactions": interactions,
                "cross_section_data_sets": row.cross_section_data_sets,
                "model_counts": row.model_counts,
                "physics_list": row.physics_list,
                "geant4_version": row.geant4_version,
                "stopping_condition": row.stopping_condition,
            })
        curve = pd.DataFrame(rows).sort_values("energy_MeV")
        curve["sigma_unc_stat_mb"] = pd.to_numeric(
            curve.sigma_unc_stat_mb, errors="coerce")
        curve["sigma_upper_95_mb"] = pd.to_numeric(
            curve.sigma_upper_95_mb, errors="coerce")
        curve_path = destination / f"{channel.channel_id}_curve.csv"
        curve.to_csv(curve_path, index=False, lineterminator="\n")

        experimental = pd.read_csv(fit_dir / f"{channel.channel_id}_curve.csv")
        experimental_at_grid = np.interp(
            curve.energy_MeV, experimental.energy_MeV,
            experimental.sigma_nominal_mb)
        production = experimental_at_grid >= 0.05 * experimental.sigma_nominal_mb.max()
        zero_support = production & (curve.residual_count == 0)
        positive = curve.residual_count > 0
        relative_stat = np.divide(
            curve.sigma_unc_stat_mb.fillna(0.0), curve.sigma_mb,
            out=np.zeros(len(curve)), where=curve.sigma_mb > 0)
        summaries.append({
            "channel_id": channel.channel_id,
            "grid_points": len(curve),
            "experimental_production_grid_points": int(production.sum()),
            "zero_g4_support_in_production_region": int(zero_support.sum()),
            "zero_support_energies_MeV": ";".join(
                f"{energy:g}" for energy in curve.loc[
                    zero_support, "energy_MeV"]),
            "first_positive_g4_energy_MeV": (
                float(curve.loc[positive, "energy_MeV"].min())
                if positive.any() else ""),
            "total_pilot_residual_count": int(curve.residual_count.sum()),
            "median_relative_statistical_uncertainty_positive": (
                float(np.median(relative_stat[positive]))
                if positive.any() else ""),
            "maximum_relative_statistical_uncertainty_positive": (
                float(np.max(relative_stat[positive]))
                if positive.any() else ""),
        })
    summary = pd.DataFrame(summaries)
    summary.to_csv(destination / "support_summary.csv", index=False,
                   lineterminator="\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts", type=Path,
        default=Path("data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP/"
                     "denominator_pilot_counts.csv"))
    parser.add_argument("--fit-dir", type=Path,
                        default=Path("data/xsections/fits"))
    parser.add_argument(
        "--output", type=Path,
        default=Path("data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP/"
                     "denominator_pilot"))
    args = parser.parse_args()
    summary = build(args.counts, args.fit_dir, args.output)
    print(f"Built {len(summary)} pilot denominator curves and support audits")


if __name__ == "__main__":
    main()
