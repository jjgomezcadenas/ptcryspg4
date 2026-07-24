"""Plot EXFOR fits, replicas, distances and evaluated-curve comparisons."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .channels import CHANNELS
from .config import load
from .fit_exfor import load_exfor_points, load_pending_exfor_points


EVALUATIONS = {
    "JENDL-4.0/HE": ("#0072B2", "-"),
    "LANL ENDF/B-VII.1": ("#D55E00", "--"),
}


def _evaluation_curves(repo, channel):
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    selected = catalog[
        (catalog.target == channel.target)
        & (catalog.residual == channel.residual)
        & catalog.library.isin(EVALUATIONS)
    ]
    return [(row, pd.read_csv(repo / row["point_file"]))
            for row in selected.to_dict("records")]


def _fit_mask(points, threshold):
    uncertainty = points[["sigma_unc_minus_mb", "sigma_unc_plus_mb"]].mean(axis=1)
    return ((points.energy_MeV > threshold) & uncertainty.notna()
            & (uncertainty > 0))


def _comparison_metrics(channel, curve, evaluations):
    rows = []
    stable = curve.sigma_nominal_mb >= 0.05 * curve.sigma_nominal_mb.max()
    for metadata, evaluation in evaluations:
        common = (stable
                  & curve.energy_MeV.between(
                      evaluation.energy_MeV.min(), evaluation.energy_MeV.max()))
        energy = curve.loc[common, "energy_MeV"].to_numpy()
        nominal = curve.loc[common, "sigma_nominal_mb"].to_numpy()
        evaluated = np.interp(energy, evaluation.energy_MeV, evaluation.sigma_mb)
        ratio = evaluated / nominal
        rows.append({
            "channel_id": channel.channel_id,
            "evaluation": metadata["library"],
            "common_grid_points": len(energy),
            "median_ratio_to_exfor_fit": float(np.median(ratio)),
            "median_absolute_fractional_difference": float(
                np.median(np.abs(ratio - 1.0))),
            "maximum_absolute_fractional_difference": float(
                np.max(np.abs(ratio - 1.0))),
        })
    return rows


def _plot(repo, channel, config, figure_path):
    fits = repo / "data/xsections/fits"
    curve = pd.read_csv(fits / f"{channel.channel_id}_curve.csv")
    representatives = pd.read_csv(
        fits / f"{channel.channel_id}_representatives.csv")
    histogram = pd.read_csv(
        fits / f"{channel.channel_id}_distance_histogram.csv")
    points, _ = load_exfor_points(
        repo, channel, config["energy_min_MeV"], config["energy_max_MeV"])
    pending = load_pending_exfor_points(
        repo, channel, config["energy_min_MeV"], config["energy_max_MeV"])
    evaluations = _evaluation_curves(repo, channel)
    used = _fit_mask(points, config["threshold_MeV"][channel.channel_id])

    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.4))
    data_axis, replica_axis, distance_axis, comparison_axis = axes.flat

    fitted = points.loc[used]
    yerr = fitted[["sigma_unc_minus_mb", "sigma_unc_plus_mb"]].to_numpy().T
    data_axis.errorbar(fitted.energy_MeV, fitted.sigma_mb, yerr=yerr,
                       linestyle="none", marker="o", markersize=2.6,
                       color="0.25", alpha=0.62, elinewidth=0.5,
                       label="accepted EXFOR points", zorder=2)
    pending_styles = {
        "incident_energy_field_unclassified": (
            "s", "#CC79A7", "B0095 pending: energy field unclassified"),
        "shared_external_normalization_not_modelled": (
            "x", "#009E73", "Masuda shape comparison (published scale)"),
        "abundance_weighted_quantity": (
            "^", "#E69F00", "abundance-weighted pending series"),
    }
    if not pending.empty:
        for reason, rows in pending.groupby("reason_code"):
            marker, color, label = pending_styles.get(
                reason, ("d", "#999999", "other pending series"))
            data_axis.scatter(
                rows.energy_MeV, rows.sigma_mb, marker=marker, s=18,
                facecolors="none" if marker not in ("x", "+") else color,
                edgecolors=color if marker not in ("x", "+") else None,
                linewidths=0.8, alpha=0.8, label=label, zorder=2.5)
    data_axis.fill_between(curve.energy_MeV, curve.sigma_lower_16_mb,
                           curve.sigma_upper_84_mb, color="#56B4E9",
                           alpha=0.35, label="16--84% replicas", zorder=3)
    data_axis.plot(curve.energy_MeV, curve.sigma_nominal_mb, color="black",
                   linewidth=1.7, label="nominal median", zorder=4)
    sensitivity_path = fits / f"{channel.channel_id}_sensitivity.csv"
    sensitivity = None
    if sensitivity_path.exists():
        sensitivity = pd.read_csv(sensitivity_path)
        data_axis.plot(
            sensitivity.energy_MeV, sensitivity.sigma_sensitivity_mb,
            color="#CC79A7", linestyle="--", linewidth=1.5,
            label="fit including B0095", zorder=4)
    data_axis.set_title("EXFOR fit and uncertainty")
    data_axis.set_ylabel("Cross section (mb)")
    data_axis.legend(fontsize=7.3)

    colors = plt.cm.viridis(np.linspace(0.08, 0.92, 9))
    for color, (rank, rows) in zip(
            colors, representatives.groupby("representative_rank", sort=True)):
        quantile = rows.distance_quantile.iloc[0]
        replica_axis.plot(rows.energy_MeV, rows.sigma_mb, color=color,
                          linewidth=1.0, alpha=0.85,
                          label=f"q={quantile:.2f}")
    replica_axis.plot(curve.energy_MeV, curve.sigma_nominal_mb, color="black",
                      linewidth=1.8, label="nominal")
    replica_axis.set_title("Nine representative replicas")
    replica_axis.legend(fontsize=6.5, ncol=2)

    centres = 0.5 * (histogram.bin_left + histogram.bin_right)
    widths = histogram.bin_right - histogram.bin_left
    distance_axis.bar(centres, histogram["count"], width=widths,
                      color="#7A6FAC", edgecolor="white", linewidth=0.25)
    distance_axis.set_title("Replica distance from nominal")
    distance_axis.set_xlabel("D")
    distance_axis.set_ylabel("Replicas")

    comparison_axis.fill_between(
        curve.energy_MeV, curve.sigma_lower_16_mb, curve.sigma_upper_84_mb,
        color="#999999", alpha=0.25, label="EXFOR 16--84%")
    comparison_axis.plot(curve.energy_MeV, curve.sigma_nominal_mb,
                         color="black", linewidth=1.7, label="EXFOR nominal")
    if sensitivity is not None:
        comparison_axis.plot(
            sensitivity.energy_MeV, sensitivity.sigma_sensitivity_mb,
            color="#CC79A7", linestyle="--", linewidth=1.5,
            label="EXFOR fit including B0095")
    for metadata, evaluation in evaluations:
        color, linestyle = EVALUATIONS[metadata["library"]]
        comparison_axis.plot(evaluation.energy_MeV, evaluation.sigma_mb,
                             color=color, linestyle=linestyle, linewidth=1.5,
                             label=metadata["library"])
    comparison_axis.set_title("Evaluations as external comparisons")
    comparison_axis.set_xlabel("Proton energy (MeV)")
    comparison_axis.legend(fontsize=7.3)

    for axis in (data_axis, replica_axis, comparison_axis):
        axis.set_xlim(config["energy_min_MeV"], config["energy_max_MeV"])
        axis.set_ylim(bottom=0)
        axis.grid(alpha=0.18)
    distance_axis.set_xlim(left=0)
    distance_axis.set_ylim(bottom=0)
    distance_axis.grid(alpha=0.18)
    replica_axis.set_ylabel("Cross section (mb)")
    comparison_axis.set_ylabel("Cross section (mb)")
    figure.suptitle(channel.title, fontsize=13)
    figure.tight_layout()
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_path)
    plt.close(figure)
    return _comparison_metrics(channel, curve, evaluations)


def _write_generated(repo, summary, comparisons):
    output = repo / "docs/generated/xsection_fit"
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(comparisons).to_csv(
        output / "evaluation_comparison.csv", index=False, lineterminator="\n")
    thresholds = pd.read_csv(repo / "data/xsections/reaction_thresholds.csv")
    thresholds.to_csv(
        output / "reaction_thresholds.csv", index=False, lineterminator="\n")
    with (output / "reaction_thresholds.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrr}\n\\toprule\n")
        stream.write("Lowest channel & $Q$ (MeV) & $E_{\\mathrm{th,lab}}$ (MeV) \\\\\n")
        stream.write("\\midrule\n")
        latex_reactions = {
            "p_O16_x_O15": "$^{16}\\mathrm{O}(p,d)^{15}\\mathrm{O}$",
            "p_O16_x_C11": (
                "$^{16}\\mathrm{O}(p,d+\\alpha)^{11}\\mathrm{C}$"),
        }
        for row in thresholds.to_dict("records"):
            reaction = latex_reactions[row["channel_id"]]
            stream.write(
                f"{reaction} & {row['q_value_MeV']:.3f} & "
                f"{row['laboratory_threshold_MeV']:.3f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")

    sensitivity_path = repo / "data/xsections/fits/sensitivity_summary.csv"
    if sensitivity_path.exists():
        sensitivity = pd.read_csv(sensitivity_path)
        sensitivity.to_csv(
            output / "curation_sensitivity.csv", index=False,
            lineterminator="\n")

    distance_rows = []
    for channel in CHANNELS:
        replicas = pd.read_csv(
            repo / "data/xsections/fits"
            / f"{channel.channel_id}_replicas.csv")
        grouped = replicas.groupby("smoothing_lambda").distance_D.agg(
            ["count", "min", "median", "max"]).reset_index()
        grouped.insert(0, "channel_id", channel.channel_id)
        distance_rows.extend(grouped.to_dict("records"))
    pd.DataFrame(distance_rows).to_csv(
        output / "distance_by_smoothing.csv", index=False,
        lineterminator="\n")
    columns = [
        "channel_id", "points_in_range", "points_used",
        "independent_campaigns_used", "campaign_fractional_spread",
        "peak_energy_MeV", "peak_sigma_mb", "peak_lower_16_mb",
        "peak_upper_84_mb", "median_relative_half_width_production",
        "p90_relative_half_width_production", "chi2_per_dof",
    ]
    summary[columns].to_csv(
        output / "fit_summary.csv", index=False, lineterminator="\n")
    curation = pd.read_csv(repo / "data/xsections/curation.csv")
    curation_rows = []
    for channel in CHANNELS:
        selected = curation[curation.channel_id == channel.channel_id]
        counts = selected.state.value_counts().to_dict()
        curation_rows.append({
            "channel_id": channel.channel_id,
            "accepted_series": counts.get("accepted", 0),
            "pending_series": counts.get("pending", 0),
            "excluded_series": counts.get("excluded", 0),
            "fit_points": int(selected.loc[
                selected.state == "accepted", "points_passing_point_rules"].sum()),
        })
    pd.DataFrame(curation_rows).to_csv(
        output / "curation_summary.csv", index=False, lineterminator="\n")
    with (output / "curation_summary.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
        stream.write("Channel & Accepted & Pending & Excluded & Fit points \\\\\n")
        stream.write("\\midrule\n")
        for row in curation_rows:
            title = next(c.title for c in CHANNELS if c.channel_id == row["channel_id"])
            stream.write(
                f"{title} & {row['accepted_series']} & {row['pending_series']} & "
                f"{row['excluded_series']} & {row['fit_points']} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")

    point_curation = pd.read_csv(repo / "data/xsections/point_curation.csv")
    energy_limits = point_curation.groupby("dataset_id").energy_MeV.agg(
        ["min", "max"])
    with (output / "curation_detail.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{longtable}{llllrr}\n\\toprule\n")
        stream.write("\\caption{Dataset-level EXFOR curation audit.}"
                     "\\label{tab:curation-detail}\\\\\n")
        stream.write("Accession & Channel & Energy (MeV) & Decision & Points & Point-rule pass \\\\\n")
        stream.write("\\midrule\n\\endfirsthead\n\\toprule\n")
        stream.write("Accession & Channel & Energy (MeV) & Decision & Points & Point-rule pass \\\\\n")
        stream.write("\\midrule\n\\endhead\n")
        for row in curation.to_dict("records"):
            title = next(c.title for c in CHANNELS if c.channel_id == row["channel_id"])
            limits = energy_limits.loc[row["dataset_id"]]
            energy = f"{limits['min']:.2f}--{limits['max']:.2f}"
            decision = row["state"].replace("pending", "pending verification")
            stream.write(
                f"{row['accession']} & {title} & {energy} & {decision} & "
                f"{row['points_total']} & {row['points_passing_point_rules']} \\\\\n")
        stream.write("\\bottomrule\n\\end{longtable}\n")

    point_categories = [
        ("accepted", "Used"),
        ("no_positive_reported_uncertainty", "No quoted uncertainty"),
        ("outside_energy_range", "Outside range"),
        ("at_or_below_threshold", "At/below threshold"),
        ("nonpositive_cross_section", "Nonpositive"),
        ("dataset_pending", "Dataset pending"),
    ]
    dataset_channel = curation.set_index("dataset_id").channel_id
    point_curation["channel_id"] = point_curation.dataset_id.map(dataset_channel)
    with (output / "point_curation_summary.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
        stream.write("Channel & Used & No error & Range & Threshold & Nonpositive & Pending \\\\\n")
        stream.write("\\midrule\n")
        for channel in CHANNELS:
            selected = point_curation[point_curation.channel_id == channel.channel_id]
            counts = selected.reason_code.value_counts().to_dict()
            values = [counts.get(code, 0) for code, _ in point_categories]
            stream.write(
                f"{channel.title} & " + " & ".join(str(value) for value in values)
                + " \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")
    with (output / "fit_summary.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
        stream.write("Channel & $N_{\\mathrm{camp}}$ & $N_{\\mathrm{pt}}$ & "
                     "$e^{s_c}-1$ & Median $u$ & Peak (mb) & "
                     "$\\chi^2/\\nu$ \\\\\n")
        stream.write("\\midrule\n")
        for row in summary.to_dict("records"):
            title = next(c.title for c in CHANNELS if c.channel_id == row["channel_id"])
            peak = (f"{row['peak_sigma_mb']:.1f} "
                    f"[{row['peak_lower_16_mb']:.1f}, {row['peak_upper_84_mb']:.1f}]")
            stream.write(
                f"{title} & {row['independent_campaigns_used']} & {row['points_used']} & "
                f"{100.0 * row['campaign_fractional_spread']:.0f}\\% & "
                f"{100.0 * row['median_relative_half_width_production']:.0f}\\% & {peak} & "
                f"{row['chi2_per_dof']:.2f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")


def generate(repo: Path, config_path: Path):
    config = load(config_path)
    summary = pd.read_csv(repo / "data/xsections/fits/fit_summary.csv")
    comparisons = []
    for channel in CHANNELS:
        comparisons.extend(_plot(
            repo, channel, config,
            repo / "docs/figures/xsection_fit" / f"{channel.channel_id}.pdf"))
    _write_generated(repo, summary, comparisons)
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--config", type=Path,
                        default=Path("config/xsection_fit.toml"))
    args = parser.parse_args()
    comparisons = generate(args.repo.resolve(), args.config.resolve())
    print(f"Generated five EXFOR-fit figures and {len(comparisons)} evaluation comparisons")


if __name__ == "__main__":
    main()
