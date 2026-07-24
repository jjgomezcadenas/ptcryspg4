"""Generate cross-section comparison figures, tables and metrics."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .channels import CHANNELS


EVALUATION_STYLES = {
    "JENDL-4.0/HE": ("#0072B2", "-"),
    "LANL ENDF/B-VII.1": ("#D55E00", "--"),
    "IAEA recommended": ("#009E73", "-."),
}


def load_catalog(repo: Path):
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    datasets = []
    for row in catalog.to_dict("records"):
        points = pd.read_csv(repo / row["point_file"])
        datasets.append((row, points))
    return datasets


def _g4_curve(repo, channel_id):
    path = repo / "data/xsections/models" / f"g4_qgsp_bic_hp_{channel_id}.csv"
    return pd.read_csv(path) if path.exists() else None


def _coverage(datasets):
    rows = []
    for channel in CHANNELS:
        selected = [(metadata, points) for metadata, points in datasets
                    if metadata["target"] == channel.target and metadata["residual"] == channel.residual]
        for library in sorted({metadata["library"] for metadata, _ in selected}):
            members = [(metadata, points) for metadata, points in selected if metadata["library"] == library]
            energies = np.concatenate([points.energy_MeV.to_numpy() for _, points in members])
            rows.append({"channel_id": channel.channel_id, "source": library,
                         "datasets": len(members),
                         "points": sum(len(points) for _, points in members),
                         "energy_min_MeV": energies.min(), "energy_max_MeV": energies.max()})
    return rows


def _metrics(repo, datasets, ratio_fraction):
    rows = []
    for channel in CHANNELS:
        g4 = _g4_curve(repo, channel.channel_id)
        selected = [(metadata, points) for metadata, points in datasets
                    if metadata["target"] == channel.target and metadata["residual"] == channel.residual]
        for metadata, points in selected:
            peak = points.iloc[points.sigma_mb.to_numpy().argmax()]
            row = {"channel_id": channel.channel_id, "dataset_id": metadata["dataset_id"],
                   "source": metadata["library"], "peak_sigma_mb": peak.sigma_mb,
                   "peak_energy_MeV": peak.energy_MeV, "median_ratio_to_g4": "",
                   "common_points": 0}
            if g4 is not None and len(g4) >= 2 and g4.sigma_mb.max() > 0:
                mask = ((points.energy_MeV >= g4.energy_MeV.min())
                        & (points.energy_MeV <= g4.energy_MeV.max()))
                comparison = points.loc[mask]
                if len(comparison):
                    denominator = np.interp(comparison.energy_MeV, g4.energy_MeV, g4.sigma_mb)
                    stable = denominator >= ratio_fraction * g4.sigma_mb.max()
                    ratios = comparison.sigma_mb.to_numpy()[stable] / denominator[stable]
                    if len(ratios):
                        row["median_ratio_to_g4"] = float(np.median(ratios))
                        row["common_points"] = len(ratios)
            rows.append(row)
    return rows


def _plot_channel(repo, channel, datasets, output, ratio_fraction):
    selected = [(metadata, points) for metadata, points in datasets
                if metadata["target"] == channel.target and metadata["residual"] == channel.residual]
    g4 = _g4_curve(repo, channel.channel_id)
    figure, (axis, ratio_axis) = plt.subplots(
        2, 1, figsize=(7.2, 6.2), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    exfor_label = True
    markers = ("o", "s", "^", "v", "D", "P", "X", "<", ">")
    exfor_index = 0
    for metadata, points in selected:
        if metadata["library"] == "EXFOR":
            yerr = points[["sigma_unc_minus_mb", "sigma_unc_plus_mb"]].to_numpy().T
            xerr = points[["energy_unc_minus_MeV", "energy_unc_plus_MeV"]].to_numpy().T
            axis.errorbar(points.energy_MeV, points.sigma_mb, yerr=yerr, xerr=xerr,
                          linestyle="none", marker=markers[exfor_index % len(markers)],
                          markersize=3.2, alpha=0.52, color="0.35",
                          elinewidth=0.55, label="EXFOR campaigns" if exfor_label else None)
            exfor_label = False; exfor_index += 1
        else:
            color, line_style = EVALUATION_STYLES[metadata["library"]]
            axis.plot(points.energy_MeV, points.sigma_mb, color=color,
                      linestyle=line_style, linewidth=1.7, label=metadata["library"])
            if g4 is not None and len(g4) >= 2:
                mask = ((points.energy_MeV >= g4.energy_MeV.min())
                        & (points.energy_MeV <= g4.energy_MeV.max()))
                energy = points.energy_MeV.to_numpy()[mask]
                denominator = np.interp(energy, g4.energy_MeV, g4.sigma_mb)
                stable = denominator >= ratio_fraction * g4.sigma_mb.max()
                ratio_axis.plot(energy[stable], points.sigma_mb.to_numpy()[mask][stable] / denominator[stable],
                                color=color, linestyle=line_style, linewidth=1.4)
    if g4 is not None:
        axis.errorbar(g4.energy_MeV, g4.sigma_mb, yerr=g4.sigma_unc_stat_mb,
                      color="black", marker=".", linewidth=1.6, label="Geant4 QGSP_BIC_HP")
    axis.set_ylabel("Production cross section (mb)")
    axis.set_title(channel.title)
    axis.set_xlim(5, 120)
    axis.set_ylim(bottom=0)
    axis.grid(alpha=0.2)
    axis.legend(fontsize=8, ncol=2)
    ratio_axis.axhline(1.0, color="black", linewidth=0.7)
    ratio_axis.set_ylabel("Eval./G4")
    ratio_axis.set_xlabel("Proton energy (MeV)")
    ratio_axis.grid(alpha=0.2)
    if g4 is None or len(g4) < 2:
        ratio_axis.text(0.5, 0.5, "Geant4 production scan pending",
                        transform=ratio_axis.transAxes, ha="center", va="center")
        ratio_axis.set_ylim(0, 2)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)


def _write_tex_tables(output, coverage, metrics):
    output.mkdir(parents=True, exist_ok=True)
    with (output / "coverage_table.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{llrrcc}\n\\toprule\nChannel & Source & Curves & Points & $E_{min}$ & $E_{max}$ \\\\\n\\midrule\n")
        for row in coverage:
            channel = row["channel_id"].replace("_", "\\_")
            source = row["source"].replace("_", "\\_")
            stream.write(f"{channel} & {source} & {row['datasets']} & {row['points']} & "
                         f"{row['energy_min_MeV']:.2f} & {row['energy_max_MeV']:.2f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")
    summaries = []
    for channel in CHANNELS:
        members = [row for row in metrics if row["channel_id"] == channel.channel_id]
        summaries.append((channel, len(members), sum(row["common_points"] for row in members)))
    with (output / "comparison_metrics.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrr}\n\\toprule\nChannel & Input curves & Ratio points \\\\\n\\midrule\n")
        for channel, curves, points in summaries:
            stream.write(f"{channel.title} & {curves} & {points} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")


def _channel_summaries(repo):
    summaries = []
    for channel in CHANNELS:
        curve = _g4_curve(repo, channel.channel_id)
        total_path = repo / "data/xsections/models" / f"g4_qgsp_bic_hp_p_{channel.target}_inelastic.csv"
        maximum_fraction = ""
        if curve is not None and total_path.exists():
            total = pd.read_csv(total_path)
            denominator = np.interp(curve.energy_MeV, total.energy_MeV, total.sigma_mb)
            valid = denominator > 0
            if valid.any(): maximum_fraction = float(np.max(curve.sigma_mb.to_numpy()[valid] / denominator[valid]))
        summaries.append({"channel_id": channel.channel_id,
                          "g4_points": 0 if curve is None else len(curve),
                          "maximum_fraction_of_total_inelastic": maximum_fraction})
    return summaries


def _write_scan_tables(repo, output):
    pilot_path = repo / "data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP/pilot_counts.csv"
    counts_path = pilot_path.parent / "thin_target_counts.csv"
    if pilot_path.exists():
        pilot = pd.read_csv(pilot_path)
        selected = pilot[pilot.selected == 1]
        thickness = selected.pivot(index="energy_MeV", columns="target", values="areal_mg_cm2")
        with (output / "selected_thickness.tex").open("w", encoding="utf-8") as stream:
            stream.write("\\begin{longtable}{rrrr}\n\\toprule\n")
            stream.write("Energy (MeV) & C-12 & N-14 & O-16 \\\\\n\\midrule\n\\endfirsthead\n")
            stream.write("\\toprule\nEnergy (MeV) & C-12 & N-14 & O-16 \\\\\n\\midrule\n\\endhead\n")
            for energy, row in thickness.iterrows():
                stream.write(f"{energy:g} & {row['C12']:g} & {row['N14']:g} & {row['O16']:g} \\\\\n")
            stream.write("\\bottomrule\n\\end{longtable}\n")
    if counts_path.exists():
        counts = pd.read_csv(counts_path)
        summary = []
        for target, rows in counts.groupby("target"):
            summary.append({
                "target": target,
                "points": len(rows),
                "maximum_mean_energy_loss_MeV": rows.mean_continuous_loss_MeV.max(),
                "maximum_inelastic_probability": (rows.n_inelastic / rows.n_protons).max(),
                "incident_protons": int(rows.n_protons.sum()),
                "C11_residuals": int(rows.n_c11.sum()),
                "N13_residuals": int(rows.n_n13.sum()),
                "O15_residuals": int(rows.n_o15.sum()),
            })
        with (output / "scan_summary.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=summary[0].keys(), lineterminator="\n")
            writer.writeheader(); writer.writerows(summary)


def generate(repo: Path, ratio_fraction=0.05):
    datasets = load_catalog(repo)
    coverage = _coverage(datasets)
    metrics = _metrics(repo, datasets, ratio_fraction)
    figure_dir = repo / "docs/figures/xsections"
    for channel in CHANNELS:
        _plot_channel(repo, channel, datasets, figure_dir / f"{channel.channel_id}.pdf", ratio_fraction)
    generated = repo / "docs/generated/xsections"
    _write_tex_tables(generated, coverage, metrics)
    _write_scan_tables(repo, generated)
    channel_summaries = _channel_summaries(repo)
    with (generated / "coverage.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=coverage[0].keys(), lineterminator="\n"); writer.writeheader(); writer.writerows(coverage)
    with (generated / "comparison_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=metrics[0].keys(), lineterminator="\n"); writer.writeheader(); writer.writerows(metrics)
    with (generated / "channel_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=channel_summaries[0].keys(), lineterminator="\n")
        writer.writeheader(); writer.writerows(channel_summaries)
    metadata = {"ratio_minimum_peak_fraction": ratio_fraction,
                "normalized_dataset_count": len(datasets), "channel_count": len(CHANNELS)}
    (generated / "comparison_meta.json").write_text(json.dumps(metadata, indent=2) + "\n")
    g4_curves = sum(_g4_curve(repo, channel.channel_id) is not None for channel in CHANNELS)
    with (generated / "scan_status.tex").open("w", encoding="utf-8") as stream:
        if g4_curves == len(CHANNELS):
            run_meta_path = repo / "data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP/run_meta.json"
            run_meta = json.loads(run_meta_path.read_text()) if run_meta_path.exists() else {}
            runtime_cap = run_meta.get("runtime_maximum_protons")
            configured_cap = run_meta.get("maximum_protons")
            if runtime_cap and configured_cap and runtime_cap < configured_cap:
                stream.write(
                    f"The stored all-grid validation scan supplies all five channel curves with a runtime cap "
                    f"of \\num{{{runtime_cap}}} protons per point.  These preliminary curves validate the "
                    "pipeline; the precision scan uses the configured residual-count stopping rule.\n")
            else:
                stream.write("The stored Geant4 scan supplies all five channel curves.\n")
        elif g4_curves:
            stream.write(f"The stored pilot scan currently supplies {g4_curves} of the five channel curves.\n")
        else:
            stream.write("The external-data comparison is complete; the Geant4 production scan is pending.\n")
    return coverage, metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--ratio-minimum-peak-fraction", type=float, default=0.05)
    args = parser.parse_args()
    coverage, metrics = generate(args.repo.resolve(), args.ratio_minimum_peak_fraction)
    print(f"Generated five figures from {len(metrics)} source curves and {len(coverage)} coverage groups")


if __name__ == "__main__":
    main()
