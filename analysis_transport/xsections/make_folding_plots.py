#!/usr/bin/env python3
"""Generate the four diagnostic figures for cross-section folding outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROFILE_ORDER = ("C11", "O15", "N13", "all_production", "all_d120s300")
PROFILE_TITLES = {
    "C11": "C-11 production",
    "O15": "O-15 production",
    "N13": "N-13 production",
    "all_production": "All production",
    "all_d120s300": "In-room measured decays",
}
COLORS = {
    "C11": "#0072B2",
    "O15": "#D55E00",
    "N13": "#009E73",
    "all_production": "#7A6FAC",
    "all_d120s300": "#CC79A7",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_columns(frame: pd.DataFrame, fields, label: str) -> None:
    missing = sorted(set(fields) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(missing)}")


def _load_products(folding_directory: Path):
    paths = {
        "nominal": folding_directory / "nominal_isotope_profiles.csv",
        "bands": folding_directory / "profile_bands.csv",
        "summary": folding_directory / "production_summary.csv",
        "uncertainty": folding_directory / "uncertainty_summary.csv",
    }
    missing = [path.name for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "folding output is incomplete: " + ", ".join(missing))
    nominal = pd.read_csv(paths["nominal"])
    bands = pd.read_csv(paths["bands"])
    summary = pd.read_csv(paths["summary"])
    uncertainty = pd.read_csv(paths["uncertainty"])
    _require_columns(
        nominal,
        ("profile_label", "depth_mm", "expected_count_run"),
        "nominal profiles",
    )
    _require_columns(
        bands,
        ("profile_label", "depth_mm", "nominal_run", "q16_run", "q84_run"),
        "profile bands",
    )
    _require_columns(
        summary,
        (
            "model", "replica_id", "profile_label", "expected_count_run",
            "R50_prod_mm", "R50_shift_mm",
        ),
        "production summary",
    )
    _require_columns(
        uncertainty,
        (
            "profile_label", "nominal_yield_run", "yield_half_width_run",
            "R50_shift_half_width_mm",
        ),
        "uncertainty summary",
    )
    labels = tuple(label for label in PROFILE_ORDER if label in set(nominal.profile_label))
    if labels != PROFILE_ORDER:
        absent = sorted(set(PROFILE_ORDER) - set(labels))
        raise ValueError("folding profiles missing labels: " + ", ".join(absent))
    for label in labels:
        rows = bands.loc[bands.profile_label == label].sort_values("depth_mm")
        if np.any(rows.q16_run > rows.q84_run):
            raise ValueError(f"{label}: q16 exceeds q84 in profile band")
        if not np.all(np.diff(rows.depth_mm) > 0):
            raise ValueError(f"{label}: profile depths must increase")
    return paths, nominal, bands, summary, uncertainty


def _figure_title(title: str, context_label: str | None) -> str:
    return f"{title}\n{context_label}" if context_label else title


def _save_figure(figure, output: Path) -> None:
    figure.savefig(
        output,
        metadata={
            "Creator": "ptcryspg4 make_folding_plots.py",
            "CreationDate": None,
            "ModDate": None,
        },
    )


def _profiles_figure(
    bands: pd.DataFrame,
    summary: pd.DataFrame,
    output: Path,
    context_label: str | None,
) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(10.6, 6.4), sharex=True)
    axes = axes.flat
    for axis, label in zip(axes, PROFILE_ORDER, strict=False):
        rows = bands.loc[bands.profile_label == label].sort_values("depth_mm")
        depth = rows.depth_mm.to_numpy(dtype=float)
        nominal = rows.nominal_run.to_numpy(dtype=float)
        scale = float(np.max(nominal))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"{label}: nominal profile has no positive support")
        axis.fill_between(
            depth,
            rows.q16_run.to_numpy(dtype=float) / scale,
            rows.q84_run.to_numpy(dtype=float) / scale,
            color=COLORS[label],
            alpha=0.27,
            label="replica 16--84%",
        )
        axis.plot(
            depth, nominal / scale, color=COLORS[label], linewidth=1.8,
            label="nominal")
        nominal_row = summary.loc[
            (summary.model == "nominal") & (summary.profile_label == label)]
        if len(nominal_row) != 1:
            raise ValueError(f"{label}: expected one nominal summary row")
        edge = float(nominal_row.R50_prod_mm.iloc[0])
        if np.isfinite(edge):
            axis.axvline(edge, color="black", linestyle="--", linewidth=1.0,
                        label=f"R50 = {edge:.2f} mm")
        axis.set_title(PROFILE_TITLES[label], fontsize=10)
        axis.set_ylim(bottom=0)
        axis.grid(alpha=0.18)
        axis.legend(fontsize=7, loc="upper right")
    axes[-1].axis("off")
    for axis in axes[:5]:
        axis.set_xlabel("Depth (mm)")
        axis.set_ylabel("Profile / nominal maximum")
    figure.suptitle(_figure_title(
        "Production-depth profiles and cross-section uncertainty", context_label))
    figure.tight_layout()
    _save_figure(figure, output)
    plt.close(figure)


def _distribution_figure(
    summary: pd.DataFrame,
    output: Path,
    *,
    quantity: str,
    context_label: str | None,
) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(10.6, 6.2))
    axes = axes.flat
    for axis, label in zip(axes, PROFILE_ORDER, strict=False):
        nominal = summary.loc[
            (summary.model == "nominal") & (summary.profile_label == label)]
        replicas = summary.loc[
            (summary.model == "replica") & (summary.profile_label == label)]
        if len(nominal) != 1 or replicas.empty:
            raise ValueError(f"{label}: incomplete nominal or replica summary")
        if quantity == "R50":
            values = replicas.R50_shift_mm.to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            xlabel = "R50(replica) - R50(nominal) (mm)"
            reference = 0.0
        elif quantity == "yield":
            denominator = float(nominal.expected_count_run.iloc[0])
            if denominator <= 0:
                raise ValueError(f"{label}: nominal yield must be positive")
            values = replicas.expected_count_run.to_numpy(dtype=float) / denominator
            xlabel = "Replica yield / nominal yield"
            reference = 1.0
        else:
            raise ValueError(f"unknown distribution quantity '{quantity}'")
        if len(values) == 0:
            raise ValueError(f"{label}: distribution contains no finite values")
        bins = min(35, max(8, int(np.sqrt(len(values)))))
        display_low, display_high = np.quantile(values, [0.005, 0.995])
        if display_high <= display_low:
            display_low = float(np.min(values)) - 0.5
            display_high = float(np.max(values)) + 0.5
        padding = 0.04 * (display_high - display_low)
        display_low -= padding
        display_high += padding
        outside = int(np.count_nonzero(
            (values < display_low) | (values > display_high)))
        axis.hist(values, bins=bins, range=(display_low, display_high),
                  color=COLORS[label], alpha=0.82,
                  edgecolor="white", linewidth=0.35)
        q16, q50, q84 = np.quantile(values, [0.16, 0.50, 0.84])
        axis.axvline(reference, color="black", linewidth=0.8)
        axis.axvline(q50, color="black", linestyle="--", linewidth=1.1,
                    label=f"q50={q50:.3g}")
        axis.axvspan(q16, q84, color="black", alpha=0.08,
                    label=f"16--84%: {q16:.3g}, {q84:.3g}")
        axis.set_title(PROFILE_TITLES[label], fontsize=10)
        axis.set_xlabel(xlabel)
        axis.set_ylabel("Replicas")
        axis.grid(alpha=0.18)
        axis.legend(fontsize=7)
        if outside:
            axis.text(
                0.02, 0.96, f"{outside} replicas outside display range",
                transform=axis.transAxes, ha="left", va="top", fontsize=6.5)
    axes[-1].axis("off")
    title = (
        "Replica-induced production-edge displacement"
        if quantity == "R50" else "Replica yield relative to nominal")
    figure.suptitle(_figure_title(title, context_label))
    figure.tight_layout()
    _save_figure(figure, output)
    plt.close(figure)


def _convergence_figure(
    convergence: pd.DataFrame,
    output: Path,
    context_label: str | None,
) -> None:
    fields = (
        "candidate_width_MeV", "profile_label",
        "max_paired_yield_change_run", "yield_replica_half_width_run",
        "max_paired_R50_change_mm", "R50_replica_half_width_mm",
    )
    _require_columns(convergence, fields, "convergence table")
    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.7))
    for label in PROFILE_ORDER:
        rows = convergence.loc[
            convergence.profile_label == label].sort_values("candidate_width_MeV")
        if rows.empty:
            raise ValueError(f"convergence table missing profile '{label}'")
        width = rows.candidate_width_MeV.to_numpy(dtype=float)
        yield_denominator = rows.yield_replica_half_width_run.to_numpy(dtype=float)
        r50_denominator = rows.R50_replica_half_width_mm.to_numpy(dtype=float)
        yield_ratio = np.divide(
            rows.max_paired_yield_change_run.to_numpy(dtype=float),
            yield_denominator,
            out=np.full(len(rows), np.nan),
            where=yield_denominator > 0,
        )
        r50_ratio = np.divide(
            rows.max_paired_R50_change_mm.to_numpy(dtype=float),
            r50_denominator,
            out=np.full(len(rows), np.nan),
            where=r50_denominator > 0,
        )
        axes[0].plot(width, yield_ratio, marker="o", markersize=3.5,
                     color=COLORS[label], label=PROFILE_TITLES[label])
        axes[1].plot(width, r50_ratio, marker="o", markersize=3.5,
                     color=COLORS[label])
        axes[2].plot(
            width, rows.max_paired_R50_change_mm, marker="o", markersize=3.5,
            color=COLORS[label])
    axes[0].axhline(0.10, color="black", linestyle="--", linewidth=1.0)
    axes[1].axhline(0.10, color="black", linestyle="--", linewidth=1.0)
    axes[2].axhline(0.10, color="black", linestyle="--", linewidth=1.0)
    axes[1].set_yscale("symlog", linthresh=0.05, linscale=0.8)
    axes[2].set_yscale("symlog", linthresh=0.01, linscale=0.8)
    axes[0].set_ylabel("Maximum yield change / replica half-width")
    axes[1].set_ylabel("Maximum R50 change / replica half-width")
    axes[2].set_ylabel("Maximum R50 change (mm)")
    for axis in axes:
        axis.set_xlabel("Candidate energy-bin width (MeV)")
        axis.set_ylim(bottom=0)
        axis.grid(alpha=0.18)
    axes[0].legend(fontsize=6.6)
    figure.suptitle(_figure_title(
        "Offline energy-grid convergence", context_label))
    figure.tight_layout()
    _save_figure(figure, output)
    plt.close(figure)


def _quantitative_summary(
    summary: pd.DataFrame,
    uncertainty: pd.DataFrame,
) -> pd.DataFrame:
    nominal = summary.loc[summary.model == "nominal", [
        "profile_label", "expected_count_run", "R50_prod_mm"]]
    result = nominal.merge(
        uncertainty[[
            "profile_label", "yield_half_width_run",
            "R50_shift_q16_mm", "R50_shift_q50_mm", "R50_shift_q84_mm",
            "R50_shift_half_width_mm",
        ]],
        on="profile_label",
        validate="one_to_one",
    )
    result["yield_relative_half_width"] = (
        result.yield_half_width_run / result.expected_count_run)
    order = {label: position for position, label in enumerate(PROFILE_ORDER)}
    result["_order"] = result.profile_label.map(order)
    return result.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def _write_summary_tex(summary: pd.DataFrame, path: Path) -> None:
    labels = {
        "C11": "$^{11}$C",
        "O15": "$^{15}$O",
        "N13": "$^{13}$N",
        "all_production": "All production",
        "all_d120s300": "All in-room",
    }
    with path.open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrrr}\n\\toprule\n")
        stream.write(
            "Profile & $R_{50}$ (mm) & Yield half-width (\\%) & "
            "$q_{50}(b)$ (mm) & $u_R$ (mm) \\\\\n")
        stream.write("\\midrule\n")
        for row in summary.to_dict("records"):
            stream.write(
                f"{labels[row['profile_label']]} & {row['R50_prod_mm']:.2f} & "
                f"{100.0 * row['yield_relative_half_width']:.1f} & "
                f"{row['R50_shift_q50_mm']:.3f} & "
                f"{row['R50_shift_half_width_mm']:.3f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")


def generate(
    folding_directory: Path,
    output_directory: Path,
    *,
    convergence_csv: Path | None = None,
    generated_directory: Path | None = None,
    context_label: str | None = None,
    source_label: str | None = None,
) -> list[Path]:
    """Validate folding products and generate all available diagnostics."""

    folding_directory = Path(folding_directory).resolve()
    output_directory = Path(output_directory).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    paths, _, bands, summary, uncertainty = _load_products(folding_directory)
    products = [
        output_directory / "production_profiles.pdf",
        output_directory / "r50_shifts.pdf",
        output_directory / "yield_ratios.pdf",
    ]
    _profiles_figure(bands, summary, products[0], context_label)
    _distribution_figure(
        summary, products[1], quantity="R50", context_label=context_label)
    _distribution_figure(
        summary, products[2], quantity="yield", context_label=context_label)

    convergence_path = Path(convergence_csv).resolve() if convergence_csv else None
    if convergence_path is not None:
        convergence = pd.read_csv(convergence_path)
        convergence_figure = output_directory / "energy_convergence.pdf"
        _convergence_figure(convergence, convergence_figure, context_label)
        products.append(convergence_figure)

    generated_directory = (
        Path(generated_directory).resolve()
        if generated_directory is not None else output_directory)
    generated_directory.mkdir(parents=True, exist_ok=True)
    quantitative = _quantitative_summary(summary, uncertainty)
    quantitative.to_csv(
        generated_directory / "folding_plot_summary.csv", index=False,
        lineterminator="\n")
    _write_summary_tex(
        quantitative, generated_directory / "folding_plot_summary.tex")
    metadata = {
        "schema_version": 1,
        "folding_source": source_label or str(folding_directory),
        "inputs_sha256": {name: _sha256(path) for name, path in paths.items()},
        "convergence_source": (
            source_label if convergence_path is not None and source_label
            else str(convergence_path) if convergence_path else None),
        "convergence_sha256": (
            _sha256(convergence_path) if convergence_path is not None else None),
        "figures": [path.name for path in products],
        "context_label": context_label,
    }
    (generated_directory / "folding_plot_meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return products


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folding_directory", type=Path)
    parser.add_argument("--convergence-csv", type=Path)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("docs/figures/xsection_folding"))
    parser.add_argument(
        "--generated-dir", type=Path,
        default=Path("docs/generated/xsection_folding"))
    parser.add_argument("--context-label")
    parser.add_argument("--source-label")
    arguments = parser.parse_args()
    products = generate(
        arguments.folding_directory,
        arguments.output_dir,
        convergence_csv=arguments.convergence_csv,
        generated_directory=arguments.generated_dir,
        context_label=arguments.context_label,
        source_label=arguments.source_label,
    )
    print(f"generated {len(products)} folding diagnostic figures")


if __name__ == "__main__":
    main()
