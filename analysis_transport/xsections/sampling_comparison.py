"""Compare data-driven sampled production with the native Geant4 tally.

Reads one sampling run of stageA_xsection_ensemble (sampled_productions.csv,
native_route_counts.csv, exposure_meta.json, the folding products, and
emitters.csv when the emitter transport has run) and writes the paired
comparison of docs/sampling_xsections.tex: depth profiles, per-Gy yields,
the distal-edge displacements Delta = R50(data) - R50(G4), the
fold-versus-sampled closure, and the positron-range check.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from math import log
from common.isotopes import ISOTOPES as ISOTOPE_TABLE, NAME_TO_ID
from decay_sampling.scenarios import DEFAULT_SCENARIO_CONFIG, resolve_scenario

from .exposure_folding import distal_r50, load_exposure_metadata

ISOTOPES = ("O15", "C11", "N13")
COLORS = {"O15": "#0072B2", "C11": "#D55E00", "N13": "#009E73"}


def _histogram(depths, edges):
    counts, _ = np.histogram(depths, bins=edges)
    return counts.astype(float)


def _delta_r50_error(centres, data, g4, n_boot=400, seed=5):
    """Poisson-bootstrap standard deviation of R50(data) - R50(g4)."""
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(n_boot):
        a = distal_r50(centres, rng.poisson(data).astype(float))
        b = distal_r50(centres, rng.poisson(g4).astype(float))
        if np.isfinite(a) and np.isfinite(b):
            deltas.append(a - b)
    return float(np.std(deltas)) if deltas else float("nan")


def analyze(run_dir: Path, repo: Path):
    metadata = load_exposure_metadata(run_dir / "exposure_meta.json",
                                      run_dir / "proton_exposure.csv")
    half_extent = metadata.depth_edges_mm[-1] / 2.0
    edges = metadata.depth_edges_mm
    centres = 0.5 * (edges[:-1] + edges[1:])

    sampled = pd.read_csv(run_dir / "sampled_productions.csv")
    sampled["depth_mm"] = sampled.prod_z_mm + half_extent
    native = pd.read_csv(run_dir / "native_route_counts.csv")
    native["depth_mm"] = 0.5 * (native.depth_low_mm + native.depth_high_mm)
    folded = pd.read_csv(run_dir / "folding/nominal_isotope_profiles.csv")

    rows = []
    profiles = {}
    for isotope in ISOTOPES:
        data = _histogram(sampled.loc[sampled.residual == isotope, "depth_mm"],
                          edges)
        g4_rows = native[native.residual == isotope]
        g4 = np.zeros(len(centres))
        for _, row in g4_rows.iterrows():
            index = int(row.depth_low_mm / (edges[1] - edges[0]))
            if 0 <= index < len(g4):
                g4[index] += row.production_count
        fold_rows = folded[(folded.profile_label == isotope)
                           & (folded.quantity == "production_nuclei")]
        fold = np.interp(centres, fold_rows.depth_mm,
                         fold_rows.expected_count_run, left=0., right=0.)
        profiles[isotope] = (data, g4, fold)
        r50_data = distal_r50(centres, data)
        r50_g4 = distal_r50(centres, g4)
        rows.append({
            "isotope": isotope,
            "sampled_count": int(data.sum()),
            "native_count": int(g4.sum()),
            "folded_expected": float(fold.sum()),
            "sampled_per_Gy": data.sum() / metadata.target_dose_Gy,
            "native_per_Gy": g4.sum() / metadata.target_dose_Gy,
            "yield_ratio_data_over_g4": (data.sum() / g4.sum())
            if g4.sum() else np.inf,
            "sampled_over_folded": data.sum() / fold.sum(),
            "R50_data_mm": r50_data,
            "R50_g4_mm": r50_g4,
            "delta_R50_mm": r50_data - r50_g4,
            "delta_R50_err_mm": _delta_r50_error(centres, data, g4),
            "R50_fold_mm": distal_r50(centres, fold),
        })

    combined_data = sum(profiles[i][0] for i in ISOTOPES)
    combined_g4 = sum(profiles[i][1] for i in ISOTOPES)
    combined_fold = sum(profiles[i][2] for i in ISOTOPES)
    rows.append({
        "isotope": "combined",
        "sampled_count": int(combined_data.sum()),
        "native_count": int(combined_g4.sum()),
        "folded_expected": float(combined_fold.sum()),
        "sampled_per_Gy": combined_data.sum() / metadata.target_dose_Gy,
        "native_per_Gy": combined_g4.sum() / metadata.target_dose_Gy,
        "yield_ratio_data_over_g4": combined_data.sum() / combined_g4.sum(),
        "sampled_over_folded": combined_data.sum() / combined_fold.sum(),
        "R50_data_mm": distal_r50(centres, combined_data),
        "R50_g4_mm": distal_r50(centres, combined_g4),
        "delta_R50_mm": distal_r50(centres, combined_data)
        - distal_r50(centres, combined_g4),
        "delta_R50_err_mm": _delta_r50_error(centres, combined_data,
                                             combined_g4),
        "R50_fold_mm": distal_r50(centres, combined_fold),
    })
    # Scenario-weighted profiles: the handoff factors of each named
    # acquisition scenario applied to both tallies. The compensation seen in
    # raw production depends on the isotope weights; these rows measure how
    # the paired displacement moves with the acquisition timing.
    for scenario_name in ("inroom", "fast", "offline"):
        scenario = resolve_scenario(scenario_name, DEFAULT_SCENARIO_CONFIG)
        factors = {
            isotope: scenario.measured_fraction(
                log(2.0) / ISOTOPE_TABLE[NAME_TO_ID[isotope]].half_life_s)
            for isotope in ISOTOPES
        }
        weighted_data = sum(factors[i] * profiles[i][0] for i in ISOTOPES)
        weighted_g4 = sum(factors[i] * profiles[i][1] for i in ISOTOPES)
        rng = np.random.default_rng(9)
        boot = []
        for _ in range(400):
            a = distal_r50(centres,
                           sum(factors[i] * rng.poisson(profiles[i][0])
                               for i in ISOTOPES))
            b = distal_r50(centres,
                           sum(factors[i] * rng.poisson(profiles[i][1])
                               for i in ISOTOPES))
            if np.isfinite(a) and np.isfinite(b):
                boot.append(a - b)
        rows.append({
            "isotope": f"combined_{scenario_name}",
            "sampled_count": float(weighted_data.sum()),
            "native_count": float(weighted_g4.sum()),
            "folded_expected": float("nan"),
            "sampled_per_Gy": weighted_data.sum() / metadata.target_dose_Gy,
            "native_per_Gy": weighted_g4.sum() / metadata.target_dose_Gy,
            "yield_ratio_data_over_g4":
                weighted_data.sum() / weighted_g4.sum(),
            "sampled_over_folded": float("nan"),
            "R50_data_mm": distal_r50(centres, weighted_data),
            "R50_g4_mm": distal_r50(centres, weighted_g4),
            "delta_R50_mm": distal_r50(centres, weighted_data)
            - distal_r50(centres, weighted_g4),
            "delta_R50_err_mm": float(np.std(boot)) if boot else float("nan"),
            "R50_fold_mm": float("nan"),
        })
    summary = pd.DataFrame(rows)

    # Positron-range check against the stageA native reference.
    ranges = []
    emitters_path = run_dir / "emitters.csv"
    reference_path = repo / "data/runs/uniform_headep_sobp_1e7/emitters.csv"
    if emitters_path.exists() and reference_path.exists():
        ours = pd.read_csv(emitters_path)
        reference = pd.read_csv(reference_path)
        for iid, isotope in enumerate(ISOTOPES):
            a = ours[ours.isotope_id == iid]
            b = reference[reference.isotope_id == iid]
            ra = np.sqrt((a.anh_x_mm - a.prod_x_mm) ** 2
                         + (a.anh_y_mm - a.prod_y_mm) ** 2
                         + (a.anh_z_mm - a.prod_z_mm) ** 2)
            rb = np.sqrt((b.anh_x_mm - b.prod_x_mm) ** 2
                         + (b.anh_y_mm - b.prod_y_mm) ** 2
                         + (b.anh_z_mm - b.prod_z_mm) ** 2)
            ranges.append({
                "isotope": isotope, "n_sampled": len(a), "n_reference": len(b),
                "median_range_mm": float(np.median(ra)),
                "median_reference_mm": float(np.median(rb)),
                "q90_range_mm": float(np.quantile(ra, 0.9)),
                "q90_reference_mm": float(np.quantile(rb, 0.9)),
            })

    # Figure: three isotope panels plus the combined profile.
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 6.8), sharex=True)
    panels = list(ISOTOPES) + ["combined"]
    for axis, label in zip(axes.flat, panels):
        if label == "combined":
            data, g4, fold = combined_data, combined_g4, combined_fold
            color = "black"
        else:
            data, g4, fold = profiles[label]
            color = COLORS[label]
        axis.step(centres, g4, where="mid", color="0.45", linewidth=1.2,
                  label="native Geant4")
        axis.plot(centres, fold, color=color, linewidth=1.1, alpha=0.7,
                  label="folded nominal")
        axis.step(centres, data, where="mid", color=color, linewidth=1.4,
                  linestyle="--", label="sampled from fits")
        axis.set_title(label, fontsize=10)
        axis.set_xlim(0, 130)
        axis.set_ylim(bottom=0)
        axis.grid(alpha=0.18)
        axis.legend(fontsize=7.5)
        axis.set_ylabel("productions / bin", fontsize=9)
    for axis in axes[1]:
        axis.set_xlabel("Depth (mm)", fontsize=9)
    figure.suptitle("Production profiles: fitted cross sections vs native model",
                    fontsize=12)
    figure.tight_layout()
    figure_dir = repo / "docs/figures/sampling_xsections"
    figure_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_dir / "production_comparison.pdf")
    plt.close(figure)

    generated = repo / "docs/generated/sampling_xsections"
    generated.mkdir(parents=True, exist_ok=True)
    summary.to_csv(generated / "sampling_summary.csv", index=False,
                   lineterminator="\n")
    if ranges:
        pd.DataFrame(ranges).to_csv(generated / "positron_ranges.csv",
                                    index=False, lineterminator="\n")
    # Cross-section systematic of the data-driven edge, from the replica fold.
    uncertainty = pd.read_csv(run_dir / "folding/uncertainty_summary.csv")
    band_labels = {"O15": "$^{15}$O", "C11": "$^{11}$C", "N13": "$^{13}$N",
                   "all_production": "combined (production)",
                   "all_inroom": "combined (in-room weighted)"}
    with (generated / "systematic_band.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrr}\n\\toprule\n")
        stream.write("Profile & $R_{50}$ (mm) & "
                     "$u_{\\mathrm{xs}}$ (mm) \\\\\n\\midrule\n")
        for key in ("O15", "C11", "N13", "all_production", "all_inroom"):
            row = uncertainty[uncertainty.profile_label == key].iloc[0]
            stream.write(
                f"{band_labels[key]} & {row.R50_nominal_mm:.2f} & "
                f"$\\pm${row.R50_shift_half_width_mm:.2f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")

    with (generated / "delta_r50.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrrrr}\n\\toprule\n")
        stream.write("Isotope & $N^{\\mathrm{data}}$/Gy & $N^{\\mathrm{G4}}$/Gy"
                     " & ratio & $R_{50}^{\\mathrm{data}}$ (mm) & "
                     "$\\Delta$ (mm) \\\\\n\\midrule\n")
        for row in summary.to_dict("records"):
            label = row["isotope"].replace("_", " ")
            stream.write(
                f"{label} & {row['sampled_per_Gy']:.3e} & "
                f"{row['native_per_Gy']:.3e} & "
                f"{row['yield_ratio_data_over_g4']:.3f} & "
                f"{row['R50_data_mm']:.2f} & "
                f"{row['delta_R50_mm']:+.2f} $\\pm$ "
                f"{row['delta_R50_err_mm']:.2f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")
    with (generated / "sampling_meta.json").open("w", encoding="utf-8") as stream:
        json.dump({"run_id": metadata.run_id,
                   "n_protons": metadata.n_protons,
                   "target_dose_Gy": metadata.target_dose_Gy,
                   "sampling_meta": json.loads(
                       (run_dir / "sampling_meta.json").read_text())},
                  stream, indent=2)
        stream.write("\n")
    return summary, ranges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    summary, ranges = analyze(args.run_dir.resolve(), args.repo.resolve())
    print(summary[["isotope", "yield_ratio_data_over_g4", "sampled_over_folded",
                   "R50_data_mm", "R50_g4_mm", "delta_R50_mm"]]
          .to_string(index=False))
    for row in ranges:
        print(f"positron range {row['isotope']}: median "
              f"{row['median_range_mm']:.3f} vs reference "
              f"{row['median_reference_mm']:.3f} mm")


if __name__ == "__main__":
    main()
