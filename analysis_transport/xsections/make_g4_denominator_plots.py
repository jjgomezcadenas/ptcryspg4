"""Plot the direct Geant4 denominator pilot against the experimental fits."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .channels import CHANNELS


def generate(repo: Path):
    pilot_dir = (repo / "data/xsections/g4/geant4-11.4.1/QGSP_BIC_HP"
                 / "denominator_pilot")
    summary = pd.read_csv(pilot_dir / "support_summary.csv").fillna("")
    figure, axes = plt.subplots(3, 2, figsize=(10.2, 11.0), sharex=True)
    for axis, channel in zip(axes.flat, CHANNELS):
        experimental = pd.read_csv(
            repo / "data/xsections/fits" / f"{channel.channel_id}_curve.csv")
        g4 = pd.read_csv(pilot_dir / f"{channel.channel_id}_curve.csv")
        axis.fill_between(
            experimental.energy_MeV, experimental.sigma_lower_16_mb,
            experimental.sigma_upper_84_mb, color="#56B4E9", alpha=0.3,
            label="EXFOR 16--84%")
        axis.plot(experimental.energy_MeV, experimental.sigma_nominal_mb,
                  color="black", linewidth=1.5, label="EXFOR nominal")
        positive = g4.residual_count > 0
        axis.errorbar(
            g4.loc[positive, "energy_MeV"], g4.loc[positive, "sigma_mb"],
            yerr=g4.loc[positive, "sigma_unc_stat_mb"], linestyle="none",
            marker="o", markersize=3.0, color="#D55E00", elinewidth=0.7,
            label="G4 direct pilot")
        zero = ~positive
        if zero.any():
            limits = g4.loc[zero, "sigma_upper_95_mb"].to_numpy()
            axis.errorbar(
                g4.loc[zero, "energy_MeV"], limits,
                yerr=np.maximum(0.25 * limits, 0.05), uplims=True,
                linestyle="none", marker="v", markersize=3.0,
                color="#CC79A7", elinewidth=0.7, label="G4 95% upper limit")
        axis.set_title(channel.title)
        axis.set_xlim(5, 150)
        axis.set_ylim(bottom=0)
        axis.set_ylabel("Cross section (mb)")
        axis.grid(alpha=0.18)
    axes.flat[-1].axis("off")
    axes.flat[0].legend(fontsize=7.5)
    for axis in axes[-1, :]:
        if axis.get_visible():
            axis.set_xlabel("Proton energy (MeV)")
    axes[1, 0].set_xlabel("Proton energy (MeV)")
    axes[1, 1].set_xlabel("Proton energy (MeV)")
    figure.suptitle("QGSP_BIC_HP direct final-state denominator pilot")
    figure.tight_layout()
    output = repo / "docs/figures/xsection_denominator/pilot_support.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)

    generated = repo / "docs/generated/xsection_denominator"
    generated.mkdir(parents=True, exist_ok=True)
    summary.to_csv(generated / "support_summary.csv", index=False,
                   lineterminator="\n")
    with (generated / "support_summary.tex").open(
            "w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrl}\n\\toprule\n")
        stream.write("Channel & Zero-support points & First positive (MeV) & "
                     "Zero-support energies (MeV) \\\\\n")
        stream.write("\\midrule\n")
        for row in summary.to_dict("records"):
            title = next(channel.title for channel in CHANNELS
                         if channel.channel_id == row["channel_id"])
            energies = str(row["zero_support_energies_MeV"]).replace(";", ", ")
            stream.write(
                f"{title} & {row['zero_g4_support_in_production_region']} & "
                f"{row['first_positive_g4_energy_MeV']} & {energies or '--'} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    output = generate(args.repo.resolve())
    print(f"Wrote Geant4 denominator pilot figure to {output}")


if __name__ == "__main__":
    main()
