"""Quantify the hadronic-transport term of the distal-edge error budget.

Production is exposure times fitted cross section, so the hadronic model
enters only through the exposure map: the secondary protons it adds and the
nonelastic attenuation of the primary fluence. This script folds the same
nominal fitted curves with exposure maps transported under alternate
final-state models (QGSP_BERT_HP, QGSP_INCLXX_HP versus the QGSP_BIC_HP
reference) and, separately, tilts the reference profiles by a +-3 percent
variation of the nonelastic attenuation. The R50 displacements of the folded
profiles are the transport term u_transport, quoted next to the
cross-section term u_xs of the same profiles.

Writes docs/generated/sampling_xsections/transport_envelope.{csv,tex}.
"""

from __future__ import annotations

import argparse
from math import log
from pathlib import Path

import numpy as np
import pandas as pd

from common.isotopes import ISOTOPES as ISOTOPE_TABLE, NAME_TO_ID
from common.phantom_material import MATERIALS, _N_A
from decay_sampling.scenarios import DEFAULT_SCENARIO_CONFIG, resolve_scenario

from .channels import CHANNELS
from .exposure_folding import (
    MB_TO_CM2,
    CrossSectionEnsemble,
    distal_r50,
    load_exposure_metadata,
    load_exposure_table,
)

ISOTOPES = ("O15", "C11", "N13")
SCENARIOS = ("d120s300", "d120s120", "d180s300", "d180s120", "d300s300")
REFERENCE_LIST = "QGSP_BIC_HP"
ALTERNATE_LISTS = ("QGSP_BERT_HP", "QGSP_INCLXX_HP")
ATTENUATION_FRACTION = 0.03  # relative variation of the nonelastic cross section


def brain_macroscopic_inelastic_cm() -> float:
    """Macroscopic proton nonelastic cross section of the brain phantom (1/cm).

    Per-element sigma_inel from the geometric estimate 45*A^(2/3) mb, which
    reproduces the measured ~100 MeV plateau values within ~15 percent —
    ample for a +-3 percent sensitivity of an optical depth below 0.1.
    Hydrogen has no proton nonelastic channel at therapy energies.
    """
    material = MATERIALS["G4_BRAIN_ICRP"]
    total = 0.0
    for symbol, element in material.composition.items():
        if symbol == "H":
            continue
        sigma_cm2 = 45.0 * element.A ** (2.0 / 3.0) * MB_TO_CM2
        number_density = material.density_g_cm3 * _N_A * element.mass_fraction / element.A
        total += number_density * sigma_cm2
    return total


def fold_nominal_profiles(run_dir: Path, ensemble: CrossSectionEnsemble):
    """Per-isotope nominal production-depth profiles of one exposure run."""
    metadata = load_exposure_metadata(run_dir / "exposure_meta.json",
                                      run_dir / "proton_exposure.csv")
    table = load_exposure_table(run_dir / "proton_exposure.csv")
    edges = metadata.depth_edges_mm
    centres = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    profiles = {isotope: np.zeros(len(centres)) for isotope in ISOTOPES}
    for channel in CHANNELS:
        rows = table.loc[table["target"] == channel.target]
        if rows.empty:
            continue
        sigma = ensemble.nominal(channel.channel_id,
                                 rows["energy_mean_MeV"].to_numpy(float))
        yields = rows["target_exposure_cm2_inv"].to_numpy(float) * sigma * MB_TO_CM2
        indices = np.rint(rows["depth_low_mm"].to_numpy(float) / width).astype(int)
        np.add.at(profiles[channel.residual], indices, yields)
    for isotope in ISOTOPES:
        profiles[isotope] /= metadata.n_protons
    return centres, profiles, metadata


def profile_set(profiles: dict) -> dict:
    """The reported profiles: isotopes, raw sum, and the CBS-weighted sums."""
    out = {isotope: profiles[isotope] for isotope in ISOTOPES}
    out["combined"] = sum(profiles[i] for i in ISOTOPES)
    for name in SCENARIOS:
        scenario = resolve_scenario(name, DEFAULT_SCENARIO_CONFIG)
        factors = {
            isotope: scenario.measured_fraction(
                log(2.0) / ISOTOPE_TABLE[NAME_TO_ID[isotope]].half_life_s)
            for isotope in ISOTOPES
        }
        out[f"combined_{name}"] = sum(factors[i] * profiles[i] for i in ISOTOPES)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exposure-base", type=Path,
                        default=Path("data/xsections/exposure"))
    parser.add_argument("--fit-dir", type=Path, default=Path("data/xsections/fits"))
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--protons-tag", default="1e7")
    parser.add_argument("--stat-tag", default="1e8",
                        help="reference-list run at higher statistics; its R50 "
                             "difference from the working run is the "
                             "statistical scale of the comparison")
    args = parser.parse_args()

    ensemble = CrossSectionEnsemble.from_fit_directory(args.fit_dir)

    def run_dir(physics_list: str, tag: str) -> Path:
        return args.exposure_base / f"uniform_headep_sobp_{physics_list}_{tag}"

    centres, reference_raw, reference_meta = fold_nominal_profiles(
        run_dir(REFERENCE_LIST, args.protons_tag), ensemble)
    reference = profile_set(reference_raw)

    alternates = {}
    for physics_list in ALTERNATE_LISTS:
        c, raw, meta = fold_nominal_profiles(run_dir(physics_list, args.protons_tag),
                                             ensemble)
        if not np.array_equal(c, centres):
            raise ValueError(f"{physics_list}: depth grid differs from reference")
        alternates[physics_list] = (profile_set(raw), meta)

    stat_c, stat_raw, _ = fold_nominal_profiles(
        run_dir(REFERENCE_LIST, args.stat_tag), ensemble)
    if not np.array_equal(stat_c, centres):
        raise ValueError("statistics run: depth grid differs from reference")
    stat = profile_set(stat_raw)

    # Attenuation tilt: the primary fluence at depth z carries the nonelastic
    # optical depth Sigma*z; a relative variation delta of Sigma multiplies
    # every production profile by exp(-+delta*Sigma*z).
    sigma_macroscopic = brain_macroscopic_inelastic_cm()
    z_cm = centres / 10.0
    tilt_low = np.exp(-ATTENUATION_FRACTION * sigma_macroscopic * z_cm)
    tilt_high = np.exp(+ATTENUATION_FRACTION * sigma_macroscopic * z_cm)

    rows = []
    for label, profile in reference.items():
        r50 = distal_r50(centres, profile)
        row = {
            "profile": label,
            "R50_reference_mm": r50,
            "yield_reference_per_proton": float(profile.sum()),
        }
        deltas = []
        for physics_list in ALTERNATE_LISTS:
            alt_profiles, _ = alternates[physics_list]
            delta = distal_r50(centres, alt_profiles[label]) - r50
            row[f"delta_R50_{physics_list}_mm"] = delta
            row[f"yield_ratio_{physics_list}"] = (
                alt_profiles[label].sum() / profile.sum())
            deltas.append(abs(delta))
        delta_atten = max(
            abs(distal_r50(centres, profile * tilt_low) - r50),
            abs(distal_r50(centres, profile * tilt_high) - r50))
        row["delta_R50_attenuation_mm"] = delta_atten
        deltas.append(delta_atten)
        row["u_transport_mm"] = max(deltas)
        row["delta_R50_statistics_mm"] = distal_r50(centres, stat[label]) - r50
        rows.append(row)
    summary = pd.DataFrame(rows)

    generated = args.repo / "docs/generated/sampling_xsections"
    generated.mkdir(parents=True, exist_ok=True)
    summary.to_csv(generated / "transport_envelope.csv", index=False,
                   lineterminator="\n")

    labels = {"O15": "$^{15}$O", "C11": "$^{11}$C", "N13": "$^{13}$N",
              "combined": "combined (production)",
              "combined_d120s300": "combined (reference scan)"}
    with (generated / "transport_envelope.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{tabular}{lrrrrr}\n\\toprule\n")
        stream.write(
            "Profile & $R_{50}$ (mm) & $\\Delta_{\\mathrm{BERT}}$ & "
            "$\\Delta_{\\mathrm{INCL}}$ & $\\Delta_{\\Sigma\\pm3\\%}$ & "
            "$u_{\\mathrm{transport}}$ (mm) \\\\\n\\midrule\n")
        for key, label in labels.items():
            row = summary[summary.profile == key].iloc[0]
            stream.write(
                f"{label} & {row.R50_reference_mm:.2f} & "
                f"{row.delta_R50_QGSP_BERT_HP_mm:+.2f} & "
                f"{row.delta_R50_QGSP_INCLXX_HP_mm:+.2f} & "
                f"$\\pm${row.delta_R50_attenuation_mm:.2f} & "
                f"{row.u_transport_mm:.2f} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")

    print(f"reference: {reference_meta.run_id}")
    print(f"brain macroscopic nonelastic cross section: "
          f"{sigma_macroscopic:.4f} /cm "
          f"(optical depth at 86 mm: {sigma_macroscopic * 8.6:.3f})")
    columns = ["profile", "R50_reference_mm", "delta_R50_QGSP_BERT_HP_mm",
               "delta_R50_QGSP_INCLXX_HP_mm", "delta_R50_attenuation_mm",
               "u_transport_mm", "delta_R50_statistics_mm"]
    print(summary[columns].to_string(index=False))


if __name__ == "__main__":
    main()
