"""Validate a source bank: unbiasedness closure and effective sample size.

Checks, per docs/sampling_xsections.tex (Sec. source bank):
1. Closure — for the nominal curve and every replica, the bank estimate of
   each channel's total production (sum of n*l*sigma/q) must agree with the
   folded expectation within the bank's own sampling error.
2. ESS floors — the effective sample size (sum w)^2 / sum w^2 per channel
   and replica, total and in the distal window, must exceed the configured
   floors.

Writes bank_validation.json into the run directory and exits non-zero on
any violation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .exposure_folding import CrossSectionEnsemble, load_exposure_metadata

ESS_FLOOR_TOTAL = 5.0e3       # per channel, worst replica
ESS_FLOOR_DISTAL = 5.0e2      # per channel, worst replica, distal window
DISTAL_WINDOW_MM = 15.0       # window below the channel's 99th-percentile
                              # production depth (nominal weights)
CLOSURE_SIGMAS = 4.0          # allowed pulls for the closure test


def validate(run_dir: Path, fit_dir: Path):
    bank = pd.read_csv(run_dir / "source_bank.csv")
    metadata = load_exposure_metadata(run_dir / "exposure_meta.json",
                                      run_dir / "proton_exposure.csv")
    half_extent = metadata.depth_edges_mm[-1] / 2.0
    bank["depth_mm"] = bank.z_mm + half_extent
    ensemble = CrossSectionEnsemble.from_fit_directory(fit_dir)

    # Folded per-replica channel totals: recompute directly from the exposure
    # table so the closure reference is independent of the bank.
    exposure = pd.read_csv(run_dir / "proton_exposure.csv")

    report = {"bank_entries": int(len(bank)), "channels": {}}
    failures = []
    for channel_id, entries in bank.groupby("channel_id"):
        rows = exposure[exposure.target == entries.target.iloc[0]]
        energies = entries.proton_energy_MeV.to_numpy()
        base = entries.exposure_cm2_inv.to_numpy() / \
            entries.keep_probability.to_numpy()

        nominal_sigma = ensemble.nominal(channel_id, energies)
        replica_sigma = ensemble.replicas(channel_id, energies)

        fold_nominal = float((rows.target_exposure_cm2_inv * np.asarray(
            ensemble.nominal(channel_id, rows.energy_mean_MeV.to_numpy()))
        ).sum() * 1.0e-27)
        fold_replicas = (rows.target_exposure_cm2_inv.to_numpy()[None, :]
                         * np.asarray(ensemble.replicas(
                             channel_id, rows.energy_mean_MeV.to_numpy()))
                         ).sum(axis=1) * 1.0e-27

        weights = replica_sigma * base[None, :] * 1.0e-27
        nominal_weights = nominal_sigma * base * 1.0e-27
        bank_totals = weights.sum(axis=1)
        bank_nominal = float(nominal_weights.sum())

        # Bank sampling error of a total: sqrt(sum w^2 (1-q)) per replica.
        one_minus_q = 1.0 - entries.keep_probability.to_numpy()
        error_nominal = float(np.sqrt(
            (nominal_weights ** 2 * one_minus_q).sum()))
        errors = np.sqrt((weights ** 2 * one_minus_q[None, :]).sum(axis=1))
        pulls = (bank_totals - fold_replicas) / np.maximum(errors, 1e-300)
        pull_nominal = (bank_nominal - fold_nominal) / max(error_nominal,
                                                           1e-300)

        ess_total = (weights.sum(axis=1) ** 2 /
                     np.maximum((weights ** 2).sum(axis=1), 1e-300))
        # Distal window anchored to the production edge: the 99th percentile
        # of the nominal-weighted depth, not the deepest stray entry.
        order = np.argsort(entries.depth_mm.to_numpy())
        cumulative = np.cumsum(nominal_weights[order])
        d99 = float(entries.depth_mm.to_numpy()[order][
            np.searchsorted(cumulative, 0.99 * cumulative[-1])])
        distal = entries.depth_mm > d99 - DISTAL_WINDOW_MM
        wd = weights[:, distal.to_numpy()]
        ess_distal = (wd.sum(axis=1) ** 2 /
                      np.maximum((wd ** 2).sum(axis=1), 1e-300))

        summary = {
            "entries": int(len(entries)),
            "bank_over_fold_nominal": bank_nominal / fold_nominal,
            "closure_pull_nominal": float(pull_nominal),
            "worst_closure_pull_replicas": float(np.max(np.abs(pulls))),
            "min_ess_total": float(ess_total.min()),
            "min_ess_distal": float(ess_distal.min()),
        }
        report["channels"][channel_id] = summary
        if abs(pull_nominal) > CLOSURE_SIGMAS:
            failures.append(f"{channel_id}: nominal closure pull "
                            f"{pull_nominal:.2f}")
        if np.max(np.abs(pulls)) > CLOSURE_SIGMAS:
            failures.append(f"{channel_id}: replica closure pull "
                            f"{np.max(np.abs(pulls)):.2f}")
        if ess_total.min() < ESS_FLOOR_TOTAL:
            failures.append(f"{channel_id}: ESS total {ess_total.min():.0f}"
                            f" < {ESS_FLOOR_TOTAL:.0f}")
        if ess_distal.min() < ESS_FLOOR_DISTAL:
            failures.append(f"{channel_id}: ESS distal {ess_distal.min():.0f}"
                            f" < {ESS_FLOOR_DISTAL:.0f}")

    report["ess_floor_total"] = ESS_FLOOR_TOTAL
    report["ess_floor_distal"] = ESS_FLOOR_DISTAL
    report["distal_window_mm"] = DISTAL_WINDOW_MM
    report["closure_sigmas"] = CLOSURE_SIGMAS
    report["failures"] = failures
    (run_dir / "bank_validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--fit-dir", type=Path,
                        default=Path("data/xsections/fits"))
    args = parser.parse_args()
    report = validate(args.run_dir.resolve(), args.fit_dir.resolve())
    for channel, summary in report["channels"].items():
        print(f"{channel}: entries {summary['entries']}, "
              f"bank/fold {summary['bank_over_fold_nominal']:.4f}, "
              f"worst pull {summary['worst_closure_pull_replicas']:.2f}, "
              f"ESS total {summary['min_ess_total']:.0f}, "
              f"distal {summary['min_ess_distal']:.0f}")
    if report["failures"]:
        raise SystemExit("; ".join(report["failures"]))
    print("bank validation passed")


if __name__ == "__main__":
    main()
