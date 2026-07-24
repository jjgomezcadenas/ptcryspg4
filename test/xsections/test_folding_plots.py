import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_transport.xsections.exposure_folding import distal_r50
from analysis_transport.xsections.make_folding_plots import PROFILE_ORDER, generate


def write_synthetic_products(directory: Path) -> None:
    depth = np.arange(0.5, 8.0, 1.0)
    base = np.asarray([0.20, 0.55, 1.00, 0.92, 0.70, 0.43, 0.19, 0.05])
    scales = {
        "C11": 1.0,
        "O15": 1.3,
        "N13": 0.35,
        "all_production": 2.65,
        "all_inroom": 1.15,
    }
    nominal_rows = []
    band_rows = []
    summary_rows = []
    uncertainty_rows = []
    for label_position, label in enumerate(PROFILE_ORDER):
        nominal_profile = scales[label] * np.roll(base, label_position % 2)
        nominal_r50 = distal_r50(depth, nominal_profile)
        nominal_yield = float(nominal_profile.sum())
        quantity = (
            "measured_decays" if label == "all_inroom" else "production_nuclei")
        replicas = []
        replica_r50 = []
        for replica_id in range(30):
            tilt = (replica_id - 14.5) / 145.0
            profile = nominal_profile * (
                1.0 + tilt * np.linspace(-0.6, 1.0, len(depth)))
            replicas.append(profile)
            edge = distal_r50(depth, profile)
            replica_r50.append(edge)
            summary_rows.append({
                "model": "replica",
                "replica_id": replica_id,
                "profile_label": label,
                "expected_count_run": float(profile.sum()),
                "R50_prod_mm": edge,
                "R50_shift_mm": edge - nominal_r50,
            })
        replicas = np.asarray(replicas)
        replica_yields = replicas.sum(axis=1)
        shifts = np.asarray(replica_r50) - nominal_r50
        q_profile = np.quantile(replicas, [0.16, 0.50, 0.84], axis=0)
        for depth_index, depth_mm in enumerate(depth):
            nominal_rows.append({
                "profile_label": label,
                "quantity": quantity,
                "depth_mm": depth_mm,
                "expected_count_run": nominal_profile[depth_index],
            })
            band_rows.append({
                "profile_label": label,
                "quantity": quantity,
                "depth_mm": depth_mm,
                "nominal_run": nominal_profile[depth_index],
                "q16_run": q_profile[0, depth_index],
                "q50_run": q_profile[1, depth_index],
                "q84_run": q_profile[2, depth_index],
            })
        summary_rows.append({
            "model": "nominal",
            "replica_id": "",
            "profile_label": label,
            "expected_count_run": nominal_yield,
            "R50_prod_mm": nominal_r50,
            "R50_shift_mm": 0.0,
        })
        yield_q = np.quantile(replica_yields, [0.16, 0.50, 0.84])
        shift_q = np.quantile(shifts, [0.16, 0.50, 0.84])
        uncertainty_rows.append({
            "profile_label": label,
            "nominal_yield_run": nominal_yield,
            "yield_half_width_run": 0.5 * (yield_q[2] - yield_q[0]),
            "R50_shift_q16_mm": shift_q[0],
            "R50_shift_q50_mm": shift_q[1],
            "R50_shift_q84_mm": shift_q[2],
            "R50_shift_half_width_mm": 0.5 * (shift_q[2] - shift_q[0]),
        })
    pd.DataFrame(nominal_rows).to_csv(
        directory / "nominal_isotope_profiles.csv", index=False)
    pd.DataFrame(band_rows).to_csv(directory / "profile_bands.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(
        directory / "production_summary.csv", index=False)
    pd.DataFrame(uncertainty_rows).to_csv(
        directory / "uncertainty_summary.csv", index=False)


def write_synthetic_convergence(path: Path) -> None:
    rows = []
    for width in (0.5, 1.0, 2.0, 5.0):
        for label_position, label in enumerate(PROFILE_ORDER):
            scale = 1.0 + 0.05 * label_position
            rows.append({
                "candidate_width_MeV": width,
                "profile_label": label,
                "max_paired_yield_change_run": width * 0.001 * scale,
                "yield_replica_half_width_run": 0.1 * scale,
                "max_paired_R50_change_mm": width * 0.008 * scale,
                "R50_replica_half_width_mm": 0.5 * scale,
            })
    pd.DataFrame(rows).to_csv(path, index=False)


class FoldingPlotTests(unittest.TestCase):
    def test_generate_writes_four_figures_and_quantitative_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folding = root / "folding"
            figures = root / "figures"
            generated = root / "generated"
            folding.mkdir()
            write_synthetic_products(folding)
            convergence = root / "convergence.csv"
            write_synthetic_convergence(convergence)
            products = generate(
                folding,
                figures,
                convergence_csv=convergence,
                generated_directory=generated,
            )
            self.assertEqual(len(products), 4)
            self.assertTrue(all(path.stat().st_size > 1000 for path in products))
            summary = pd.read_csv(generated / "folding_plot_summary.csv")
            self.assertEqual(summary.profile_label.tolist(), list(PROFILE_ORDER))
            self.assertTrue((summary.yield_relative_half_width > 0).all())
            metadata = json.loads(
                (generated / "folding_plot_meta.json").read_text())
            self.assertEqual(len(metadata["inputs_sha256"]), 4)
            self.assertEqual(len(metadata["figures"]), 4)

    def test_incomplete_folding_directory_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "incomplete"):
                generate(Path(directory), Path(directory) / "figures")


if __name__ == "__main__":
    unittest.main()
