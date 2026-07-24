import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_transport.xsections.make_folding_reference import analytic_exposure
from analysis_transport.xsections.make_folding_plots import PROFILE_ORDER


REPO = Path(__file__).resolve().parents[2]


class FoldingReferenceTests(unittest.TestCase):
    def test_analytic_exposure_has_declared_grids_and_distal_zero_energy_bins(self):
        exposure, energy_edges, depth_edges = analytic_exposure()
        self.assertEqual(energy_edges[0], 5.0)
        self.assertEqual(energy_edges[-1], 120.0)
        self.assertEqual(depth_edges[0], 0.0)
        self.assertEqual(depth_edges[-1], 122.0)
        self.assertTrue(np.isfinite(exposure.target_exposure_cm2_inv).all())
        self.assertEqual(set(exposure.target), {"C12", "N14", "O16"})
        distal = exposure.loc[exposure.depth_low_mm >= 112.0]
        self.assertTrue((distal.energy_high_MeV == 5.5).all())

    def test_generated_reference_is_complete_and_identified_as_analytic(self):
        generated = REPO / "docs/generated/xsection_folding/reference"
        figures = REPO / "docs/figures/xsection_folding/reference"
        expected_generated = {
            "exposure_convergence.csv",
            "folding_plot_meta.json",
            "folding_plot_summary.csv",
            "folding_plot_summary.tex",
            "nominal_profiles.csv",
            "production_summary.csv",
            "profile_bands.csv",
            "reference_definition.json",
            "uncertainty_summary.csv",
        }
        self.assertTrue(expected_generated.issubset(
            {path.name for path in generated.iterdir()}))
        self.assertEqual(
            {path.name for path in figures.glob("*.pdf")},
            {
                "production_profiles.pdf", "r50_shifts.pdf",
                "yield_ratios.pdf", "energy_convergence.pdf",
            },
        )
        self.assertTrue(all(path.stat().st_size > 1000 for path in figures.glob("*.pdf")))
        definition = json.loads((generated / "reference_definition.json").read_text())
        self.assertEqual(definition["replicas"], 1000)
        self.assertIn("not a Geant4 treatment result", definition["parameters"]["description"])
        summary = pd.read_csv(generated / "folding_plot_summary.csv")
        self.assertEqual(summary.profile_label.tolist(), list(PROFILE_ORDER))
        metadata = json.loads((generated / "folding_plot_meta.json").read_text())
        self.assertEqual(metadata["folding_source"], "analytic_folding_reference_v1")


if __name__ == "__main__":
    unittest.main()
