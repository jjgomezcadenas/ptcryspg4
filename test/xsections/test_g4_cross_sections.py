import unittest

from analysis_transport.xsections.g4_cross_sections import (
    ZERO_COUNT_95_UPPER, estimate, estimate_factorized,
)
from analysis_transport.xsections.run_g4_denominator import active_residuals, merge_rows as merge_denominator_rows
from analysis_transport.xsections.run_g4_scan import merge_rows, select_thickness


class G4CrossSectionTests(unittest.TestCase):
    def test_known_estimator(self):
        sigma, uncertainty, upper = estimate(100, 1_000_000, 2.0e20)
        self.assertAlmostEqual(sigma, 500.0)
        self.assertAlmostEqual(uncertainty, 50.0)
        self.assertEqual(upper, "")

    def test_zero_count_limit(self):
        sigma, uncertainty, upper = estimate(0, 1_000_000, 2.0e20)
        self.assertEqual(sigma, 0.0)
        self.assertEqual(uncertainty, "")
        self.assertAlmostEqual(upper, ZERO_COUNT_95_UPPER * 5.0)

    def test_factorized_estimator(self):
        sigma, uncertainty, upper = estimate_factorized(100, 1000, 400.0)
        self.assertAlmostEqual(sigma, 40.0)
        self.assertAlmostEqual(uncertainty, 3.794733192, places=8)
        self.assertEqual(upper, "")

    def test_factorized_zero_count_limit(self):
        sigma, uncertainty, upper = estimate_factorized(0, 1000, 400.0)
        self.assertEqual(sigma, 0.0)
        self.assertEqual(uncertainty, "")
        self.assertAlmostEqual(upper, ZERO_COUNT_95_UPPER * 0.4)

    def test_denominator_batch_merge(self):
        base = {
            "n_interactions": "100", "n_c11": "10", "n_n13": "2",
            "n_o15": "5", "n_secondaries": "300", "n_nuclei": "100",
            "sigma_inelastic_mb": "400", "target_z": "8", "target_a": "16",
            "cross_section_data_sets": "xs", "physics_list": "BIC",
            "geant4_version": "11.4.1", "seed": "1",
        }
        other = dict(base, n_c11="20", seed="2")
        merged = merge_denominator_rows([base, other])
        self.assertEqual(merged["n_interactions"], 200)
        self.assertEqual(merged["n_c11"], 30)
        self.assertEqual(merged["seeds"], "1;2")

    def test_active_residual_thresholds(self):
        config = {"target": {"O16": {
            "active_residuals": ["C11", "N13", "O15"],
            "threshold_MeV": {"C11": 25.17, "N13": 5.5, "O15": 14.29},
        }}}
        self.assertEqual(active_residuals(config, "O16", 10.0), ["N13"])

    def test_thickness_selection(self):
        config = {"maximum_mean_energy_loss_MeV": 0.5,
                  "maximum_inelastic_probability": 1.0e-3}
        rows = [
            {"n_protons": "20000", "n_inelastic": "30", "mean_continuous_loss_MeV": "0.3",
             "areal_mg_cm2": "80"},
            {"n_protons": "20000", "n_inelastic": "5", "mean_continuous_loss_MeV": "0.2",
             "areal_mg_cm2": "40"},
        ]
        self.assertEqual(select_thickness(rows, config), 40.0)


if __name__ == "__main__":
    unittest.main()
