import unittest

import pandas as pd

from analysis_transport.xsections.native_routes import classify_native_routes


class NativeRouteTests(unittest.TestCase):
    def test_selected_proton_routes_are_modeled_and_others_are_not(self):
        routes = pd.DataFrame({
            "projectile": ["proton", "neutron", "proton", "proton"],
            "target": ["C12", "C12", "O16", "O16"],
            "residual": ["C11", "C11", "C10", "O15"],
            "depth_low_mm": [0.0, 0.0, 0.0, 0.0],
            "depth_high_mm": [1.0, 1.0, 1.0, 1.0],
            "production_count": [80.0, 5.0, 10.0, 5.0],
        })
        detailed, summary = classify_native_routes(routes)
        self.assertEqual(
            detailed["represented_by_fold"].tolist(),
            [True, False, False, True],
        )
        production = summary.loc[
            summary["profile_label"] == "all_production"].iloc[0]
        self.assertAlmostEqual(production["unmodeled_fraction"], 0.15)
        inroom = summary.loc[summary["profile_label"] == "all_inroom"].iloc[0]
        self.assertLess(inroom["unmodeled_fraction"], production["unmodeled_fraction"])

    def test_negative_counts_are_rejected(self):
        routes = pd.DataFrame({
            "projectile": ["proton"],
            "target": ["C12"],
            "residual": ["C11"],
            "depth_low_mm": [0.0],
            "depth_high_mm": [1.0],
            "production_count": [-1.0],
        })
        with self.assertRaisesRegex(ValueError, "non-negative"):
            classify_native_routes(routes)


if __name__ == "__main__":
    unittest.main()
