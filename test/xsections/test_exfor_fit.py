import unittest

import numpy as np
import pandas as pd

from analysis_transport.xsections.fit_exfor import (
    PSpline, fit_coefficients, inverse_covariance,
)


class ExforFitTests(unittest.TestCase):
    def setUp(self):
        self.model = PSpline.create(
            threshold=5.0, energy_max=40.0, threshold_power=0.5,
            internal_knots=[8.0, 12.0, 18.0, 25.0, 32.0])

    def test_threshold_and_positivity(self):
        coefficients = np.zeros(self.model.coefficient_count)
        prediction = self.model.predict(coefficients, [4.0, 5.0, 6.0, 20.0])
        np.testing.assert_array_equal(prediction[:2], [0.0, 0.0])
        self.assertTrue((prediction[2:] > 0).all())

    def test_campaign_spread_creates_off_diagonal_covariance(self):
        frame = pd.DataFrame({"campaign_id": ["A", "A", "B"]})
        inverse = inverse_covariance(frame, np.full(3, 0.1), 0.2)
        self.assertLess(inverse[0, 1], 0.0)
        self.assertEqual(inverse[0, 2], 0.0)
        self.assertTrue(np.linalg.eigvalsh(inverse).min() > 0)

    def test_synthetic_curve_fit_is_finite_and_positive(self):
        energy = np.linspace(6.0, 38.0, 24)
        sigma = np.sqrt(energy - 5.0) * np.exp(
            1.2 - 0.002 * (energy - 20.0) ** 2)
        frame = pd.DataFrame({
            "energy_MeV": energy,
            "campaign_id": np.repeat(["A", "B", "C"], 8),
        })
        transformed = self.model.transform(energy, sigma)
        coefficients = fit_coefficients(
            self.model, frame, transformed, np.full(len(frame), 0.05),
            smoothing=1.0, campaign_log_spread=0.1)
        prediction = self.model.predict(coefficients, energy)
        self.assertTrue(np.isfinite(prediction).all())
        self.assertTrue((prediction > 0).all())


if __name__ == "__main__":
    unittest.main()
