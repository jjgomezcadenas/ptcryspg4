import unittest

import numpy as np
import pandas as pd

from analysis_transport.xsections.channels import CHANNELS
from analysis_transport.xsections.exposure_convergence import (
    coarsen_exposure,
    compare_folding_results,
)
from analysis_transport.xsections.exposure_folding import (
    ChannelCurves,
    CrossSectionEnsemble,
    fold_exposure,
)


def exposure_row(energy_low, energy_high, energy_mean, exposure, depth_low=0.0):
    return {
        "target": "C12",
        "energy_low_MeV": energy_low,
        "energy_high_MeV": energy_high,
        "energy_mean_MeV": energy_mean,
        "depth_low_mm": depth_low,
        "depth_high_mm": depth_low + 1.0,
        "depth_mean_mm": depth_low + 0.5,
        "target_exposure_cm2_inv": exposure,
    }


def quadratic_ensemble():
    energy = np.linspace(0.0, 100.0, 1001)
    channels = {}
    for channel in CHANNELS:
        nominal = energy**2
        replicas = np.asarray([
            energy**2,
            energy**2 * (1.0 + 0.004 * energy),
            energy**2 * (1.0 - 0.003 * energy),
        ])
        channels[channel.channel_id] = ChannelCurves(
            threshold_MeV=0.0,
            nominal_energy_MeV=energy,
            nominal_sigma_mb=nominal,
            replica_energy_MeV=energy,
            replica_sigma_mb=replicas,
        )
    return CrossSectionEnsemble(channels, np.arange(3))


class ExposureCoarseningTests(unittest.TestCase):
    def test_coarsening_preserves_exposure_and_first_moments(self):
        fine = pd.DataFrame([
            exposure_row(0.0, 10.0, 4.0, 2.0e20),
            exposure_row(10.0, 20.0, 16.0, 3.0e20),
        ])
        coarse = coarsen_exposure(fine, [0.0, 20.0])
        self.assertEqual(len(coarse), 1)
        self.assertAlmostEqual(coarse.loc[0, "target_exposure_cm2_inv"], 5.0e20)
        self.assertAlmostEqual(coarse.loc[0, "energy_mean_MeV"], 11.2)
        self.assertAlmostEqual(coarse.loc[0, "depth_mean_mm"], 0.5)

    def test_fine_bin_straddling_coarse_boundary_is_rejected(self):
        fine = pd.DataFrame([exposure_row(5.0, 15.0, 10.0, 1.0)])
        with self.assertRaisesRegex(ValueError, "does not lie wholly"):
            coarsen_exposure(fine, [0.0, 10.0, 20.0])

    def test_nonlinear_cross_section_exposes_mean_energy_bias(self):
        ensemble = quadratic_ensemble()
        reference = pd.DataFrame([
            exposure_row(0.0, 10.0, 5.0, 1.0e27),
            exposure_row(10.0, 20.0, 15.0, 1.0e27),
        ])
        coarse = coarsen_exposure(reference, [0.0, 20.0])
        reference_yield = fold_exposure(
            reference, ensemble).nominal_channel_contributions["expected_nuclei_run"].sum()
        coarse_yield = fold_exposure(
            coarse, ensemble).nominal_channel_contributions["expected_nuclei_run"].sum()
        self.assertGreater(abs(reference_yield - coarse_yield), 0.0)

    def test_identical_foldings_pass_replica_relative_convergence(self):
        ensemble = quadratic_ensemble()
        rows = []
        for depth, exposure in enumerate((1.0e27, 1.0e27, 0.3e27)):
            rows.append(exposure_row(
                10.0 + depth, 11.0 + depth, 10.5 + depth,
                exposure, depth_low=float(depth)))
        result = fold_exposure(pd.DataFrame(rows), ensemble)
        comparison = compare_folding_results(result, result)
        self.assertTrue(comparison["pass"].all())


if __name__ == "__main__":
    unittest.main()
