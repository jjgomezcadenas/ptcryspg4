import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_transport.xsections.exposure_folding import (
    accumulate_step_exposure,
)

REPO = Path(__file__).resolve().parents[2]
EXECUTABLE = REPO / "stageA_xsection_ensemble/build/xsection_ensemble"


@unittest.skipUnless(EXECUTABLE.exists(),
                     "xsection_ensemble executable not built")
class StepAccumulatorEqualityTests(unittest.TestCase):
    """The C++ scorer and the Python reference accumulator must produce the
    same exposure table from the same step list (phase 1a acceptance 2)."""

    def test_cxx_matches_python_reference(self):
        rng = np.random.default_rng(11)
        n = 2000
        steps = pd.DataFrame({
            "target": rng.choice(["C12", "N14", "O16"], n),
            "proton_weight": rng.uniform(0.5, 2.0, n),
            "target_number_density_cm3": rng.uniform(1e22, 9e22, n),
            "step_length_cm": rng.uniform(1e-3, 0.1, n),
            "energy_MeV": rng.uniform(0.1, 129.9, n),
            "depth_mm": rng.uniform(0.0, 199.9, n),
        })
        with tempfile.TemporaryDirectory() as directory:
            steps_path = Path(directory) / "steps.csv"
            out_path = Path(directory) / "cxx.csv"
            steps.to_csv(steps_path, index=False)
            subprocess.run(
                [str(EXECUTABLE), "--steps-csv", str(steps_path),
                 "--steps-out", str(out_path), "--zmax-mm", "200",
                 "--ebin-width-mev", "0.5", "--emax-mev", "130",
                 "--zbin-width-mm", "1.0"],
                check=True, capture_output=True)
            cxx = pd.read_csv(out_path)
        python = accumulate_step_exposure(
            steps, np.arange(0.0, 130.001, 0.5), np.arange(0.0, 200.001, 1.0))
        key = ["target", "energy_low_MeV", "depth_low_mm"]
        merged = python.merge(cxx, on=key, suffixes=("_py", "_cxx"))
        self.assertEqual(len(merged), len(python))
        self.assertEqual(len(merged), len(cxx))
        for column in ("target_exposure_cm2_inv", "energy_mean_MeV",
                       "depth_mean_mm"):
            np.testing.assert_allclose(
                merged[f"{column}_py"], merged[f"{column}_cxx"],
                rtol=1.0e-10,
                err_msg=f"{column} differs between C++ and Python")


if __name__ == "__main__":
    unittest.main()
