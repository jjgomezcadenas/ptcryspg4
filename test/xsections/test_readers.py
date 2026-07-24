import tempfile
import unittest
from pathlib import Path

from analysis_transport.xsections import endf6, exfor, iaea, jendl
from analysis_transport.xsections.normalize import collect


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data/xsections/raw"


class ReaderTests(unittest.TestCase):
    def test_endf_float(self):
        self.assertEqual(endf6.endf_float(" 1.250000+3"), 1250.0)
        self.assertAlmostEqual(endf6.endf_float("-2.500000-2"), -0.025)

    def test_exfor_units_and_identity(self):
        path = next((RAW / "exfor/2026-06-29/cx4/p/O16/x/sig").glob("*prodO15_Akagi*.cx4"))
        dataset = exfor.read(path)
        self.assertEqual((dataset["target"], dataset["residual"]), ("O16", "O15"))
        self.assertAlmostEqual(dataset["points"][0]["energy_MeV"], 7.2)
        self.assertAlmostEqual(dataset["points"][0]["sigma_mb"], 1.2)

    def test_jendl_block(self):
        dataset = jendl.read(RAW / "jendl/4.0he/O016.txt", "N13")
        self.assertEqual(dataset["target"], "O16")
        self.assertGreater(max(point["sigma_mb"] for point in dataset["points"]), 1.0)

    def test_iaea_table(self):
        dataset = iaea.read(RAW / "iaea/medical/o6p13nt.txt")
        self.assertEqual(dataset["points"][0]["energy_MeV"], 6.5)
        self.assertEqual(dataset["points"][0]["sigma_unc_plus_mb"], 0.003)

    def test_lanl_endf_residual(self):
        dataset = endf6.read(RAW / "tendl/2023/p-O016.tendl", "O16", "N13", 7, 13)
        self.assertEqual(dataset["residual"], "N13")
        self.assertGreater(max(point["sigma_mb"] for point in dataset["points"]), 1.0)

    def test_all_dataset_ids_unique(self):
        datasets = collect(REPO)
        identifiers = [dataset["dataset_id"] for dataset in datasets]
        self.assertEqual(len(identifiers), len(set(identifiers)))
        self.assertEqual(sum(dataset["library"] == "EXFOR" for dataset in datasets), 85)


if __name__ == "__main__":
    unittest.main()
