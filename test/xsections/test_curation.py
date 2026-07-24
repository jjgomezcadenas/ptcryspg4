import unittest
from pathlib import Path

import pandas as pd

from analysis_transport.xsections.curate_exfor import reaction_suffix


REPO = Path(__file__).resolve().parents[2]


class CurationTests(unittest.TestCase):
    def test_reaction_suffix(self):
        self.assertEqual(
            reaction_suffix("6-C-12(P,X)6-C-11,,SIG"), ",,SIG")
        self.assertEqual(
            reaction_suffix("6-C-12(P,X)6-C-11,IND,SIG,,A,EXP"),
            ",IND,SIG,,A,EXP")

    def test_every_exfor_series_has_one_decision(self):
        catalog = pd.read_csv(REPO / "data/xsections/normalized/datasets.csv")
        expected = set(catalog.loc[catalog.library == "EXFOR", "dataset_id"])
        curation = pd.read_csv(REPO / "data/xsections/curation.csv")
        self.assertFalse(curation.dataset_id.duplicated().any())
        self.assertEqual(set(curation.dataset_id), expected)

    def test_abundance_weighted_series_are_pending(self):
        curation = pd.read_csv(REPO / "data/xsections/curation.csv")
        weighted = curation[curation.reaction.str.contains(",A,", regex=False)]
        self.assertEqual(len(weighted), 2)
        self.assertTrue((weighted.state == "pending").all())
        self.assertFalse(
            curation.loc[curation.state == "accepted", "reaction"]
            .str.contains(",A,", regex=False).any())

    def test_named_holdouts_are_pending(self):
        curation = pd.read_csv(REPO / "data/xsections/curation.csv")
        expected = {
            "B0095.002": "incident_energy_field_unclassified",
            "E2568.002": "shared_external_normalization_not_modelled",
            "E2568.003": "shared_external_normalization_not_modelled",
            "E2568.004": "shared_external_normalization_not_modelled",
        }
        selected = curation[curation.accession.isin(expected)]
        self.assertEqual(set(selected.accession), set(expected))
        self.assertTrue((selected.state == "pending").all())
        self.assertEqual(
            dict(zip(selected.accession, selected.reason_code)), expected)


if __name__ == "__main__":
    unittest.main()
