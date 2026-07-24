import unittest
from pathlib import Path

from analysis_transport.xsections.channels import CHANNELS
from analysis_transport.xsections.make_comparison import _coverage, load_catalog
from analysis_transport.xsections.validate_comparison import validate


REPO = Path(__file__).resolve().parents[2]


class ComparisonTests(unittest.TestCase):
    def test_every_channel_has_external_coverage(self):
        coverage = _coverage(load_catalog(REPO))
        covered = {row["channel_id"] for row in coverage}
        self.assertEqual(covered, {channel.channel_id for channel in CHANNELS})

    def test_generated_products_validate(self):
        self.assertEqual(validate(REPO), 96)


if __name__ == "__main__":
    unittest.main()
