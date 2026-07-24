import unittest

from analysis_transport.xsections.reaction_thresholds import (
    REACTIONS, laboratory_threshold_mev, q_value_mev,
)


class ReactionThresholdTests(unittest.TestCase):
    def test_oxygen_15_lowest_channel_threshold(self):
        reaction = next(r for r in REACTIONS if r.channel_id == "p_O16_x_O15")
        self.assertAlmostEqual(q_value_mev(reaction), -13.4394, places=3)
        self.assertAlmostEqual(laboratory_threshold_mev(reaction), 14.292, places=3)

    def test_carbon_11_lowest_channel_threshold(self):
        reaction = next(r for r in REACTIONS if r.channel_id == "p_O16_x_C11")
        self.assertAlmostEqual(q_value_mev(reaction), -23.6581, places=3)
        self.assertAlmostEqual(laboratory_threshold_mev(reaction), 25.168, places=3)


if __name__ == "__main__":
    unittest.main()
