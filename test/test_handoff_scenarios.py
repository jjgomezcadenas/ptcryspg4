import unittest

from decay_sampling.scenarios import resolve_scenario


class HandoffScenarioTests(unittest.TestCase):
    def test_named_inroom_scenario_has_frozen_times(self):
        scenario = resolve_scenario("inroom")
        self.assertEqual(
            (scenario.t_irr_s, scenario.t_del_s, scenario.t_meas_s),
            (60.0, 120.0, 1200.0),
        )
        self.assertEqual(len(scenario.config_sha256), 64)

    def test_unknown_scenario_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unknown handoff scenario"):
            resolve_scenario("not-a-scenario")


if __name__ == "__main__":
    unittest.main()
