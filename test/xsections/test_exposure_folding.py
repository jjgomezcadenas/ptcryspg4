import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_transport.xsections.channels import CHANNELS
from analysis_transport.xsections.exposure_folding import (
    ChannelCurves,
    CrossSectionEnsemble,
    ExposureMetadata,
    accumulate_step_exposure,
    distal_r50,
    fold_exposure,
    interpolate_curves,
    load_exposure_metadata,
    validate_exposure_table,
    validate_threshold_coverage,
    write_result,
)
from common.isotopes import ISOTOPES as ISOTOPE_DATA, NAME_TO_ID
from decay_sampling.scenarios import resolve_scenario


def constant_ensemble(nominal_mb=10.0, replica_scales=(0.5, 1.5)):
    energy = np.array([0.0, 100.0])
    channels = {}
    for channel in CHANNELS:
        channels[channel.channel_id] = ChannelCurves(
            threshold_MeV=0.0,
            nominal_energy_MeV=energy,
            nominal_sigma_mb=np.full(2, nominal_mb),
            replica_energy_MeV=energy,
            replica_sigma_mb=np.asarray([
                np.full(2, nominal_mb * scale) for scale in replica_scales
            ]),
        )
    return CrossSectionEnsemble(
        channels=channels,
        replica_ids=np.arange(len(replica_scales), dtype=int),
    )


def exposure_row(
    target,
    exposure,
    *,
    energy_low=19.0,
    energy_high=21.0,
    energy_mean=20.0,
    depth_low=0.0,
    depth_high=1.0,
):
    return {
        "target": target,
        "energy_low_MeV": energy_low,
        "energy_high_MeV": energy_high,
        "energy_mean_MeV": energy_mean,
        "depth_low_mm": depth_low,
        "depth_high_mm": depth_high,
        "depth_mean_mm": 0.5 * (depth_low + depth_high),
        "target_exposure_cm2_inv": exposure,
    }


class ExposureSchemaTests(unittest.TestCase):
    def test_valid_table_is_sorted(self):
        frame = pd.DataFrame([
            exposure_row("O16", 2.0, depth_low=1.0, depth_high=2.0),
            exposure_row("C12", 1.0),
        ])
        result = validate_exposure_table(frame)
        self.assertEqual(list(result["target"]), ["C12", "O16"])

    def test_unknown_target_is_rejected(self):
        frame = pd.DataFrame([exposure_row("H1", 1.0)])
        with self.assertRaisesRegex(ValueError, "unknown targets"):
            validate_exposure_table(frame)

    def test_negative_exposure_is_rejected(self):
        frame = pd.DataFrame([exposure_row("C12", -1.0)])
        with self.assertRaisesRegex(ValueError, "non-negative"):
            validate_exposure_table(frame)

    def test_overlapping_depth_bins_are_rejected(self):
        frame = pd.DataFrame([
            exposure_row("C12", 1.0, depth_low=0.0, depth_high=1.0),
            exposure_row(
                "C12", 1.0, depth_low=0.5, depth_high=1.5,
                energy_low=21.0, energy_high=23.0, energy_mean=22.0),
        ])
        with self.assertRaisesRegex(ValueError, "depth bins overlap"):
            validate_exposure_table(frame)

    def test_overlapping_energy_bins_are_rejected_within_depth_bin(self):
        frame = pd.DataFrame([
            exposure_row(
                "C12", 1.0, energy_low=10.0,
                energy_high=20.0, energy_mean=15.0),
            exposure_row(
                "C12", 1.0, energy_low=19.0,
                energy_high=25.0, energy_mean=22.0),
        ])
        with self.assertRaisesRegex(ValueError, "energy bins overlap"):
            validate_exposure_table(frame)

    def test_step_accumulator_uses_weight_density_and_length(self):
        steps = pd.DataFrame({
            "target": ["O16", "O16"],
            "proton_weight": [1.0, 2.0],
            "target_number_density_cm3": [3.0e22, 3.0e22],
            "step_length_cm": [0.1, 0.2],
            "energy_MeV": [20.0, 22.0],
            "depth_mm": [5.2, 5.8],
        })
        result = accumulate_step_exposure(steps, [15.0, 25.0], [5.0, 6.0])
        expected = 1.0 * 3.0e22 * 0.1 + 2.0 * 3.0e22 * 0.2
        self.assertAlmostEqual(result.loc[0, "target_exposure_cm2_inv"], expected)
        self.assertAlmostEqual(
            result.loc[0, "energy_mean_MeV"],
            (3.0e21 * 20.0 + 1.2e22 * 22.0) / expected,
        )


class CrossSectionInterpolationTests(unittest.TestCase):
    def test_threshold_and_linear_interpolation(self):
        result = interpolate_curves(
            [5.0, 10.0, 15.0],
            [0.0, 20.0, 40.0],
            [4.0, 7.0, 10.0, 12.5],
            threshold_MeV=6.0,
            label="test",
        )
        np.testing.assert_allclose(result, [0.0, 8.0, 20.0, 30.0])

    def test_replica_interpolation_preserves_rows(self):
        result = interpolate_curves(
            [0.0, 10.0],
            [[0.0, 10.0], [0.0, 20.0]],
            [2.5, 7.5],
            threshold_MeV=0.0,
            label="replicas",
        )
        np.testing.assert_allclose(result, [[2.5, 7.5], [5.0, 15.0]])

    def test_extrapolation_above_fit_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "exceeds fitted limit"):
            interpolate_curves(
                [0.0, 10.0], [1.0, 1.0], [10.1],
                threshold_MeV=0.0, label="test")

    def test_repository_fit_loads_all_replicas(self):
        repository = Path(__file__).resolve().parents[2]
        ensemble = CrossSectionEnsemble.from_fit_directory(
            repository / "data/xsections/fits")
        self.assertEqual(len(ensemble.replica_ids), 1000)
        self.assertEqual(set(ensemble.channels), {
            channel.channel_id for channel in CHANNELS})
        self.assertEqual(
            ensemble.nominal("p_O16_x_N13", [5.0, 5.55]).tolist(),
            [0.0, 0.0],
        )

    def test_grid_starting_above_threshold_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "above the 5 MeV reaction threshold"):
            validate_threshold_coverage([6.0, 10.0], 5.0, "test")


class ExposureMetadataTests(unittest.TestCase):
    def test_companion_metadata_checks_hash_and_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            exposure_path = root / "exposure.csv"
            pd.DataFrame([exposure_row("C12", 1.0)]).to_csv(
                exposure_path, index=False)
            import hashlib
            digest = hashlib.sha256(exposure_path.read_bytes()).hexdigest()
            document = {
                "schema_version": 1,
                "run_id": "test",
                "exposure_file": "exposure.csv",
                "exposure_sha256": digest,
                "n_protons": 100.0,
                "target_dose_Gy": 2.0,
                "Np_per_Gy": 50.0,
                "physics_list": "QGSP_BIC_HP",
                "geant4_version": "test",
                "random_seed": 7,
                "software_revision": "abc",
                "beam_axis": "z",
                "depth_origin": "phantom entrance",
                "depth_unit": "mm",
                "energy_edges_MeV": [19.0, 21.0],
                "depth_edges_mm": [0.0, 1.0],
            }
            meta_path = root / "exposure_meta.json"
            meta_path.write_text(json.dumps(document))
            metadata = load_exposure_metadata(meta_path, exposure_path)
            self.assertEqual(metadata.Np_per_Gy, 50.0)
            document["Np_per_Gy"] = 40.0
            meta_path.write_text(json.dumps(document))
            with self.assertRaisesRegex(ValueError, "inconsistent"):
                load_exposure_metadata(meta_path, exposure_path)


class ExposureFoldingTests(unittest.TestCase):
    def setUp(self):
        self.ensemble = constant_ensemble()

    def test_uniform_slab_matches_number_density_path_length_formula(self):
        number_density_cm3 = 2.0e22
        total_weighted_path_cm = 5.0
        exposure = number_density_cm3 * total_weighted_path_cm
        frame = pd.DataFrame([exposure_row("C12", exposure)])
        result = fold_exposure(frame, self.ensemble)
        contribution = result.nominal_channel_contributions.iloc[0]
        expected = exposure * 10.0e-27
        self.assertAlmostEqual(contribution["expected_nuclei_run"], expected)

    def test_two_target_channels_add_to_c11(self):
        frame = pd.DataFrame([
            exposure_row("C12", 1.0e27),
            exposure_row("O16", 2.0e27),
        ])
        result = fold_exposure(frame, self.ensemble)
        profile = result.nominal_isotope_profiles
        c11 = profile.loc[
            profile["profile_label"] == "C11", "expected_count_run"].iloc[0]
        o15 = profile.loc[
            profile["profile_label"] == "O15", "expected_count_run"].iloc[0]
        n13 = profile.loc[
            profile["profile_label"] == "N13", "expected_count_run"].iloc[0]
        self.assertAlmostEqual(c11, 30.0)
        self.assertAlmostEqual(o15, 20.0)
        self.assertAlmostEqual(n13, 20.0)
        self.assertEqual(
            set(profile["profile_label"]),
            {"C11", "O15", "N13", "all_production", "all_d120s300"},
        )

    def test_scenario_aggregate_uses_named_handoff_factors(self):
        frame = pd.DataFrame([
            exposure_row("C12", 1.0e27),
            exposure_row("O16", 2.0e27),
        ])
        result = fold_exposure(frame, self.ensemble)
        profile = result.nominal_isotope_profiles.set_index("profile_label")
        scenario = resolve_scenario("d120s300")
        expected = sum(
            profile.loc[name, "expected_count_run"]
            * scenario.measured_fraction(ISOTOPE_DATA[NAME_TO_ID[name]].lam)
            for name in ("C11", "O15", "N13")
        )
        self.assertAlmostEqual(
            profile.loc["all_d120s300", "expected_count_run"], expected)
        self.assertEqual(profile.loc["all_d120s300", "quantity"], "measured_decays")

    def test_replicas_scale_same_frozen_exposure(self):
        rows = []
        for depth, exposure in enumerate([1.0e27, 1.0e27, 0.2e27]):
            rows.append(exposure_row(
                "C12", exposure, depth_low=float(depth),
                depth_high=float(depth + 1)))
        result = fold_exposure(pd.DataFrame(rows), self.ensemble)
        summary = result.production_summary
        nominal = summary[(summary["model"] == "nominal")
                          & (summary["profile_label"] == "C11")].iloc[0]
        replicas = summary[(summary["model"] == "replica")
                           & (summary["profile_label"] == "C11")]
        np.testing.assert_allclose(
            replicas["expected_count_run"],
            nominal["expected_count_run"] * np.array([0.5, 1.5]),
        )
        np.testing.assert_allclose(replicas["R50_shift_mm"], [0.0, 0.0])

    def test_constant_cross_section_is_invariant_under_energy_bin_refinement(self):
        coarse = pd.DataFrame([exposure_row("C12", 2.0e27)])
        fine = pd.DataFrame([
            exposure_row(
                "C12", 0.8e27, energy_low=19.0,
                energy_high=20.0, energy_mean=19.5),
            exposure_row(
                "C12", 1.2e27, energy_low=20.0,
                energy_high=21.0, energy_mean=20.5),
        ])
        coarse_result = fold_exposure(coarse, self.ensemble)
        fine_result = fold_exposure(fine, self.ensemble)
        coarse_yield = coarse_result.nominal_channel_contributions[
            "expected_nuclei_run"].sum()
        fine_yield = fine_result.nominal_channel_contributions[
            "expected_nuclei_run"].sum()
        self.assertAlmostEqual(coarse_yield, fine_yield)

    def test_binned_fold_matches_direct_step_sum_for_linear_curve(self):
        steps = pd.DataFrame({
            "target": ["C12", "C12"],
            "proton_weight": [1.0, 1.0],
            "target_number_density_cm3": [2.0e22, 2.0e22],
            "step_length_cm": [0.1, 0.3],
            "energy_MeV": [10.0, 30.0],
            "depth_mm": [0.4, 0.6],
        })
        exposure = accumulate_step_exposure(
            steps, [0.0, 40.0], [0.0, 1.0])
        energy = np.array([0.0, 100.0])
        channels = {}
        for channel in CHANNELS:
            channels[channel.channel_id] = ChannelCurves(
                threshold_MeV=0.0,
                nominal_energy_MeV=energy,
                nominal_sigma_mb=2.0 * energy,
                replica_energy_MeV=energy,
                replica_sigma_mb=np.asarray([2.0 * energy]),
            )
        ensemble = CrossSectionEnsemble(channels, np.array([0]))
        result = fold_exposure(exposure, ensemble)
        folded = result.nominal_channel_contributions.loc[
            result.nominal_channel_contributions["channel_id"]
            == "p_C12_x_C11", "expected_nuclei_run"].sum()
        direct = np.sum(
            steps["proton_weight"]
            * steps["target_number_density_cm3"]
            * steps["step_length_cm"]
            * (2.0 * steps["energy_MeV"])
            * 1.0e-27
        )
        self.assertAlmostEqual(folded, direct)

    def test_repository_replicas_fold_on_synthetic_exposure(self):
        repository = Path(__file__).resolve().parents[2]
        ensemble = CrossSectionEnsemble.from_fit_directory(
            repository / "data/xsections/fits")
        rows = []
        for target in ("C12", "N14", "O16"):
            for depth, scale in enumerate((1.0, 0.8, 0.2)):
                rows.append(exposure_row(
                    target, scale * 1.0e24,
                    energy_low=49.0, energy_high=51.0, energy_mean=50.0,
                    depth_low=float(depth), depth_high=float(depth + 1),
                ))
        result = fold_exposure(pd.DataFrame(rows), ensemble)
        self.assertEqual(len(result.replica_isotope_profiles), 1000 * 5 * 3)
        self.assertEqual(len(result.production_summary), (1000 + 1) * 5)
        self.assertTrue(
            np.isfinite(result.production_summary["expected_count_run"]).all())
        self.assertEqual(len(result.profile_bands), 5 * 3)
        self.assertEqual(len(result.uncertainty_summary), 5)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            write_result(result, output, cross_sections=ensemble)
            metadata = json.loads((output / "folding_meta.json").read_text())
            self.assertEqual(metadata["n_replicas"], 1000)
            self.assertEqual(len(metadata["fit_files_sha256"]), 10)

    def test_result_writer_creates_complete_product(self):
        frame = pd.DataFrame([
            exposure_row("C12", 1.0e27, depth_low=0.0, depth_high=1.0),
            exposure_row("C12", 0.2e27, depth_low=1.0, depth_high=2.0),
        ])
        result = fold_exposure(frame, self.ensemble)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            write_result(result, output, cross_sections=self.ensemble)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {
                    "nominal_channel_contributions.csv",
                    "nominal_isotope_profiles.csv",
                    "replica_isotope_profiles.csv.gz",
                    "profile_bands.csv",
                    "production_summary.csv",
                    "uncertainty_summary.csv",
                    "folding_meta.json",
                },
            )

    def test_result_writer_can_omit_replica_profiles(self):
        result = fold_exposure(
            pd.DataFrame([exposure_row("C12", 1.0e27)]), self.ensemble)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            write_result(
                result, output, cross_sections=self.ensemble,
                write_replica_profiles=False)
            self.assertFalse((output / "replica_isotope_profiles.csv.gz").exists())

    def test_result_writer_records_scenario_and_native_route_fraction(self):
        result = fold_exposure(
            pd.DataFrame([exposure_row("C12", 1.0e27)]), self.ensemble)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            route_path = output / "native_route_summary.csv"
            pd.DataFrame([
                {
                    "profile_label": "all_production",
                    "modeled_count": 90.0,
                    "unmodeled_count": 10.0,
                    "total_count": 100.0,
                    "unmodeled_fraction": 0.10,
                },
                {
                    "profile_label": "all_d120s300",
                    "modeled_count": 45.0,
                    "unmodeled_count": 2.0,
                    "total_count": 47.0,
                    "unmodeled_fraction": 2.0 / 47.0,
                },
            ]).to_csv(route_path, index=False)
            scenario = resolve_scenario("d120s300")
            write_result(
                result,
                output / "fold",
                cross_sections=self.ensemble,
                scenario=scenario,
                native_route_summary_path=route_path,
                write_replica_profiles=False,
            )
            metadata = json.loads(
                (output / "fold/folding_meta.json").read_text())
            self.assertEqual(metadata["handoff_scenario"]["name"], "d120s300")
            self.assertFalse(
                metadata["native_route_diagnostic"]["affects_folded_source"])
            self.assertEqual(
                len(metadata["native_route_diagnostic"]["fractions"]), 2)

    def test_run_per_proton_and_per_gy_normalizations_are_consistent(self):
        metadata = ExposureMetadata.synthetic()
        metadata = ExposureMetadata(
            **{
                **metadata.__dict__,
                "n_protons": 10.0,
                "target_dose_Gy": 2.0,
                "Np_per_Gy": 5.0,
            }
        )
        result = fold_exposure(
            pd.DataFrame([exposure_row("C12", 1.0e27)]),
            self.ensemble,
            metadata,
        )
        row = result.nominal_isotope_profiles.loc[
            result.nominal_isotope_profiles["profile_label"] == "C11"].iloc[0]
        self.assertAlmostEqual(row["expected_count_per_proton"], row["expected_count_run"] / 10.0)
        self.assertAlmostEqual(row["expected_count_per_Gy"], row["expected_count_run"] / 2.0)


class DistalEdgeTests(unittest.TestCase):
    def test_linear_distal_crossing(self):
        edge = distal_r50([0.5, 1.5, 2.5], [1.0, 1.0, 0.2])
        self.assertAlmostEqual(edge, 2.125)

    def test_profile_without_distal_crossing_returns_nan(self):
        self.assertTrue(np.isnan(distal_r50([0.5, 1.5], [1.0, 1.0])))


if __name__ == "__main__":
    unittest.main()
