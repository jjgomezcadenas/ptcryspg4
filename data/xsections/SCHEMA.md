# Cross-section data schema

## Stable identifiers

A normalized dataset has a stable identifier of the form

```text
<library>_<accession>_<subentry>_p_<target>_x_<residual>
```

Examples:

```text
exfor_E2449_002_p_O16_x_O15
tendl2023_p_O16_x_O15
jendl40he_p_O16_x_O15
```

## `sources.csv`

One row identifies an upstream snapshot or downloaded target file:

| Column | Meaning |
|---|---|
| `source_id` | Stable local source identifier |
| `library` | EXFOR, TENDL or JENDL |
| `release` | Upstream release or snapshot date |
| `scope` | Reactions or target covered by the file |
| `source_url` | Retrieval URL |
| `retrieval_date` | Local retrieval date in ISO format |
| `revision` | Upstream commit or evaluation revision |
| `sha256` | SHA-256 digest for a single downloaded file |
| `local_path` | Path relative to `data/xsections/` |
| `terms` | Upstream licence or distribution terms |

For the selected EXFOR tree, the upstream Git commit fixes the contents of all
individual CX4 files.

## Normalized point table

Each file in `normalized/` represents one experimental or evaluated curve:

| Column | Unit | Meaning |
|---|---:|---|
| `point_id` | -- | Row identifier within the dataset |
| `energy_MeV` | MeV | Incident proton kinetic energy |
| `energy_unc_minus_MeV` | MeV | Lower energy uncertainty |
| `energy_unc_plus_MeV` | MeV | Upper energy uncertainty |
| `sigma_mb` | mb | Residual-production cross section |
| `sigma_unc_stat_mb` | mb | Statistical uncertainty |
| `sigma_unc_sys_mb` | mb | Systematic uncertainty |
| `sigma_unc_minus_mb` | mb | Lower total uncertainty |
| `sigma_unc_plus_mb` | mb | Upper total uncertainty |

An unavailable uncertainty component is empty. Original units and uncertainty
labels remain in the raw file and are recorded by the normalization metadata.

## EXFOR curation tables

`curation.csv` contains one decision for every normalized EXFOR series. Its
principal columns are:

| Column | Meaning |
|---|---|
| `dataset_id` | Stable normalized dataset identifier |
| `channel_id` | Fitted production channel |
| `reaction`, `quantity` | Original EXFOR interpretation fields |
| `state` | `accepted`, `excluded`, or `pending` |
| `reason_code` | Machine-readable decision reason |
| `verification_basis` | Metadata or publication field supporting the decision |
| `notes` | Required follow-up or explanatory detail |
| `points_passing_point_rules` | Points eligible for the weighted fit |

`point_curation.csv` records the corresponding decision for every experimental
point. `include_in_fit` is one only when the dataset is accepted and the point
passes the configured energy, threshold, positivity, and quoted-uncertainty
rules. `reason_code` records the first failed rule. `curation_meta.json` pins
the curation policy, fit configuration, normalized catalog, and row counts by
SHA-256 digest.

## EXFOR fit products

Files under `fits/` are keyed by `channel_id`. A `*_curve.csv` file contains
the dense-grid median, 16th and 84th percentiles, and relative band half-width.
The common-grid values and all replica samples are in `*_table.csv` and
`*_replicas.csv`. `*_representatives.csv` stores the nine selected full curves;
`*_distance_histogram.csv` stores the binned covariance-normalized distance.
`fit_summary.csv` and `fit_meta.json` record channel metrics and configuration
and curation digests.

`p_C12_x_C11_sensitivity.csv` is the central fit obtained by adding pending
B0095.002 to the accepted C-12 data. `sensitivity_summary.csv` records the
added dataset, fitted-point count, selected hyperparameters, peak, and
fractional difference from the nominal fit. It is a curation sensitivity and
is not part of the replica ensemble.

`reaction_thresholds.csv` records the lowest production channel, Q value,
exact laboratory-frame threshold, and AME2020 mass provenance for the two
oxygen-channel threshold corrections. The fit configuration stores rounded
values and the generator checks them against this calculation.

## Covariance table

Evaluated covariance data use a long-form table:

| Column | Unit | Meaning |
|---|---:|---|
| `point_i` | -- | First normalized grid index |
| `point_j` | -- | Second normalized grid index |
| `energy_i_MeV` | MeV | First grid energy |
| `energy_j_MeV` | MeV | Second grid energy |
| `covariance_mb2` | mb2 | Cross-section covariance |

## `models/models.csv`

One row defines an interpolation-ready production curve:

| Column | Meaning |
|---|---|
| `model_id` | Stable model identifier |
| `channel_id` | Projectile, target and residual identifier |
| `curve_file` | Model point-table path |
| `parent_dataset_ids` | Source datasets used to construct the curve |
| `construction` | Selection, fit, join or envelope rule |
| `energy_min_MeV` | Lower valid energy |
| `energy_max_MeV` | Upper valid energy |
| `interpolation` | Interpolation rule within the valid range |
| `extrapolation` | Rule outside the valid range |
| `version` | Local model version |

The model catalog records every scientific choice used to turn measurements or
evaluations into an interpolation-ready production curve.

## Target-resolved proton exposure table

The dedicated Geant4 transport run accumulates one row for every target
isotope, proton-energy bin and beam-depth bin. The table is the transport input
to the direct cross-section folding calculation.

| Column | Unit | Meaning |
|---|---:|---|
| `target` | -- | `C12`, `N14`, or `O16` |
| `energy_low_MeV` | MeV | Lower edge of the proton-energy bin |
| `energy_high_MeV` | MeV | Upper edge of the proton-energy bin |
| `energy_mean_MeV` | MeV | Exposure-weighted mean proton energy in the bin |
| `depth_low_mm` | mm | Lower edge of the beam-depth bin |
| `depth_high_mm` | mm | Upper edge of the beam-depth bin |
| `depth_mean_mm` | mm | Exposure-weighted mean depth in the bin |
| `target_exposure_cm2_inv` | cm^-2 | Sum of `proton_weight * target_number_density_cm3 * step_length_cm` |

Multiplication of `target_exposure_cm2_inv` by a channel cross section in
square centimetres gives the expected number of residual nuclei produced in
that bin. The table contains the summed exposure of the simulated run; the run
metadata provide the primary-proton count, dose, physics list, Geant4 version,
random seed, energy binning, depth binning and coordinate definition.

The first implementation is longitudinal because its immediate observable is
the production-depth edge. A separate proposal sample of proton step segments
will provide three-dimensional production sites for positron transport.

## Direct-folding products

`analysis_transport/xsections/exposure_folding.py` combines the exposure table
with the EXFOR nominal curves and all cross-section replicas. It writes:

| File | Meaning |
|---|---|
| `nominal_channel_contributions.csv` | Nominal cross section and expected production for every exposure row and channel |
| `nominal_isotope_profiles.csv` | Nominal depth profiles for C-11, N-13, O-15 and their sum |
| `replica_isotope_profiles.csv` | Corresponding profile for every replica |
| `production_summary.csv` | Isotope budgets, production R50 and replica displacement from nominal |
| `folding_meta.json` | Schema version, unit conversion, input digests and fit provenance |

The conversion used by the folding calculation is `1 mb = 1e-27 cm2`.

## Geant4 thin-target tables

The versioned physics-list directory under g4 contains pilot_counts.csv,
which records every tested target thickness and marks the selected row.
thin_target_counts.csv contains one aggregate production row per target and
incident energy. Its counters, geometry fields, mean interaction energies,
seed, thread count, physics list and Geant4 version are direct simulation
outputs. run_meta.json records the scan-configuration digest and runtime
limits. The conversion into mb is performed by
analysis_transport/xsections/g4_cross_sections.py.

## Geant4 direct denominator tables

`denominator_pilot_counts.csv` contains one row per pure target and incident
energy. Each row records a fixed number of forced proton-inelastic final-state
samples, the residual counts, the queried inelastic cross section, the Geant4
cross-section data-set name, final-state model counts, version, physics list,
and seed. `denominator_pilot_meta.json` pins the executable and configuration
digests.

The `denominator_pilot/` subdirectory contains the factorized channel curves,
where `sigma_mb` is the queried inelastic cross section multiplied by the
sampled residual probability. Zero counts carry a 95% binomial upper limit.
`support_summary.csv` audits nonzero Geant4 support over the production-bearing
interval of each experimental fit.
