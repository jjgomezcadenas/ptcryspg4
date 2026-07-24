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

One row defines an interpolation-ready curve used by the reweighting code:

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
evaluations into a reweighting curve.
