# Interface contracts

**This file is the authoritative contract** for everything shared between the
stages of the simulation chain: the isotope encoding, the units convention, and
the layout of the interface files. The code mirrors it (`common/Isotopes.hh` for
C++, `common/isotopes.py` for Python) and every reader/writer must conform. If
code and this file disagree, the code is wrong.

Spec reference: `docs/simulate_pt_pet.tex` (§2 Stage A, §3 handoff, §4 Stage B).

**File format: CSV.** Interface files are flat, columnar **CSV** (one row per
record, one scalar per column, header row, units as `_mm`/`_keV`/`_ns` column
suffixes). Each data file has a companion `*_meta.csv` (a single wide row)
holding run-level metadata. CSV was chosen over HDF5 because the files are small
at realistic statistics and it needs no extra dependency (pandas-native,
human-readable). **HDF5 is deferred** — if size or read speed ever demands it,
switch to HighFive (C++) / h5py (Python); the columnar schema below is identical,
so only the container changes.

---

## Units convention

| Quantity | Unit | Notes |
|----------|------|-------|
| position | mm | |
| energy | keV | beam energy is quoted in MeV |
| time | ns | |
| dose | Gy | |

Columns carry the unit as a suffix (e.g. `prod_x_mm`); the companion
`*_meta.csv` records the rest.

---

## Isotope encoding

| `isotope_id` | isotope | Z | A | T½ (s) | β⁺ endpoint (MeV) | prompt γ |
|---|---|---|---|---|---|---|
| 0 | ¹⁵O | 8 | 15 | 122.24 | 1.74 | no |
| 1 | ¹¹C | 6 | 11 | 1223.4 | 0.96 | no |
| 2 | ¹³N | 7 | 13 | 597.9 | 1.19 | no |
| 3 | ¹⁰C | 6 | 10 | 19.29 | 1.91 | yes (718 keV) |
| 4 | ¹⁴O | 8 | 14 | 70.62 | 1.81 | yes |

- 100 % β⁺ branching assumed in the bookkeeping (spec §3); Geant4 samples the
  true β⁺/EC split, so captured positrons ≈ production × β⁺-BR (≳99 %).
- Half-lives: tabulated nuclear data, consistent with spec Table 2 / Parodi 2008.
- Prompt γ's are **discarded** in Stage A (spec §2.4, "backgrounds discarded").
- Per-isotope arrays (e.g. `Nj_budget`) expand to 5 columns suffixed by
  `isotope_id` (`_0`…`_4`), ordered by the table above.

---

## File 1: Stage A output — `emitters.csv` (+ `run_meta.csv`, `depth_dose.csv`)

Written once by `stageA_transport`. The detector-independent source.

### `emitters.csv` — one row per β⁺ emitter

| column | type | meaning |
|--------|------|---------|
| `event_id` | int | primary (event) that produced it — diagnostic, multiplicity |
| `isotope_id` | int8 | see encoding table |
| `prod_x_mm`, `prod_y_mm`, `prod_z_mm` | float32 | production/decay point — the **truth** activity map |
| `anh_x_mm`, `anh_y_mm`, `anh_z_mm` | float32 | positron annihilation point — the detector source |

Positron range per event = `anh − prod` (no separate kernel needed). Positrons
that would escape the phantom into air are killed at the boundary, so their
`anh` is pinned to the phantom surface (keep-at-surface; ~0.6 % of rows).

### `run_meta.csv` — one wide row of run metadata

| column | type | meaning |
|--------|------|---------|
| `n_protons` | int | primaries run |
| `beam_energy_MeV` | float | primary proton kinetic energy |
| `beam_sigma_mm` | float | pencil-beam Gaussian σ at entrance |
| `phantom_material` | string | `"G4_PLEXIGLASS"` (PMMA) |
| `phantom_diameter_mm`, `phantom_length_mm` | float | 200.0 baseline |
| `phantom_mass_g` | float | scoring-volume mass (for dose) |
| `edep_total_MeV` | float | total energy deposited in the phantom |
| `dose_total_Gy` | float | whole-phantom dose for `n_protons` |
| `geant4_version` | string | e.g. `"11.4.1"` |
| `physics_list` | string | e.g. `"QGSP_BIC_HP"` |
| `random_seed` | int | master seed of the run |

`Pj_per_proton` is derivable as per-isotope count / `n_protons`. Dose is
whole-phantom; a Bragg-region "dose at target" definition is TBD (spec §2.5).

### `depth_dose.csv` — Bragg profile (1 mm z-bins)

| column | type | meaning |
|--------|------|---------|
| `z_mm` | float | bin centre along the beam axis |
| `edep_total_MeV` | float | energy deposit in the bin (all particles) |
| `edep_primary_MeV` | float | energy deposit by the primary proton only |

### Invariants
- Generated **once**; never modified by downstream stages.
- Encodes the production **shape** P_j(r) at high statistics; the absolute
  measured count is set downstream (File 2), never here.

---

## File 2: handoff output — `sampling_budget_<scenario>.csv`

Written by `decay_sampling/budget.py` (Python). The annihilation events are
**not materialized** — the handoff writes only the per-isotope counts (the
budget); Stage B draws that many annihilation points from File 1 on the fly,
with seed `master_seed + realization` (so every detector gets the identical
source). Method: `docs/handoff.tex`.

### Columns — one row per (realization, isotope)

| column | type | meaning |
|--------|------|---------|
| `realization` | int | realization index (0 … Z−1) |
| `isotope_id` | int8 | which species (encoding above) |
| `N_expected` | float | measured decays from Eq. 1, `P_j(D)·survival_j` |
| `N_poisson` | int | this realization's draw `M_j ~ Poisson(N_expected)` |

### Companion `sampling_budget_<scenario>_meta.csv` (one wide row)

`scenario`, `source_file` (the `emitters.csv` drawn from), `dose_Gy`,
`t_irr_s`, `t_del_s`, `t_meas_s`, `n_realizations` (Z), `master_seed`,
`target_dose_Gy`.

### Invariants
- N_j sampling (the budget) is set **only here** (CLAUDE.md invariant 3).
- Every detector config consumes the **identical** source — guaranteed by the
  same `emitters.csv` + budget + seed (CLAUDE.md invariant 4).

---

## File 3: Stage B output — `coincidences_<config>.csv` (the deliverable)

Written by `stageB_detector`, one file per (detector config, realization).

### Columns — one row per accepted coincidence

| column | type | meaning |
|--------|------|---------|
| `hit1_x_mm`, `hit1_y_mm`, `hit1_z_mm` | float32 | DOI-resolved 3-D hit position |
| `hit2_x_mm`, `hit2_y_mm`, `hit2_z_mm` | float32 | DOI-resolved 3-D hit position |
| `e1_keV`, `e2_keV` | float32 | deposited energies |
| `t1_ns`, `t2_ns` | float64 | timestamps — for TOF |
| `truth` | int8 | 0 = true, 1 = phantom scatter, 2 = random, −1 = unclassified |

### Companion `coincidences_<config>_meta.csv` (one wide row)

`detector_config` (e.g. `CRYSP_baseline`/`LYSO`/`BGO`), geometry params (ring Ø,
AFOV, crystal dims, σ_E, σ_t, DOI res), `energy_window_lo_keV`,
`energy_window_hi_keV`, `coinc_time_window_ns`, `source_file`, `Nj_budget_0..4`,
`realization_index`, `random_seed`, `geant4_version`.

### Invariants
- Reconstruction (Stage C) consumes **only this file** — never Geant4 truth
  (CLAUDE.md invariant 5). The `truth` column is for our diagnostics; a
  reconstructor must work with it removed.
