# File formats

The columns of every file shared between the stages: the isotope encoding, the
units, and the layout of each file. The code mirrors it (`common/Isotopes.hh` for
C++, `common/isotopes.py` for Python). If code and this file disagree, the code is
wrong. The detector stage (Stage B) is a separate downstream repo; each frozen run
carries a short copy of these formats.

Spec reference: the LaTeX docs in `latex/` — `01_user_guide` for the pipeline
overview, `04_source_reference` for the annotated, canonical form of this file.

**File format: CSV.** Files are flat, columnar **CSV** (one row per record, one
scalar per column, header row, units as `_mm`/`_keV`/`_ns` column suffixes). Each
data file has a companion `*_meta.csv` (a single wide row) of run-level metadata.
CSV was chosen over HDF5 because the files are small and it needs no extra
dependency (pandas-native, readable). **HDF5 is deferred** — if size or read speed
ever demands it, switch to HighFive (C++) / h5py (Python); the columns below are
the same, so only the container changes.

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
| `beam_energy_MeV` | float | nominal proton energy; for a SOBP run this is the fallback single energy (real spectrum in `sobp_layers.csv`) |
| `beam_sigma_mm` | float | pencil-beam Gaussian σ at entrance |
| `phantom_material` | string | Geant4 NIST name, e.g. `"G4_BRAIN_ICRP"` |
| `phantom_diameter_mm`, `phantom_length_mm` | float | phantom cylinder size (standard 160 × 160) |
| `phantom_mass_g` | float | phantom mass (for dose) |
| `edep_total_MeV` | float | total energy deposited in the phantom |
| `dose_total_Gy` | float | whole-phantom dose for `n_protons` |
| `target_dose_Gy` | float | dose in the target box for `n_protons` (normalization) |
| `target_mass_g` | float | target-box mass |
| `target_radius_mm` | float | target box (cylinder) radius |
| `target_prox_depth_mm`, `target_dist_depth_mm` | float | target box faces, depth from the entrance |
| `Np_per_Gy` | float | protons for 1 Gy (= `n_protons` / `target_dose_Gy`) |
| `geant4_version` | string | e.g. `"11.4.1"` |
| `physics_list` | string | e.g. `"QGSP_BIC_HP"` |
| `random_seed` | int | master seed of the run |

`Pj_per_proton` is derivable as per-isotope count / `n_protons`. The absolute
normalization scales production to a clinical dose via `target_dose_Gy`.

### Coordinate frame

All positions — `emitters.csv`, the phantom, the target box — share one frame:
the phantom is a cylinder **centred at the origin** with its axis along **+z**,
the beam direction. It spans `z ∈ [−L/2, +L/2]` and `r = √(x²+y²) ∈ [0, R]`
(L = `phantom_length_mm`, R = `phantom_diameter_mm`/2). The beam enters at the
`z = −L/2` face; depth from the entrance is `d = z + L/2` (so the target-box
depths map to `z = depth − L/2`). A consumer that places the phantom with its
centre at the origin is co-registered with the `emitters.csv` source.

**Medium / attenuation map.** The source is detector-independent, but **not
medium-independent**: the downstream sim must propagate the 511 keV annihilation
photons through the phantom (attenuation + scatter) and again needs μ(511 keV)
for reconstruction attenuation correction. For the standard scenario this is
fully specified by the three `run_meta.csv` scalars — `phantom_material` (a NIST
name that resolves to exact composition + density) and `phantom_diameter_mm` /
`phantom_length_mm` (a homogeneous cylinder on the +z axis, centred at origin,
sharing the `emitters.csv` frame). No voxel map is needed while the phantom is
homogeneous. **Heterogeneous phantoms are deferred** — a voxelized patient
(lung, bone, cavities) cannot be described by those scalars; supporting it means
freezing a companion voxel material/density map (the attenuation map, or the
GDML geometry) into the scenario snapshot and adding a pointer column to
`run_meta.csv`. Until then this schema is homogeneous-phantom-only by
construction. The composition + μ are materialized into File 3 (below) by
`common/phantom_material.py` and frozen into every scenario snapshot.

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

## File 2: handoff output — the budget (deterministic) + realizations (stochastic)

The handoff is split at the **A|B seam** into a deterministic, RNG-free budget
(this repo) and the stochastic Poisson realizations (the downstream detector
study). The annihilation events are **not materialized** — the handoff
writes only the per-isotope counts; the detector draws that many annihilation
points from File 1 on the fly, with seed `master_seed + realization` (so every
detector gets the identical source). Method: `latex/03_decay_kinetics.tex`.

### File 2a: `sampling_budget_<scenario>.csv` — Stage B0 deterministic budget

Written by `decay_sampling/budget.py`. The thin quantity crossing the A|B seam:
the measured decays N_j from Eq. 1, detector-independent and seed-free. One row
per isotope.

| column | type | meaning |
|--------|------|---------|
| `isotope_id` | int8 | which species (encoding above) |
| `N_expected` | float | measured decays from Eq. 1, `P_j(D)·survival_j` |

Companion `sampling_budget_<scenario>_meta.csv` (one wide row): `scenario`,
`source_file` (the `emitters.csv` drawn from), `dose_Gy`, `t_irr_s`, `t_del_s`,
`t_meas_s`, `target_dose_Gy`.

### File 2b: `sampling_realizations_<scenario>.csv` — stochastic Poisson draws

Written by `decay_sampling/budget_gen.py` (which **reads File 2a**), and destined
to move to the downstream detector study. All RNG lives here, never in File 2a. One row per
(realization, isotope).

| column | type | meaning |
|--------|------|---------|
| `realization` | int | realization index (0 … Z−1) |
| `isotope_id` | int8 | which species (encoding above) |
| `N_poisson` | int | this realization's draw `M_j ~ Poisson(N_expected)` |

Companion `sampling_realizations_<scenario>_meta.csv` (one wide row): `scenario`,
`source_budget` (the File 2a it drew from), `n_realizations` (Z), `master_seed`.

### Invariants
- N_j (the budget) is set **only** by File 2a (CLAUDE.md invariant 3); it is
  deterministic and detector-independent.
- All randomness lives in File 2b; changing Z or the seed never touches File 2a,
  so the source side stays reproducible without an RNG.
- Every detector config consumes the **identical** source — guaranteed by the
  same `emitters.csv` + budget + seed (CLAUDE.md invariant 4).

---

## File 3: phantom medium — `phantom_material_<name>.csv` (+ `_meta.csv`)

Written by `common/phantom_material.py` (and emitted into every snapshot by
`tools/snapshot_scenario.py`). The 511 keV photon-transport / attenuation-
correction data for the homogeneous phantom — the medium the source implies but
`emitters.csv` does not carry (see the "Medium / attenuation map" note above).
`<name>` is the `run_meta.csv` `phantom_material`, lower-cased and sanitized
(e.g. `g4_brain_icrp`). Compositions are the authoritative Geant4 NIST
definitions; `μ` is Compton/Klein-Nishina (coherent + photoelectric ~1–2 %, not
included).

### `phantom_material_<name>.csv` — one row per element

| column | type | meaning |
|--------|------|---------|
| `element` | string | element symbol (H, C, N, O, …) |
| `Z` | int | atomic number |
| `A_g_mol` | float | standard atomic weight |
| `mass_fraction` | float | fraction by mass (rows sum to 1) |

### `phantom_material_<name>_meta.csv` — one wide row

| column | type | meaning |
|--------|------|---------|
| `material` | string | Geant4 NIST name (e.g. `G4_BRAIN_ICRP`) |
| `energy_keV` | float | photon energy the μ is for (511, annihilation) |
| `density_g_cm3` | float | mass density ρ |
| `mean_excitation_eV` | float | Geant4 Imean (diagnostic; not used for μ) |
| `mu_rho_cm2_g` | float | mass attenuation coefficient μ/ρ |
| `mu_cm_inv` | float | linear attenuation coefficient μ |
| `mu_mm_inv` | float | μ in mm⁻¹ |
| `mean_free_path_cm` | float | 1/μ |
| `note` | string | short material description |

---

## Downstream: `coincidences_<config>.csv`

The detector output (the coincidence list) is produced and documented by the
downstream detector simulation, which reads File 1 + the budget. Its columns are
defined in that repo, not here.
