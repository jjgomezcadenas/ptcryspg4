# CLAUDE.md — Proton-therapy PET source generation (Stage A + handoff)

Orientation for any Claude Code session on this repo.

## Purpose

Produce the **positron-emitter source** that a PET detector simulation needs for
**in-room proton-therapy range verification**. A Geant4 proton run (Stage A) sends
protons through a phantom and records where β⁺ emitters are produced and where the
positrons annihilate; a short Python handoff (Stage B0) turns the production into
the number of decays a scanner would measure.

The output — `emitters.csv` (annihilation points) + `run_meta.csv` (normalization)
+ the per-isotope budget — is **detector-independent**: it knows nothing about any
detector, ring, or PET design. A separate downstream simulation reads it and models
a detector. This repo stops at the source.

## The spec

`docs/simulate_pt_pet.tex` describes the physics and the pipeline. Read it first.
If code and spec disagree, the spec wins unless a decision here supersedes it.
This file records implementation decisions; the spec records the method.

## Architecture — stages joined by CSV files

```
[A]  Geant4 transport    protons → phantom → β+ emitter → annihilation     RUNS ONCE
         ⟹  emitters.csv + run_meta.csv   (detector-independent source)
[B0] Time-decay handoff (Python)   P_j → N_j(t_del)
         ⟹  sampling_budget_<scenario>.csv   (measured decays per isotope)
       ─────── frozen as a named scenario in the ptcrysp-scenarios data repo ───────
       a downstream PET detector simulation reads the scenario (separate repo)
```

Stage A runs once; its output is frozen as a named scenario in the
`ptcrysp-scenarios` data repo (`tools/snapshot_scenario.py`). The detector study is
a separate repo and is not described here.

## Invariants — do not violate these

1. **Geant4 owns space; the Python handoff owns number and acquisition timing.**
   Let emitters decay through Geant4's standard radioactive-decay process at their
   natural lifetime — do **not** override lifetimes or force prompt decay. (Forcing
   prompt is an unnecessary intervention: the recoil nucleus moves only ~µm before
   decaying, so PROD/ANH are identical either way, and Geant4 still yields exactly
   **one annihilation per produced emitter**, production-proportional, 1:1.) The
   measured number N_j and the acquisition timing are set solely by the handoff
   (Stage B0). The positron's global-time stamp is physically real but carries no
   acquisition meaning and is never written — `emitters.csv` has no time column.
2. **Stage A runs once.** Never re-run proton transport when the detector changes.
3. **The measured count N_j is set only at the handoff.** Raw produced counts are
   *production*-proportional; the measured ¹⁵O/¹¹C mix comes from the per-isotope
   budget N_j (`budget.py`, deterministic), not from raw counts. The Poisson
   realizations of N_j are drawn later, by the downstream study.
4. **The source is detector-independent.** Every downstream consumer reads the same
   `emitters.csv` + budget; nothing detector-specific is baked into this repo, so
   the same source can feed any PET detector simulation.

## Tech stack (decided)

- **Stage A:** C++ / **Geant4** (11.4.x; CMake) — proton transport only.
- **Stage B0 (handoff) + analysis:** **Python** (`numpy`, `scipy`, `pandas`).
- **File format:** **CSV** (flat, columnar), each with a companion `*_meta.csv` for
  run-level metadata. Chosen over HDF5: files are small, no extra dependency,
  pandas-native and readable. **HDF5 is deferred** — revisit (HighFive / h5py) only
  if size or read speed demands it; the columns are the same, so the switch is cheap.
- **Stage A starting point:** we built a minimal custom app, `stageA_transport/`,
  rather than stripping the Geant4 `hadrontherapy` example. Radioactive decay is
  active in `QGSP_BIC_HP` — no prompt-decay override.

The detector (Stage B) and reconstruction (Stage C) are a separate downstream
simulation, not part of this repo.

## Units convention

Positions **mm**, energies **keV** (beam energy **MeV**), times **ns**, dose
**Gy**. State units as column-name suffixes (e.g. `prod_x_mm`) and in the
companion `*_meta.csv`.

## File formats

`common/SCHEMA.md` has the full column list; this is the summary.

### `emitters.csv` — Stage A output (one row per β⁺ emitter)
Flat columns:

| column | type | meaning |
|--------|------|---------|
| `event_id` | int | primary that produced it (diagnostic) |
| `isotope_id` | int8 | 0=¹⁵O 1=¹¹C 2=¹³N 3=¹⁰C 4=¹⁴O |
| `prod_x_mm,prod_y_mm,prod_z_mm` | float | production point (mm) — the **truth** map |
| `anh_x_mm,anh_y_mm,anh_z_mm` | float | annihilation point (mm) — detector source |

> One row holds both points (positron range = `anh − prod`). Run-level metadata
> (`n_protons`, `beam_energy_MeV`, `dose_total_Gy`, `geant4_version`,
> `physics_list`, `random_seed`, …) is in the companion `run_meta.csv`;
> `Pj_per_proton` is derivable as per-isotope count / `n_protons`. Stage A also
> writes `depth_dose.csv` (the Bragg profile).

### `sampling_budget_<scenario>.csv` — handoff output (one row per isotope)
`isotope_id`, `N_expected` (measured decays from Eq. 1, for 1 Gy at the scenario's
timing). Companion `*_meta.csv` carries the dose and timing (`t_irr`, `t_del`,
`t_meas`). See `common/SCHEMA.md` File 2.

## Fixed parameters

**Standard scenario — Parodi skull-base reference.** The reference case is
Parodi et al. 2008 (`docs/proton-therapy-ptet-doses-mc.pdf`, Fig. 2a / Table 2):
a proton **SOBP** field delivering **1 Gy** to a target inside a head, so our
absolute β⁺ yields are **checkable against their published numbers**. Everything
else (PMMA, water, single pencil, 2 Gy, spine, in-room timing) is a **variant**.

| element | standard value | source |
|---|---|---|
| site / beam | skull-base, proton **SOBP + lateral fluence** | Fig. 2a |
| phantom | homogeneous **brain** (`G4_BRAIN_ICRP`), cylinder **Ø16 cm × 16 cm** | head ≈ cylinder |
| target | box **Ø6 cm × 5 cm at ~8 cm depth** | Fig. 2a dose box |
| dose | **1 Gy** to the target | Fig. 2 caption |
| timing | in-room **t_irr=60, t_del=120, t_meas=1200 s** (cyclotron); t_del swept | — |

**Cross-checks (proton, head, 1 Gy):**
- *Production* — our G4 integral yields vs **Table 2**: ¹⁵O 9.4e7, ¹¹C 7.7e7,
  ¹³N 9.4e6, ¹⁰C 1.6e6, ¹⁴O 5.5e5 → total **~1.8e8 /Gy**.
- *Measured, no washout* — handoff Eq. 1 at the **offline** point (t_del=300,
  t_meas=1800) vs **Table 4/5**: total **~6.1e7** (the offline variant; the
  in-room baseline keeps more ¹⁵O).

**Beam (gun):** standard = the SOBP above. *Variant:* single mono-energetic
pencil (100 MeV, σ≈3 mm) — a clean range-verification testbed, but it paints a
thin track, not a target volume, so its integral yield is unrealistically small
→ use it for shape tests, not for absolute rates.

**Phantom material** is macro-selectable (`/stageA/phantom/material`, before
`/run/initialize`): **`G4_BRAIN_ICRP`** (standard) · `G4_TISSUE_SOFT_ICRP`
(oxygen-rich tissue) · `G4_PLEXIGLASS` (PMMA, carbon-rich benchmark) · `G4_WATER`
(pure-¹⁵O extreme). The ¹⁵O/¹¹C mix is computed from MC P_j, never assumed
(O15/C11 ≈ 1.2 Parodi head, ~1.7 our tissue, ~0.6 PMMA).

**Isotopes / half-lives:** ¹⁵O 122 s, ¹¹C 1223 s, ¹³N 598 s, ¹⁰C 19.3 s,
¹⁴O 70.6 s (100 % β⁺ assumed in the handoff).

**Acquisition / handoff (cyclotron, post-beam).** Continuous/cyclotron limit
(post-beam → pulsed structure irrelevant), so the measured decays per species
are the **three-factor expression — the *only* handoff equation:**

```
N_j = P_j · (1−e^(−λ_j·t_irr))/(λ_j·t_irr) · e^(−λ_j·t_del) · (1−e^(−λ_j·t_meas))
           └──── build-up ────┘              └ transport ┘   └──── window ────┘
```

Baseline timing (in-room): **t_irr=60, t_del=120, t_meas=1200 s** — `t_del` is
the key swept parameter (120 s is the fastest realistic in-room delay; longer
delay = less ¹⁵O). *Conservative/offline variant:* t_del=300, t_meas=1800 s
(where the Parodi Table 4 measured-decay check applies). The positron-range floor
(¹⁵O longest) matters at any of these.

**Full handoff method in `docs/handoff.tex`** — the time-decay model and the
absolute normalization `P_j(D)=count_j·D/target_dose`. The measured budget N_j is
computed by `decay_sampling/budget.py`.

## Reference material

- `docs/ptet_mc.txt` — **Parodi, Bortfeld & Haberer 2008** — the time-decay model
  (Eqs. 1–7) and integral isotope yields. *Core.*
- `docs/ptet_doses_mc.txt` (+ `proton-therapy-ptet-doses-mc.pdf`) — Parodi
  production/dose MC methodology and cross-section caveats → Stage A.
- `docs/proton_therapy.txt` (+ `Proton_Therapy_Verification_with_PET_Imaging.pdf`)
  — Zhu & El Fakhri review → modality / isotope / in-room rationale.
- `biblio.bib` keys: `parodi2008comparison`, `parodi2007patient`, `zhu2013proton`.
- External, fetch as needed: Geant4 `hadrontherapy` example README + Physics
  Reference (radioactive-decay / prompt-decay biasing); Levin & Hoffman 1999
  (positron-range validation).

## Suggested repo layout

```
stageA_transport/    # Geant4 C++: protons → emitters.csv + run_meta.csv
field_design/        # Python: SOBP beam design (sobp.py) + depth-dose plots
decay_sampling/      # Python: time-decay budget (budget.py) + realizations (budget_gen.py)
analysis_transport/  # Python: validate Stage A output (dashboard, diagnostics)
tools/               # snapshot_scenario.py: freeze a run into the scenarios repo
common/              # shared schema, units, isotope table
docs/                # spec (simulate_pt_pet.tex), sobp.tex, handoff.tex, refs
data/                # generated CSV (gitignored)
```

Frozen Stage-A runs live in the `ptcrysp-scenarios` data repo, one named directory
per scenario; the downstream detector simulation reads from there.

## Beam: SOBP field (`field_design/`, method in `docs/sobp.tex`)

The standard scenario uses a **Spread-Out Bragg Peak**, not a single pencil. The
depth field is implemented and verified:
- `field_design/sobp.py` computes the energy-layer weights (Bortfeld range-energy
  + Abel-inversion flattening weights + an `exp(µR)` attenuation correction,
  `--mu` tuned) → `data/sobp_layers.csv`.
- The gun (`BeamConfig` + `/stageA/beam/layers <file>`) samples a layer energy
  per primary (stochastic; per-layer Poisson noise ≲0.2% at ≥10⁷ protons).
- `field_design/plot_sobp.py` renders the realized G4 depth-dose
  (`data/depth_dose.csv`) → `data/sobp_g4.png`. Verified flat plateau (~7%) over
  the target depth with a sharp distal edge.

## Build / run

**Toolchain (this machine).** Geant4 **11.4.1** at `~/Software/geant4/install`,
built multithreaded (tasking) with **Qt visualization ON** (ToolsSG-Qt drivers:
`TSG_QT_ZB` software / `TSG_QT_GLES` GPU; no classic OpenGL). No ROOT. Activate
with `source $HOME/Software/geant4/install/bin/geant4.sh` (in `~/.zshrc`). Output
is CSV — no HDF5/HighFive dependency.

### Stage A — `proton_transport`

```bash
source $HOME/Software/geant4/install/bin/geant4.sh
cd stageA_transport && mkdir -p build && cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build . -j 18
```

Run:

```bash
./proton_transport            # interactive Qt viewer; runs vis.mac,
                              #   draws the phantom + shoots 1 proton (more: /run/beamOn N)
./proton_transport sobp.mac   # batch (18 threads); writes emitters/run_meta/depth_dose.csv to data/
```

SOBP run (depth field): design the layers, then run the macro that loads them:

```bash
python3 field_design/sobp.py --mu 0.025          # -> data/sobp_layers.csv
./proton_transport sobp.mac                      # macro loads the layer table
python3 field_design/plot_sobp.py                # -> data/sobp_g4.png
```

`build/compile_commands.json` (the `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON` flag)
feeds VS Code IntelliSense; `.vscode/settings.json` points the C++ extension at
it via CMake Tools.

### Analysis + handoff (Python)

System `python3` (already has numpy/scipy/pandas/matplotlib — no venv). Validate
a Stage-A run:

```bash
python3 analysis_transport/validate_transport.py   # -> data/transport_validation.png
```

Deps listed in `analysis_transport/requirements.txt` and
`decay_sampling/requirements.txt` (numpy, scipy, pandas).

### Handoff + snapshot

```bash
python3 decay_sampling/budget.py                 # -> data/sampling_budget_inroom.csv (measured N_j)
python3 decay_sampling/budget_gen.py             # -> data/sampling_realizations_inroom.csv (Poisson draws)
python3 tools/snapshot_scenario.py head_sobp_1e7 # freeze data/ into ../ptcrysp-scenarios
```

`budget.py` writes the per-isotope measured count N_j (no random numbers).
`budget_gen.py` draws the Poisson realizations of N_j — this step will move to the
downstream study. `snapshot_scenario.py` copies a finished run into the
`ptcrysp-scenarios` repo as a named scenario.

## First-session checklist

**Status:** Stage A is done and validated (custom `stageA_transport/`, MT on 18
threads; `emitters.csv` + `run_meta.csv` + `depth_dose.csv`; yields sane, Bragg
curve, endpoint-ordered positron ranges). The handoff is done (`budget.py` +
`budget_gen.py`). The standard run (`head_sobp_1e7`, 10⁷ protons) is frozen in the
`ptcrysp-scenarios` repo. The source side is complete; the detector study is a
separate downstream repo.

1. Read `docs/simulate_pt_pet.tex` end to end.
2. Stand up the Stage-A app; confirm proton range and dose in the phantom.
3. Write `emitters.csv` (+ `run_meta.csv`); rely on standard radioactive decay
   (no prompt-decay override); check yields-per-proton against the literature.
4. Compute the A→B handoff (`budget.py`, spec §3) and check the ¹⁵O/¹¹C mix vs t_del.
5. Freeze a run into `ptcrysp-scenarios` with `tools/snapshot_scenario.py`.
