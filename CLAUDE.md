# CLAUDE.md — Proton-therapy PET source generation (Stage A + handoff)

Orientation for any Claude Code session on this repo.

## Working style

Do **not** present "shopping lists" — the rigid, multiple-choice
question-with-fixed-options format (the `AskUserQuestion` tool). The user finds
them frozen and impersonal. Instead, when a decision needs the user's input,
talk it through in prose: lay out the options and trade-offs conversationally,
give a recommendation, and ask one question at a time. Options are welcome;
canned option menus are not.

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

## The documentation (`latex/`)

Our LaTeX documentation lives in `latex/`, numbered in reading order; `docs/`
holds only reference papers (`.pdf`, `.txt`). Build any of them with `pdflatex`.

- `latex/01_user_guide.tex` — the user guide: the physics (induced activity and
  its time evolution, the decay kinetics that set the measured counts, what the
  generated source looks like, with figures) and how to run the pipeline. **Read
  it first.**
- `latex/02_beam_design.tex` — how the SOBP depth field is designed
  (Bortfeld + Abel-inversion weights + attenuation correction), as built.
- `latex/03_decay_kinetics.tex` — the decay-plus-timing model that turns produced
  emitters into measured decays, and the σ(range) figure of merit.
- `latex/04_source_reference.tex` — the downstream-consumer interface contract:
  the scenario data product (files, columns, metadata, coordinate frame) and the
  recipe for turning it into a PET acquisition, with the head/brain case worked
  through. The canonical, annotated form of the per-snapshot `SCHEMA.md`.

> **Terminology — "decay kinetics" = "the handoff".** The decay-kinetics step
> (doc `03_decay_kinetics`, code `decay_sampling/`) is exactly what the
> architecture diagram below and the rest of this file call **the handoff**
> (Stage B0): the radioactive-decay + acquisition-timing model that turns
> produced emitters `P_j` into measured decays `N_j`. The two names are
> interchangeable; the doc was renamed for clarity, the concept kept its name.

This file (CLAUDE.md) records the implementation decisions, parameters, and
build/run details a coding session needs, and refers to the guide for the physics
rather than repeating it.

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

**Acquisition / handoff (cyclotron, post-beam).** The measured decays per species
are the three-factor expression:

```
N_j = P_j · (1−e^(−λ_j·t_irr))/(λ_j·t_irr) · e^(−λ_j·t_del) · (1−e^(−λ_j·t_meas))
           └──── build-up ────┘              └ transport ┘   └──── window ────┘
```

Baseline timing (in-room): **t_irr=60, t_del=120, t_meas=1200 s** — `t_del` is
the key swept parameter (120 s is the fastest realistic in-room delay; longer
delay = less ¹⁵O). *Conservative/offline variant:* t_del=300, t_meas=1800 s
(where the Parodi Table 4 measured-decay check applies). The positron-range floor
(¹⁵O longest) matters at any of these.

The factors are explained in the user guide (§decay kinetics); the full
derivation, including pulsed deliveries, is in `latex/03_decay_kinetics.tex`. Absolute
normalization `P_j(D)=count_j·D/target_dose`; the budget N_j is computed by
`decay_sampling/budget.py`.

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
analysis_transport/  # Python: validate Stage A (validate_transport.py) + make_figures.py
tools/               # snapshot_scenario.py: freeze a run into the scenarios repo
common/              # shared schema, units, isotope table
latex/               # our LaTeX docs (01_user_guide … 04_source_reference) + figures + biblio
docs/                # reference papers only (.pdf, .txt)
data/runs/<run_tag>/ # generated CSV, one self-contained dir per run (gitignored)
```

Frozen Stage-A runs live in the `ptcrysp-scenarios` data repo, one named directory
per scenario; the downstream detector simulation reads from there.

## Beam: SOBP field (`field_design/`, method in `latex/02_beam_design.tex`)

The standard scenario uses a **Spread-Out Bragg Peak**, not a single pencil. The
depth field is implemented and verified:
- `field_design/sobp.py` computes the energy-layer weights (Bortfeld range-energy
  + Abel-inversion flattening weights + an `exp(µR)` attenuation correction,
  `--mu` tuned) → `data/sobp_layers.csv` (the design staging input).
- The gun (`BeamConfig` + `/stageA/beam/layers <file>`) samples a layer energy
  per primary (stochastic; per-layer Poisson noise ≲0.2% at ≥10⁷ protons); the
  run copies the layer table into its own dir (`data/runs/<tag>/sobp_layers.csv`).
- `make_figures.py` (via `field_design/plot_sobp.py`) renders the realized G4
  depth-dose (`data/runs/<tag>/depth_dose.csv`) → the run's `figures/sobp_g4.png`.
  Verified flat plateau (~7%) over the target depth with a sharp distal edge.

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
./proton_transport sobp.mac   # batch (18 threads)
```

**Each run owns a directory.** Macros set the base `/stageA/output/dir
../../data/runs`; Stage A appends an auto-derived tag `<geometry>_<beam>_<N>`
(e.g. `cylinder_sobp_1e7`, `mird_head_pencil_1e5`) and writes
emitters/run_meta/phantom_regions/depth_dose there. The tag is built from the run
itself (geometry from the detector, beam = sobp if a layer table is loaded else
pencil, N the proton count), so it cannot disagree with what ran, and distinct
cases coexist under `data/runs/` without clobbering. A re-run of the same config
overwrites only its own directory.

SOBP run (depth field): design the layers, then run the macro that loads them:

```bash
python3 field_design/sobp.py --mu 0.025          # -> data/sobp_layers.csv
./proton_transport sobp.mac                      # -> data/runs/cylinder_sobp_1e7/
```

`build/compile_commands.json` (the `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON` flag)
feeds VS Code IntelliSense; `.vscode/settings.json` points the C++ extension at
it via CMake Tools.

### Analysis + handoff (Python)

System `python3` (already has numpy/scipy/pandas/matplotlib — no venv). The
analysis/handoff scripts all take a run directory as their first argument
(default `data/`); pass `data/runs/<run_tag>` to operate on one run.

**Validate a run before trusting it.** `analysis_transport/check_run.py <run_dir>`
recomputes the geometry/normalization invariants in plain Python from
`run_meta.csv` + `phantom_regions.csv` and exits non-zero on any violation (target
box in air, wrong target-mass density, target dose < whole-phantom dose, …) — the
gate that catches geometry bugs the units don't. `make_figures.py` also draws
`phantom.png` (`plot_phantom.py`): the medium regions + target box + beam in two
cross-sections, so a mis-placed box is obvious by eye. A C++ guard in `RunAction`
warns at run time if the target-box centre is in air.

### Handoff + figures + snapshot

```bash
RUN=data/runs/cylinder_sobp_1e7
python3 decay_sampling/budget.py     $RUN              # -> $RUN/sampling_budget_inroom.csv (measured N_j)
python3 decay_sampling/budget.py     $RUN --scenario fast --t-del 60
python3 decay_sampling/budget_gen.py $RUN              # -> $RUN/sampling_realizations_inroom.csv (Poisson draws)
python3 analysis_transport/make_figures.py $RUN        # -> $RUN/figures/ (geometry-aware control plots)
python3 tools/snapshot_scenario.py cylinder_sobp_1e7   # freeze that run into ../ptcrysp-scenarios
```

`budget.py` writes the per-isotope measured count N_j (no random numbers).
`budget_gen.py` draws the Poisson realizations of N_j — this step will move to the
downstream study. `make_figures.py` reads `run_meta.geometry` and writes the
figures that fit the case (cylinder → dashboard + Bragg/SOBP plateau; head →
dashboard + phantom/beam/emitter plot), all into `<run_dir>/figures/`.
`snapshot_scenario.py <run_tag>` freezes *that specific run* (identity, the files
it actually has, and its figures) into the `ptcrysp-scenarios` repo as a named
scenario (`--name` to override; defaults to the run tag). Python deps are listed
in `analysis_transport/requirements.txt` and `decay_sampling/requirements.txt`
(numpy, scipy, pandas).

## First-session checklist

**Status:** Stage A is done and validated (custom `stageA_transport/`, MT on 18
threads; `emitters.csv` + `run_meta.csv` + `depth_dose.csv`; yields sane, Bragg
curve, endpoint-ordered positron ranges). The handoff is done (`budget.py` +
`budget_gen.py`). The standard run (`head_sobp_1e7`, 10⁷ protons) is frozen in the
`ptcrysp-scenarios` repo. The source side is complete; the detector study is a
separate downstream repo.

1. Read the user guide `latex/01_user_guide.tex` end to end.
2. Stand up the Stage-A app; confirm proton range and dose in the phantom.
3. Write `emitters.csv` (+ `run_meta.csv`); rely on standard radioactive decay
   (no prompt-decay override); check yields-per-proton against the literature.
4. Compute the A→B handoff (`budget.py`; see the user guide) and check the
   ¹⁵O/¹¹C mix vs t_del.
5. Freeze a run into `ptcrysp-scenarios` with `tools/snapshot_scenario.py`.
