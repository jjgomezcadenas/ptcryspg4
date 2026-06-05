# CLAUDE.md — Proton-therapy PET detector-comparison simulation

> Copy this file to the **root of the new coding repository**. It is the orientation
> document for any Claude Code session working on this project.

## Purpose

Compare candidate PET detectors (the CRYSP family) for **in-room** proton-therapy
range verification. The chain has two simulation stages: a **Geant4** proton run
(Stage A, this repo) that produces the positron-emitter source, and an **analytic
detector Monte Carlo in Julia** (`PTCryspMC.jl`, a separate repo) that turns that
source into a **coincidence list**. Image reconstruction is a separate stage. The
goal is *not* to reproduce any clinical study, but to rank detectors on how well
they recover the proton range at realistic, photon-starved statistics.

## The spec

`docs/simulate_pt_pet.tex` describes the physics and the pipeline. Read it first.
If code and spec disagree, the spec wins unless a decision here supersedes it.
This file records implementation decisions; the spec records the method.

## Architecture — stages joined by CSV files

```
ptcryspg4 (this repo) — Geant4 + Python, RUNS ONCE
  [A]  Geant4 transport    protons → phantom → β+ emitter → annihilation
          ⟹  emitters.csv + run_meta.csv   (detector-independent source)
  [B0] Time-decay handoff (Python)   P_j → N_j(t_del)
          ⟹  sampling_budget_<scenario>.csv   (measured decays per isotope)
        ─────── frozen and passed via the ptcrysp-scenarios data repo ───────
PTCryspMC.jl (separate repo) — Julia + Python, RUNS PER DETECTOR
  [B]  Analytic detector MC   annihilation events → detector → coincidences
          ⟹  coincidences_<config>.csv
  [C]  Reconstruction         coincidence list → image → σ(range)   DEFERRED
```

Stage A runs once; its output is frozen as a named scenario in the
`ptcrysp-scenarios` data repo (`tools/snapshot_scenario.py`). `PTCryspMC.jl` reads
a scenario from there and runs the per-detector study.

## Invariants — do not violate these

1. **Geant4 owns space; the analytic bookkeeping owns number and acquisition
   timing.** Let emitters decay through Geant4's standard radioactive-decay
   process at their natural lifetime — do **not** override lifetimes or force
   prompt decay. (Forcing prompt is an unnecessary intervention: the recoil
   nucleus moves only ~µm before decaying, so PROD/ANH are identical either way,
   and Geant4 still yields exactly **one annihilation per produced emitter**,
   production-proportional, 1:1.) The measured number N_j and the acquisition
   timing are owned solely by the analytic bookkeeping (Stage B0). The positron's
   global-time stamp is physically real but carries no acquisition meaning and is
   never written — `emitters.csv` has no time column.
2. **Stage A runs once.** Never re-run proton transport when the detector changes.
3. **The measured count N_j is set only at the handoff.** Raw produced counts are
   *production*-proportional; the measured ¹⁵O/¹¹C mix comes from the per-isotope
   budget N_j (`budget.py`, deterministic), not from raw counts. The Poisson
   realizations of N_j are drawn later, in `PTCryspMC.jl`.
4. **Every detector config consumes the identical source** (same `emitters.csv`,
   same sampled-event set per realization) — otherwise the ranking is polluted by
   per-run statistical drift.
5. **Reconstruction uses only the coincidence list** — never the Geant4 truth or
   the detector's internal state.

## Tech stack (decided)

- **Stage A:** C++ / **Geant4** (11.4.x; CMake) — proton transport only.
- **Stage B0 (handoff) + analysis:** **Python** (`numpy`, `scipy`, `pandas`).
- **Stage B (detector) + C (reconstruction):** **Julia**, in the separate repo
  `PTCryspMC.jl`. The detector is an **analytic γ-transport Monte Carlo** (adapted
  from LXeMC), not Geant4: 511 keV photons are tracked through the phantom (treated
  as water) and the crystal by Klein–Nishina Compton + photoelectric, and the ring
  acceptance is a ray–cylinder test. Geant4 was dropped here because the geometry
  is simple, the same physics is easy to do analytically, and it avoids a
  per-detector Geant4 build. Reconstruction is **deferred** (list-mode MLEM/OSEM,
  or CASToR, or STIR; must be list-mode, DOI-aware, optionally TOF).
- **File format:** **CSV** (flat, columnar), each with a companion `*_meta.csv` for
  run-level metadata. Chosen over HDF5: files are small, no extra dependency,
  pandas-native and readable. **HDF5 is deferred** — revisit (HighFive / h5py) only
  if size or read speed demands it; the columns are the same, so the switch is cheap.
- **Stage A starting point:** we built a minimal custom app, `stageA_transport/`,
  rather than stripping the Geant4 `hadrontherapy` example. Radioactive decay is
  active in `QGSP_BIC_HP` — no prompt-decay override.

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

### `coincidences_<config>.csv` — Stage B output (PTCryspMC.jl; one row per accepted coincidence)
Flat columns:

| field      | type       | meaning                          |
|------------|------------|----------------------------------|
| `hit1_xyz` | float32[3] | DOI-resolved 3-D position (mm)   |
| `hit2_xyz` | float32[3] | DOI-resolved 3-D position (mm)   |
| `e1`,`e2`  | float32    | deposited energies (keV)         |
| `t1`,`t2`  | float64    | timestamps (ns) — for TOF        |
| `truth`    | int8       | optional: 0 true,1 scatter,2 rand|

Companion `*_meta.csv`: `detector_config`, geometry params, `energy_window_keV`,
`coinc_time_window_ns`, `Nj_budget` (per isotope), `realization_index`,
`random_seed`. (Hit positions may carry `_mm` / energies `_keV` / times `_ns`
column suffixes.)

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
→ use it for shape/detector tests, not for absolute rates.

**Phantom material** is macro-selectable (`/stageA/phantom/material`, before
`/run/initialize`): **`G4_BRAIN_ICRP`** (standard) · `G4_TISSUE_SOFT_ICRP`
(oxygen-rich tissue) · `G4_PLEXIGLASS` (PMMA, carbon-rich benchmark) · `G4_WATER`
(pure-¹⁵O extreme). The ¹⁵O/¹¹C mix is computed from MC P_j, never assumed
(O15/C11 ≈ 1.2 Parodi head, ~1.7 our tissue, ~0.6 PMMA).

**Isotopes / half-lives:** ¹⁵O 122 s, ¹¹C 1223 s, ¹³N 598 s, ¹⁰C 19.3 s,
¹⁴O 70.6 s (100 % β⁺ assumed in the bookkeeping).

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

**Full handoff method in `docs/handoff.tex`** — the time-decay model, absolute
normalization `P_j(D)=count_j·D/target_dose`, and the σ(range) figure of merit.
The measured budget N_j is computed here (`decay_sampling/budget.py`); the Poisson
realizations and σ(range) are computed in `PTCryspMC.jl`.

**CRYSP baseline detector** (from `crysp_for_ht.tex` / Soleti 2024): ring Ø
77.4 cm, AFOV 102.4 cm, monolithic crystals 48×48×37 mm, **6.3 % FWHM** energy
resolution @ 511 keV, **1.7 mm** 3-D resolution incl. DOI, NEMA sensitivity
~120 kcps/MBq, ~1 µs decay, **no TOF**, **no intrinsic radioactivity**.

**Comparators / variants:** LYSO (~10 % E-res, **¹⁷⁶Lu** intrinsic background —
inject as volumetric source), BGO (poor E-res, no intrinsic activity, high
density); CRYSP variants: room-temperature CsI(Tl), cryogenic BGO,
BGO-core/CsI-wing hybrid.

**Discriminators to sweep:** AFOV/solid-angle (acceptance ≈ L/√(L²+R²)), energy
resolution (scatter rejection on oblique LORs), DOI (parallax), intrinsic
radioactivity (LYSO only), TOF (treated as low-value — show diminishing returns).

**Figure of merit:** σ(range) in mm — spread of the fitted distal fall-off over M
Poisson realizations — vs. counts (or dose, or t_del).

## Reference material (in the companion paper repo)

- `docs/ptet_mc.txt` — **Parodi, Bortfeld & Haberer 2008** — the time-decay model
  (Eqs. 1–7) and integral isotope yields. *Core.*
- `docs/ptet_doses_mc.txt` (+ `proton-therapy-ptet-doses-mc.pdf`) — Parodi
  production/dose MC methodology and cross-section caveats → Stage A.
- `docs/proton_therapy.txt` (+ `Proton_Therapy_Verification_with_PET_Imaging.pdf`)
  — Zhu & El Fakhri review → modality / isotope / in-room rationale.
- `crysp_for_ht.tex` — the CRYSP design and detector requirements (the *why*
  behind the comparison).
- `biblio.bib` keys: `parodi2008comparison`, `parodi2007patient`,
  `zhu2013proton`, `Soleti:2024flw`, `mikhailik2015luminescence`,
  `moszynski2005energy`.
- External, fetch as needed: Soleti et al. 2024 (CRYSP detector numbers);
  Geant4 `hadrontherapy` example README + Physics Reference (radioactive-decay /
  prompt-decay biasing); Levin & Hoffman 1999 (positron-range validation);
  CASToR (Merlin 2018) / STIR if used for Stage C.

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

The detector and reconstruction (Stages B, C) live in the separate `PTCryspMC.jl`
repo. Frozen Stage-A runs live in the `ptcrysp-scenarios` data repo, one named
directory per scenario; `PTCryspMC.jl` reads from there.

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
- **Pending:** lateral fluence over the target cross-section (→ a dose *box*),
  dose normalization to 1 Gy in the target, and the Parodi Table 2 cross-check.
  Until then the gun is laterally a σ=3 mm pencil.

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
./proton_transport            # interactive Qt viewer; runs macros/vis.mac,
                              #   draws the phantom + shoots 1 proton (more: /run/beamOn N)
./proton_transport run.mac    # batch (18 threads); writes emitters/run_meta/depth_dose.csv to data/
```

SOBP run (depth field): design the layers, then load them in the macro before
`/run/beamOn` with `/stageA/beam/layers ../../data/sobp_layers.csv`:

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
`budget_gen.py` draws the Poisson realizations of N_j — this step will move to
`PTCryspMC.jl`. `snapshot_scenario.py` copies a finished run into the
`ptcrysp-scenarios` repo as a named scenario.

## First-session checklist

**Status:** Stage A is done and validated (custom `stageA_transport/`, MT on 18
threads; `emitters.csv` + `run_meta.csv` + `depth_dose.csv`; yields sane, Bragg
curve, endpoint-ordered positron ranges). The handoff is done (`budget.py` +
`budget_gen.py`). The standard run (`head_sobp_1e7`, 10⁷ protons) is frozen in the
`ptcrysp-scenarios` repo. **Next: the detector study in `PTCryspMC.jl`.**

1. Read `docs/simulate_pt_pet.tex` end to end.
2. Stand up the Stage-A app; confirm proton range and dose in the phantom.
3. Write `emitters.csv` (+ `run_meta.csv`); rely on standard radioactive decay
   (no prompt-decay override); check yields-per-proton against the literature.
4. Compute the A→B handoff (`budget.py`, spec §3) and check the ¹⁵O/¹¹C mix vs t_del.
5. Build the analytic detector (`PTCryspMC.jl`): annihilation events → detector
   response → coincidence list → σ(range).
