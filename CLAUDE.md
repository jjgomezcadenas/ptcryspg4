# CLAUDE.md — Proton-therapy PET detector-comparison simulation

> Copy this file to the **root of the new coding repository**. It is the orientation
> document for any Claude Code session working on this project.

## Purpose

Build a **Geant4** simulation that compares candidate PET detectors (the CRYSP
family) for **in-room** proton-therapy range verification. The deliverable of the
simulation is a **coincidence list**; image reconstruction is a separate,
decoupled stage. The goal is *not* to reproduce any clinical study, but to rank
detectors on how well they recover the proton range at realistic, photon-starved
statistics.

## Authoritative spec

`docs/simulate_pt_pet.tex` is the **source of truth** for the physics and the
pipeline. Read it first. If code and spec disagree, the spec wins unless a
decision here supersedes it. This file records implementation decisions; the spec
records the method.

## Architecture — three stages joined by HDF5 files

```
[A] Geant4 transport      protons → PMMA → β+ emitter → annihilation    RUNS ONCE
        ⟹  prod_anh.h5   (detector-independent source)
[B0] Time-decay bookkeeping (Python, analytic)   P_j → N_j(t_del)
        ⟹  sampled annihilation events (draw N_j per isotope)
[B] Geant4 detector       annihilation events → detector → coincidences  RUNS PER DETECTOR
        ⟹  coincidences_<config>.h5
        ───────── independent boundary ─────────
[C] Reconstruction        coincidence list → image → σ(range)            TOOLKIT DEFERRED
```

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
   never written — `prod_anh.h5` stores no time field.
2. **Stage A runs once.** Never re-run proton transport when the detector changes.
3. **N_j sampling happens only at the A→B handoff.** Raw produced counts are
   *production*-proportional; the measured ¹⁵O/¹¹C mix is set only by drawing
   `N_j` per isotope from the source file. Do not use raw counts as the budget.
4. **Every detector config consumes the identical source** (same `prod_anh.h5`,
   same sampled-event set per realization) — otherwise the ranking is polluted by
   per-run statistical drift.
5. **Reconstruction consumes only the coincidence list.** Keep it decoupled; no
   reconstruction code may depend on Geant4 truth.

## Tech stack (decided)

- **Stage A & B:** C++ / **Geant4** (target a recent release, e.g. 11.2.x; CMake).
- **File format:** **HDF5** for both interface files. In C++ use
  [HighFive](https://github.com/BlueBrain/HighFive) (header-only wrapper). In
  Python use `h5py`.
- **Stage B0 (bookkeeping/sampling) + analysis:** **Python** (`numpy`, `scipy`,
  `h5py`).
- **Reconstruction (Stage C):** **deferred.** Candidates: custom list-mode
  MLEM/OSEM, or CASToR, or STIR. Must be list-mode, DOI-aware, optionally TOF.
- **Stage A starting point:** the Geant4 advanced example **`hadrontherapy`**
  (proton beam + phantom + dose scoring + isotope production already wired);
  add the prod/anh writer. Radioactive decay is already active in
  `QGSP_BIC_HP` — no prompt-decay override. *(Superseded — we built a minimal
  custom app, `stageA_transport/`, rather than stripping the example.)*

## Units convention

Positions **mm**, energies **keV** (beam energy **MeV**), times **ns**, dose
**Gy**. State units in HDF5 attributes.

## HDF5 schemas (the interface contracts)

### `prod_anh.h5` — Stage A output (one row per β⁺ emitter)
Dataset `/emitters` (columnar or compound):

| field        | type      | meaning                                   |
|--------------|-----------|-------------------------------------------|
| `isotope_id` | int8      | 0=¹⁵O 1=¹¹C 2=¹³N 3=¹⁰C 4=¹⁴O             |
| `prod_xyz`   | float32[3]| production point (mm) — the **truth** map |
| `anh_xyz`    | float32[3]| annihilation point (mm) — detector source |

> The two-point pair is consolidated into **one row** (cleaner in HDF5 than
> separate PROD/ANH records); positron range = `anh_xyz − prod_xyz`.

Root attributes (for normalization + reproducibility): `beam_energy_MeV`,
`phantom_material`, `n_protons_simulated`, `dose_Gy_total`, `Pj_per_proton`
(per isotope), `geant4_version`, `physics_list`, `random_seed`.

### `coincidences_<config>.h5` — Stage B output (one row per accepted coincidence)
Dataset `/coincidences`:

| field      | type       | meaning                          |
|------------|------------|----------------------------------|
| `hit1_xyz` | float32[3] | DOI-resolved 3-D position (mm)   |
| `hit2_xyz` | float32[3] | DOI-resolved 3-D position (mm)   |
| `e1`,`e2`  | float32    | deposited energies (keV)         |
| `t1`,`t2`  | float64    | timestamps (ns) — for TOF        |
| `truth`    | int8       | optional: 0 true,1 scatter,2 rand|

Root attributes: `detector_config`, geometry params, `energy_window_keV`,
`coinc_time_window_ns`, `Nj_budget` (per isotope), `realization_index`,
`random_seed`.

## Fixed parameters

**Beam (gun):** proton; baseline **100 MeV** (range ≈ 7.7 cm in water),
cross-check **150 MeV**; Gaussian pencil σ ≈ 3–5 mm; along cylinder axis;
fluence normalized to **1 Gy** at target (or fixed N_p with dose reported).

**Phantom (geometry):** **PMMA** C₅H₈O₂, ρ = 1.18 g/cm³; cylinder **Ø20 cm ×
20 cm**, axis along beam; present as passive attenuator in **every** Stage-B
geometry. Water (pure ¹⁵O source) is an optional cross-check.

**Isotopes / half-lives:** ¹⁵O 122 s, ¹¹C 1223 s, ¹³N 598 s, ¹⁰C 19.3 s,
¹⁴O 70.6 s (100 % β⁺ assumed).

**In-room acquisition:** delay **t_del ≳ 120 s** (MGH mean ≈ 150 s); window from
the measured-window equation (spec §3). At this delay the source is
**¹⁵O-dominated** → positron-range floor matters.

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
stageA_transport/    # Geant4 C++: protons → prod_anh.h5
handoff/             # Python: time-decay bookkeeping + N_j sampling
stageB_detector/     # Geant4 C++: annihilation events → coincidences_*.h5
reconstruction/      # deferred (Stage C)
common/              # shared HDF5 schema defs, units, isotope table
docs/                # copy of simulate_pt_pet.tex (the spec)
data/                # generated HDF5 (gitignored)
```

## Build / run

**Toolchain (this machine).** Geant4 **11.4.1** at `~/Software/geant4/install`,
built multithreaded (tasking) with **Qt visualization ON** (ToolsSG-Qt drivers:
`TSG_QT_ZB` software / `TSG_QT_GLES` GPU; no classic OpenGL). No ROOT. Activate
with `source $HOME/Software/geant4/install/bin/geant4.sh` (in `~/.zshrc`). HDF5
via Homebrew (`/opt/homebrew/opt/hdf5`); HighFive fetched by CMake when needed.

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
                              #   draws the PMMA phantom + shoots 20 protons
./proton_transport run.mac    # batch (18 threads); macros/run.mac smoke test
```

`build/compile_commands.json` (the `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON` flag)
feeds VS Code IntelliSense; `.vscode/settings.json` points the C++ extension at
it via CMake Tools.

### Handoff + analysis (Python)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r handoff/requirements.txt   # numpy scipy h5py
```

## First-session checklist

1. Read `docs/simulate_pt_pet.tex` end to end.
2. Stand up the `hadrontherapy`-based Stage-A app; confirm proton range and dose
   in the PMMA cylinder.
3. Add the `prod_anh.h5` writer (with the metadata attributes); rely on standard
   radioactive decay (no prompt-decay override); validate yields-per-proton
   against literature order-of-magnitude.
4. Implement the Python A→B handoff (Eqs. of spec §3) and verify the ¹⁵O/¹¹C mix
   vs. t_del.
5. Only then start Stage B (detector geometry + coincidence list).
