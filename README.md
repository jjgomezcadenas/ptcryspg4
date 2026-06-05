# ptcryspg4 — PET detector comparison for proton-therapy range verification

Simulation chain that compares candidate PET detectors (the CRYSP family vs.
conventional LYSO/BGO) for **in-room** proton-therapy range verification. A Geant4
proton run (this repo) produces the positron-emitter source; an analytic detector
Monte Carlo in Julia (`PTCryspMC.jl`, a separate repo) turns it into a coincidence
list per detector. The figure of merit is σ(range) — the precision of the
recovered distal fall-off — at photon-starved in-room statistics.

## Documents

| File | Role |
|------|------|
| `docs/simulate_pt_pet.tex` | **Spec** — physics and pipeline |
| `CLAUDE.md` | Orientation + implementation decisions |
| `common/SCHEMA.md` | CSV file formats, isotope encoding, units |

## Pipeline

```
ptcryspg4 (this repo) — RUNS ONCE
  [A]  Geant4 transport    protons → phantom → β+ emitter → annihilation
           ⟹  data/emitters.csv + data/run_meta.csv
  [B0] Handoff (Python)    P_j → N_j(t_del)
           ⟹  data/sampling_budget_<scenario>.csv
PTCryspMC.jl (separate repo) — RUNS PER DETECTOR
  [B]  Analytic detector   events → detector response → coincidences
           ⟹  coincidences_<config>.csv
  [C]  Reconstruction      coincidence list → image → σ(range)     DEFERRED
```

Stage A output is frozen as a named scenario in the `ptcrysp-scenarios` data repo;
`PTCryspMC.jl` reads from there.

## Layout

```
stageA_transport/    Geant4 C++: proton transport, β+ decay capture, emitters.csv writer
field_design/        Python: SOBP beam design + depth-dose plots
decay_sampling/      Python: time-decay budget (budget.py) + realizations (budget_gen.py)
analysis_transport/  Python: validate Stage A output (dashboard + diagnostics)
tools/               snapshot_scenario.py: freeze a run into the scenarios repo
common/              schemas + isotope table (C++/Python mirrors)
docs/                Spec + reference material
data/                Generated CSV files (gitignored)
```

Detector + reconstruction (Stages B, C) are in the separate `PTCryspMC.jl` repo.

## Requirements

- **Geant4 ≥ 11.4** built with multithreading
- **CMake ≥ 3.20**, Ninja (or make)
- **Python 3** with `numpy scipy pandas matplotlib` (system Python; see
  `analysis_transport/requirements.txt`). Output is CSV — no HDF5 dependency.

## Build / run

See **CLAUDE.md → Build / run** for the current commands (build `stageA_transport`,
run `proton_transport`, validate with `analysis_transport/validate_transport.py`).
