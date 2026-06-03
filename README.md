# ptcryspg4 — PET detector comparison for proton-therapy range verification

Geant4 simulation chain that compares candidate PET detectors (the CRYSP family
vs. conventional LYSO/BGO) for **in-room** proton-therapy range verification.
The deliverable is a coincidence list per detector; the figure of merit is
σ(range) — the precision of the recovered distal fall-off — at photon-starved
in-room statistics.

## Documents

| File | Role |
|------|------|
| `docs/simulate_pt_pet.tex` | **Authoritative spec** — physics and pipeline |
| `CLAUDE.md` | Orientation + implementation decisions |
| `common/SCHEMA.md` | Interface contracts (CSV file formats, isotope encoding, units) |

## Pipeline

```
[A]  Geant4 transport    protons → PMMA → β+ emitter → annihilation     RUNS ONCE
         ⟹  data/emitters.csv + data/run_meta.csv
[B0] Handoff (Python)    P_j → N_j(t_del); sample annihilation events
         ⟹  data/sampled_events_*.csv
[B]  Geant4 detector     events → detector response → coincidences      RUNS PER DETECTOR
         ⟹  data/coincidences_<config>.csv
[C]  Reconstruction      coincidence list → image → σ(range)            DEFERRED
```

## Layout

```
stageA_transport/    Geant4 C++: proton transport, β+ decay capture, emitters.csv writer
decay_sampling/      Python: time-decay bookkeeping (spec §3) + N_j sampling (Stage B0)
stageB_detector/     Geant4 C++: detector response + coincidence sorting
reconstruction/      Stage C (deferred)
analysis_transport/  Python: validate Stage A output (dashboard + diagnostics)
common/              Interface contracts: schemas + isotope table (C++/Python mirrors)
docs/                Spec + reference material
data/                Generated CSV files (gitignored)
```

## Requirements

- **Geant4 ≥ 11.4** built with multithreading
- **CMake ≥ 3.20**, Ninja (or make)
- **Python 3** with `numpy scipy pandas matplotlib` (system Python; see
  `analysis_transport/requirements.txt`). Output is CSV — no HDF5 dependency.

## Build / run

See **CLAUDE.md → Build / run** for the current commands (build `stageA_transport`,
run `proton_transport`, validate with `analysis_transport/validate_transport.py`).
