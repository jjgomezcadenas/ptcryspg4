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
| `common/SCHEMA.md` | HDF5 interface contracts (file formats, isotope encoding, units) |

## Pipeline

```
[A]  Geant4 transport    protons → PMMA → β+ emitter → annihilation     RUNS ONCE
         ⟹  data/prod_anh.h5
[B0] Handoff (Python)    P_j → N_j(t_del); sample annihilation events
         ⟹  data/sampled_events_*.h5
[B]  Geant4 detector     events → detector response → coincidences      RUNS PER DETECTOR
         ⟹  data/coincidences_<config>.h5
[C]  Reconstruction      coincidence list → image → σ(range)            DEFERRED
```

## Layout

```
stageA_transport/    Geant4 C++: proton transport, β+ decay capture, prod_anh.h5 writer
handoff/             Python: time-decay bookkeeping (spec §3) + N_j sampling
stageB_detector/     Geant4 C++: detector response + coincidence sorting
reconstruction/      Stage C (deferred)
common/              Interface contracts: HDF5 schemas + isotope table (C++/Python mirrors)
docs/                Spec + reference material
data/                Generated HDF5 files (gitignored)
```

## Requirements

- **Geant4 ≥ 11.4** built with multithreading (no Qt/OpenGL/GDML needed)
- **CMake ≥ 3.20**, Ninja (or make)
- **HDF5** C library (`brew install hdf5`); HighFive is fetched automatically by CMake
- **Python ≥ 3.10** with `numpy scipy h5py` (see `handoff/requirements.txt`)

## Build / run

*(TBD — filled in as code lands.)*
