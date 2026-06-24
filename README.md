# ptcryspg4 — proton-therapy PET source generation

Produces the **positron-emitter source** for PET-based **in-room** proton-therapy
range verification. A Geant4 proton run sends protons through a phantom and records
where β⁺ emitters annihilate; a Python handoff turns the production into the number
of decays a scanner would measure. The output is **detector-independent** — a
separate downstream simulation reads it and models a PET detector.

## Documents

| File | Role |
|------|------|
| `latex/01_user_guide.tex` | **User guide** — physics and pipeline (read first) |
| `latex/04_source_reference.tex` | Consumer interface contract: the scenario data product |
| `CLAUDE.md` | Orientation + implementation decisions |
| `common/SCHEMA.md` | CSV file formats, isotope encoding, units |

The full doc set (numbered in reading order: `01_user_guide`, `02_beam_design`,
`03_decay_kinetics`, `04_source_reference`) lives in `latex/`.

## Pipeline

```
[A]  Geant4 transport    protons → phantom → β+ emitter → annihilation     RUNS ONCE
         ⟹  data/emitters.csv + data/run_meta.csv
[B0] Handoff (Python)    P_j → N_j(t_del)
         ⟹  data/sampling_budget_<scenario>.csv
       ─── frozen as a named scenario in ptcrysp-scenarios; read downstream ───
```

Stage A output is frozen as a named scenario in the `ptcrysp-scenarios` data repo;
a downstream PET detector simulation reads the scenario.

## Layout

```
stageA_transport/    Geant4 C++: proton transport, β+ decay capture, emitters.csv writer
field_design/        Python: SOBP beam design + depth-dose plots
decay_sampling/      Python: time-decay budget (budget.py) + realizations (budget_gen.py)
analysis_transport/  Python: validate Stage A output (dashboard + diagnostics)
tools/               snapshot_scenario.py: freeze a run into the scenarios repo
common/              schemas + isotope table (C++/Python mirrors)
latex/               LaTeX docs (01_user_guide … 04_source_reference) + figures + biblio
docs/                reference papers only (.pdf, .txt)
data/                Generated CSV files (gitignored)
```

Detector + reconstruction (Stages B, C) are a separate downstream repo.

## Requirements

- **Geant4 ≥ 11.4** built with multithreading
- **CMake ≥ 3.20**, Ninja (or make)
- **Python 3** with `numpy scipy pandas matplotlib` (system Python; see
  `analysis_transport/requirements.txt`). Output is CSV — no HDF5 dependency.

## Build / run

See **CLAUDE.md → Build / run** for the current commands (build `stageA_transport`,
run `proton_transport`, validate with `analysis_transport/validate_transport.py`).
