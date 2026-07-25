# Cross-section work — phase log

Running record of the implementation phases on the `xsection-source-bank`
branch (this repo) and the `xsection-weighted-lors` branch (PTCryspMC.jl).
One section per phase: goal, deliverables, acceptance, status, recovery notes.
Update the status line of a phase in the same commit that completes it.

Governing documents: `docs/shared_plan.tex` (cross-repository contract),
`docs/xsections_plan.tex` (Geant4 + folding implementation),
`docs/sampling_xsections.tex` (the physics of the sampling method and the
R50 comparison, kept pedagogical and bank-free).

## Phase 0 — branches and documentation

**Goal.** Both repos carry their feature branch and document it.

**Done.** PTCryspMC.jl: `position-resolution` merged into `main` (pushed),
`xsection-weighted-lors` branched from it; AGENTS.md created; CLAUDE.md points
to it (commit `f3d3892`). ptcryspg4: AGENTS.md records both branches;
CLAUDE.md points to it; weighted-route and randoms-approximation sentences
added to the two plan documents (commit `72864dc`).

**Status: complete.**

## Phase 1a — new Geant4 application, exposure mode

**Goal.** A standalone app `stageA_xsection_ensemble/` that transports the
standard SOBP field and writes the target-resolved proton exposure table:
for every (target isotope, proton-energy bin, depth bin), the sum of
weight × nuclei/cm³ × step length over all proton steps, in 1/cm².
Multiplying a bin by a cross section in cm² gives expected production there.

**Deliverables.**
- `stageA_xsection_ensemble/` with its own CMake target; geometry and SOBP
  inputs read from the same versioned files as `stageA_transport`; physics
  list selectable at run time (`QGSP_BIC_HP` default).
- Exposure scorer on the fine energy grid of
  `config/xsection_exposure_convergence.toml`.
- `exposure_meta.json` (proton count, target dose, Np_per_Gy, physics list,
  seed, bin edges — the contract `exposure_folding.py` validates).
- `depth_dose.csv`; native-route counters (all β+ residuals produced by
  Geant4's own model, by projectile / target / residual / depth — diagnostic
  only).

**Acceptance.**
1. Thin uniform slab: exposure per proton equals nuclei/cm³ × thickness,
   computable on paper.
2. Same synthetic step list into the C++ scorer and the Python reference
   accumulator (`accumulate_step_exposure`): identical tables to
   floating-point precision.
3. Matched-settings run vs `stageA_transport`: target dose and depth-dose
   shape agree statistically.
4. ~1e6-proton pilot: exposure table and metadata pass the Python validator
   and fold cleanly with the existing layer.

**Status: pending.**

## Phase 1b — in-flight nominal sampler and the R50 comparison

**Goal.** In the same app, sample emitter production directly from the fitted
nominal curves during transport: per proton step and open channel, an emitter
is created with probability weight × n × ℓ × σ(E); the isotope decays
naturally, the positron annihilates, and an ordinary `emitters.csv` is
written. The native tally counts what Geant4's model would have produced on
the identical protons. Headline result: the paired distal-edge comparison
Δ = R50(data-driven) − R50(G4), per isotope and combined, reported in
`docs/sampling_xsections.tex`.

**Deliverables.**
- Sampling action with its own random-number stream (transport histories
  stay bit-identical with the sampler on or off).
- Data-driven `emitters.csv` + provenance (fit digest, curve version);
  native-route counters from the same run.
- Tracked analysis producing both production profiles, per-Gy yields, and
  the Δ table + figure into `docs/sampling_xsections.tex`.

**Acceptance.**
1. Mono-energetic slab: sampled yield reproduces n × ℓ × σ(E) analytically.
2. Sampled production profile agrees with folding the same run's exposure
   table by the nominal curves (same physics, two routes).
3. Positron-range distributions match the per-isotope references.
4. Exposure-only and sampling runs with equal seeds produce identical
   depth-dose (transport untouched by sampling).

**Status: pending.**

## Later phases (scoped when reached)

- Phase 2: fine-grid pilot, energy-bin convergence study, frozen grid,
  production folding of all 1000 replicas (b_prod distribution).
- Phase 3: source bank writer + variant-source mode
  (`docs/xsections_plan.tex`).
- Phase 4: PTCryspMC.jl — second parent id on randoms, bank source mode,
  weights evaluator, resampler (`docs/shared_plan.tex`;
  `PTCryspMC.jl/dev/xsection_weighted_lors_plan.md`).
- Phase 5: end-to-end replica propagation and the hadronic-transport
  envelope.

## Recovery notes

- Fit products: regenerate figures/tables with
  `python3 -m analysis_transport.xsections.make_fit_plots`, verify with
  `python3 -m analysis_transport.xsections.validate_fit`.
- Folding layer: tests under `test/xsections/`
  (`python3 -m unittest discover -s test/xsections`).
- Documents: `pdflatex` from `docs/` (PDFs are build products, gitignored).
- Each phase's completing commit is recorded in its status line; `git log`
  on the branch is the authoritative sequence.
