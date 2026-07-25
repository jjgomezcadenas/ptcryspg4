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

**Status: complete.** Results: (1) water slab, 2 mm, 100 MeV — O16 exposure
per proton 1.0010 of n×l (excess = secondary-proton path). (2) 4930 common
bins, maximum relative difference 5e-12 (CSV print precision); now the
tracked test `test/xsections/test_exposure_app.py`. (3) uniform_headep SOBP
vs `data/runs/uniform_headep_sobp_1e7`: target dose per proton agrees to
0.02%, core-dose ratio over the SOBP 0.9975 (min 0.987, max 1.005), distal
R80 difference 0.024 mm. (4) pilot
`data/xsections/exposure/uniform_headep_sobp_QGSP_BIC_HP_1e6` (regenerable,
gitignored): validator clean, 1000-replica fold in ~8 s; production-stage
R50 band ±0.16 mm (all_production), N13 alone ±0.71 mm.

**Implementation note — step subdivision.** Booking a whole Geant4 step at
its midpoint aliases the depth profile with the multi-mm step length near the
entrance (a ±2x comb, found against the stageA pencil reference). Every tally
therefore subdivides steps into segments shorter than half a bin in depth
(and in energy for the exposure), interpolating position and kinetic energy
along the path. This is the "tested by subdividing" rule of
`docs/xsections_plan.tex`, adopted permanently.

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

**Design decisions.** Two-pass sampling: the proton run samples production
points from its own random engine and writes `sampled_productions.csv`; a
separate emitter-transport mode (`--emitters-in`) creates each residual ion
at rest at its point and records the positron annihilation into
`emitters.csv`. The nominal curves reach the app as one flat file,
`data/xsections/fits/sampling_curves.csv` (tracked, written by
`export_sampling_curves.py`). Run-to-run reproducibility required
pre-generating the Geant4 per-event seed table
(`SetSeedOncePerCommunication(0)`): with the default per-communication
seeding, two identical runs did not reproduce. With the fix, sampler-on and
sampler-off runs give byte-identical depth-dose and native counters, and
exposure tables equal to summation round-off (~6e-12) — acceptance 4.
Acceptance 1: slab pulls -0.5σ / +0.2σ / +1.2σ on the three O16 channels.

**Status: complete.** Acceptance 2 at 1e7 protons: sampled/folded yields
0.9997 (O15), 1.0056 (C11), 1.0004 (N13); the folded R50 matches the
sampled R50 within its Poisson error for every isotope. Acceptance 3:
median positron ranges 1.920/0.764/1.073 mm vs the stageA reference
1.918/0.751/1.105 mm (O15/C11/N13). Emitter transport captured 99.87% of
sampled productions (the rest decay without a positron).

**Headline result** (run `uniform_headep_sobp_QGSP_BIC_HP_1e7`, reported in
`docs/sampling_xsections.tex`): Delta = R50(data) − R50(G4) per isotope:
O15 −1.48 ± 0.44 mm; C11 +5.36 ± 0.95 mm (the BIC near-threshold deficit
priced in millimetres); N13 +14.3 ± 4.1 mm (native yield low by 2.5x);
combined +0.36 ± 0.54 mm — the per-isotope shifts compensate at this field
and isotope mix. Yield ratios data/G4: 0.99 / 1.13 / 2.45 / 1.10.

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
