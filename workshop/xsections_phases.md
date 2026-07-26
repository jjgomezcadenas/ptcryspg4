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

## Fit revision — level anchoring (triggered by the N13 check)

**Trigger.** The planned check of the O-16(p,x)N-13 fit found the curve
1.85x above four mutually consistent plateau campaigns, insensitive to the
smoothing and spread parameters.

**Diagnosis.** In the campaign-offset model every campaign constrains the
curve *level* with equal weight regardless of precision; where one region's
campaigns are numerous and discrepant (the N13 resonances), they set the
level everywhere through the shared spline. Three of five channels were
affected (data/fit above 40 MeV: N13-from-O16 0.54, N13-from-N14 0.78,
C11-from-O16 1.32).

**Fixes (fit v4).**
1. Campaigns with a documented normalization uncertainty anchor the curve
   level at that value (covariance, CV selection and replica draws); values
   with citations in `config/xsection_fit.toml`: Akagi 2013 (EXFOR
   ERR-ANALYS, 4.1-4.9%) and Rodriguez-Gonzalez 2022/2023 (published
   budgets, 4.6-5.5%; paper in `papers/`).
2. The Akagi anchor is withdrawn in C11-from-O16 only: its level is
   contradicted there by the monitor-validated 2023 campaign and every
   independent campaign (its points stay, unanchored).
3. O-16(p,x)N-13 is fitted in two segments joined at 22 MeV (resonance /
   plateau), blended over 2 MeV, each with its own campaign-offset model.

**Result.** All five channels sit on their high-energy data
(data/fit 1.04 / 0.97 / 0.98 / 0.98 / 0.97). Bands tightened (production
median half-widths 2.6-4.8%). Method documented in `docs/xsection_fit.tex`.

**Updated headline numbers** (1e7 rerun, `docs/sampling_xsections.tex`):
Delta R50 O15 −1.78 ± 0.46, C11 +5.02 ± 0.91, N13 +28.4 ± 4.5 (its
measured production persists to the lowest energies, reaching deepest),
combined −0.03 ± 0.62; scenario-weighted +0.11 (fast), +0.19 (in-room),
+2.67 ± 0.58 (offline). Yield ratios data/G4: 0.97 / 1.32 / 1.49 / 1.10.
Replica band on the data-driven production edge: ±0.13 mm.

**Status: complete.**

## Phase 2 — energy-grid convergence and the frozen grid

**Goal.** Decide, by measurement, how coarse the exposure energy bins may be
before the folding results change, and freeze that grid for production.

**Method.** The 1e7 fine-grid (0.5 MeV) run is coarsened offline and
losslessly (exposures and moments re-summed); each candidate (1, 2, 5 MeV)
is folded with the nominal curve and all 1000 replicas and compared to the
fine grid replica by replica, so only the binning effect remains.

**Result.** 1.0 MeV passes every criterion for every profile with ~two
orders of magnitude of margin (largest paired R50 change 0.002 mm against
tolerances of 0.011–0.052 mm). 2.0 MeV fails on the N13 edge (0.134 mm vs
0.052 mm — the resonance structure), 5 MeV fails broadly. Frozen grid:
**1.0 MeV** (`frozen_energy_width_MeV` in the convergence config); decision
record promoted to `data/xsections/convergence/`. The production-stage
b_prod replica band at this field: ±0.13 mm (all_production and in-room),
per-isotope ±0.11 (C11), ±0.17 (O15), ±0.52 mm (N13).

**Status: complete.**

## Phase 3a — scenario assembly for the data-driven source

**Goal.** A sampling run becomes a complete, detector-consumable scenario
directory with no detector-side code changes.

**Done.** `make_scenario.py` writes `run_meta.csv` (stageA columns +
additive provenance: production model, sampling-curves and fit digests),
`phantom_regions.csv`, `isotopes.csv`, and drives `budget.py` for the
inroom/fast/offline budgets. The exposure app zeroes `dose_core` on air
bins (stageA convention) via geometry-aware on-axis containment.

**Acceptance.** `check_run.py`: all 13 checks pass on the 1e7 data-driven
run. PTCryspMC `load_scenario` reads the directory unchanged: pools
127007/76833/14237 (O15/C11/N13), in-room 1 Gy budgets 3.81e7 / 2.43e7 /
6.28e6 measured decays — the +31% C11 and +50% N13 against the native
scenario reproduce the comparison's yield ratios — and a 6.8e7-event
realization materializes.

**Status: complete.**

## Phase 3b — the source bank

**Goal.** One transported sample serving all 1000 replicas: candidates kept
with probability alpha x (n x l) x sigma_env, stored with q, unbiased after
division by q (docs/sampling_xsections.tex, Sec. source bank).

**Done.** Envelope column in sampling_curves.csv (nominal x replica cover
factor per channel); bank draw in the stepping action (same private engine;
transport untouched); source_bank.csv + bank metadata; validate_bank.py
(closure + ESS floors, bank_validation.json).

**Acceptance (1e7 field, alpha=2.48, 1,000,670 candidates).** Closure:
bank/fold nominal ratios 0.999-1.004 per channel, worst pull over all 1000
replicas x 5 channels = 1.5 sigma. ESS: totals equal entry counts to <1%
(weights near-uniform under the envelope); distal-window ESS 2,200-47,000
against the 500 floor. The distal window is anchored to the 99th-percentile
production depth, not the deepest stray entry.

**Status: complete.**

## Phase 3c — the 1e8 production run and the frozen data-driven scenario

**Done.** One 1e8-proton transport on the frozen 1.0 MeV grid with sampler
and bank (alpha 2.48): 2,194,871 sampled emitters, 10,001,867 bank
candidates; emitter transport, finalize, 1000-replica fold, comparison,
scenario assembly; check_run 13/13. Bank certified at scale (chunked
validator): closure ratios 0.997-1.002, worst pull 3.5 sigma over 5,000
checks (the expected extreme), distal-window ESS 22k-465k. Scenario frozen
as `ptcrysp-scenarios/scenarios/uniform_headep_sobp_1e8_dd` (26 files,
figures included; scenarios-repo commit d6f887f).

**Final comparison at 1e8** (in docs/sampling_xsections.tex): Delta R50 =
-1.72 +- 0.17 (O15), +4.58 +- 0.38 (C11), +23.35 +- 1.45 (N13),
+0.05 +- 0.18 (combined); scenario-weighted -0.19 +- 0.16 (fast),
+0.33 +- 0.19 (in-room), +2.66 +- 0.22 (offline, 12 sigma). Fold closure
0.999-1.001. u_xs bands unchanged: +-0.13 mm production, +-0.12 mm in-room.

**Status: complete.**

## Scenario modernization — CBS acquisition protocols

**Trigger.** The fast/inroom/offline handoff scenarios (20-30 min windows)
predate the CRYSP study; cbs.tex defines short modulated scans by start
delay and duration (reference: 300 s scan starting 120 s after beam-off).

**Done.** config/handoff_scenarios.toml replaced by five named protocols
(t_irr 60 s): d120s300 (reference), d120s120, d180s300, d180s120,
d300s300. Every scenario reference repointed (folding, comparison,
budgets, native routes, plots, tests, analytic reference); the comparison
table carries a generated timing legend so labels cannot drift from the
numbers. The 1e8 products were refolded and recompared; budgets re-frozen
into uniform_headep_sobp_1e8_dd (scenarios-repo commit 4323c71).

**Physics.** Under CBS weighting the combined displacement runs from
-1.08 +- 0.16 mm (d120s120) through -0.80 +- 0.17 (reference) to
+0.32 +- 0.19 mm (d300s300) - crossing zero inside the studied range;
short windows weight O15, later starts restore C11. Scenario-weighted
yield ratios 1.02-1.11. Reference-scan u_xs band: +-0.14 mm. The note's
Results and Conclusions carry the new numbers, with the
displacement-not-error phrasing made explicit.

**Status: complete.**

## Later phases (scoped when reached)
- Phase 3d: detector-study reruns on the data-driven scenario and the
  native-vs-data sigma_R cross-check (PTCryspMC, unweighted pipeline).
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
