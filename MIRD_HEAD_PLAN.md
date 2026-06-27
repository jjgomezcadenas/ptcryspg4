# Plan — MIRD-head phantom variant

A more realistic, **heterogeneous** head phantom for the Parodi skull-base case,
as an alternative to the current homogeneous `G4_BRAIN_ICRP` cylinder. Built from
the stylized MIRD head shipped with Geant4 (`advanced/human_phantom`). The
cylinder stays the **default**; the MIRD head is a selectable **variant**.

This is a living planning doc: work proceeds in **three branches, in sequence**,
each reviewed before the next. Do not start a phase without an explicit go-ahead.

## Why

Parodi 2008 did **not** use a cylinder: the β⁺ maps come from CT-based MC of a
*real patient* head (heterogeneous — brain, skull bone, fat, cavities). The
homogeneous cylinder captures bulk yield but not the skull (bone changes proton
range, the O/C mix, and 511 keV attenuation) or fat-driven spatial structure.
The MIRD head is the middle ground: a real anatomical head (brain in a skull),
analytic, no patient CT needed, ships with Geant4.

## Fixed decisions (agreed)

- **Beam orientation:** lateral — the head is oriented so the beam (`+z` in our
  frame) enters the side and crosses skull → brain → skull (~14 cm), matching
  Parodi's lateral/oblique portals.
- **Anatomy:** head only — scalp/face (soft tissue) + skull (bone) + brain.
  **No spine** (that would be the paraspinal variant).
- **Beam type:** single **pencil first** (Phase 1 shape/feasibility), **SOBP
  later** (Phase 3).
- **Materials:** substitute NIST names (not the example's custom materials), so
  the medium stays consistent with `common/phantom_material.py`:
  - brain → `G4_BRAIN_ICRP` (matches the cylinder)
  - skull → `G4_BONE_CORTICAL_ICRP` (bone choice to confirm at Phase 1)
  - scalp/face → `G4_TISSUE_SOFT_ICRP`
- **Default unchanged:** cylinder remains the default geometry.

## MIRD head geometry (from `human_phantom` source, to replicate ~60 lines)

Three nested regions; we replicate the solids rather than import the example's
class framework. Frame in the example: `x` = left-right, `y` = anterior-posterior,
`z` = superior-inferior. We re-center the head at the origin (drop the example's
+77.75 cm body offset) and orient it for a lateral beam.

| region | solid | semi-axes (cm) | our material |
|---|---|---|---|
| Head (outer) | single scalp ellipsoid (skull-outer + ~0.4 cm) | 7.2 × 10.2 × 8.7 | `G4_TISSUE_SOFT_ICRP` |
| Skull | ellipsoid(6.8×9.8×8.3) − ellipsoid(6×9×6.5) | ~0.8 cm cranium shell | `G4_BONE_CORTICAL_ICRP` |
| Brain | ellipsoid (+1 cm crown offset) | 6 × 9 × 6.5 | `G4_BRAIN_ICRP` |

The MIRD **face/neck (elliptical tube) is dropped**: the lateral cranial field
never crosses it (≈0 emitters there), and it only served to enclose the lower
skull — so the outer head is a single soft-tissue scalp ellipsoid enclosing the
whole cranium. Skull and brain stay the verbatim MIRD ellipsoids. Concentric
scalp ⊃ skull ⊃ brain, the beam crossing scalp → skull → brain (~14 cm).

## Phased plan

### Phase 1 — branch `geometry-proof`
Geometry only; prove it works and is more realistic; **no schema changes**.

Confirmed defaults: bone = `G4_BONE_CORTICAL_ICRP`; edep scorer attached to the
**head** LV (real skull-base target + dose normalization deferred to Phase 3);
**faithfully replicate** the example's relative solid placements (re-centered +
rotated), not idealized concentric ellipsoids.

- [ ] Replicate the 3-region MIRD head in `DetectorConstruction`:
      Head `EllipticalTube(7,10,7.75)∪Ellipsoid(7,10,8.5)` = `G4_TISSUE_SOFT_ICRP`;
      Skull `Ellipsoid(6.8,9.8,8.3)−Ellipsoid(6,9,6.5)` = `G4_BONE_CORTICAL_ICRP`;
      Brain `Ellipsoid(6,9,6.5)` = `G4_BRAIN_ICRP`. Re-center at origin (drop the
      +77.75 cm body offset); keep relative skull/brain placements; rely on G4
      overlap checking.
- [ ] Lateral orientation: rotate the head so its L-R axis (±7 cm) aligns with
      the beam `+z` (`rotateY(90°)`), so the pencil crosses skull→brain→skull.
- [ ] Geometry selectable: `/stageA/phantom/geometry cylinder|mird_head`
      (PreInit, cylinder default). `Construct()` branches into
      `BuildCylinder()` / `BuildMirdHead()`.
- [ ] Single pencil, lateral entry. Generalize the detector's beam-axis
      half-extent (cylinder→`fHalfZ`; head→outer L-R semi-axis) so the gun starts
      1 mm upstream for either geometry.
- [ ] Scorer on the head LV (`fPhantomLV` → head) so the existing edep scorer +
      `PhantomMass()` work unchanged.
- [ ] `emitters.csv` is the deliverable; `TrackingAction` (Z,A capture) is
      geometry-independent so it just works.
- [ ] **Detailed plots** — `analysis_transport/plot_mird_head.py` rendering:
      (1) the phantom (3 nested regions in two orthogonal cross-sections),
      (2) the beam (lateral pencil entry + path), and
      (3) the emitter paths (production + annihilation points and the prod→anh
      segments from `emitters.csv`, coloured by isotope / region).
      Mirrors the `StageAConfig` ellipsoid constants for now (Phase 2 writes them
      into `run_meta.csv`). This is the Phase-1 verification artifact.
- Files: `StageAConfig.hh` (head semi-axes + 3 materials + default geometry),
  `DetectorConstruction.{cc,hh}`, `DetectorMessenger.{cc,hh}`,
  `PrimaryGeneratorAction.cc`, `macros/mird_head.mac`,
  `analysis_transport/plot_mird_head.py`.
- Keep `run_meta.csv` columns unchanged (set `phantom_material="MIRD_head"`,
  provisional bounding dims; full per-region medium is Phase 2). Do **not** touch
  `phantom_material.py` / `SCHEMA.md`.
- Acceptance: builds; runs `mird_head.mac` (~10⁵ protons) with no overlap
  warnings; the plots show activity tracking the proton path through the head,
  with a visible skull-vs-brain isotope-mix difference.

### Phase 2 — branch `heterogeneous-medium`
Make the heterogeneous medium first-class (the deferred "multi-material /
attenuation-map" work, in its tractable analytic form).
- [ ] `run_meta.csv`: record per-region material + geometry (not one scalar).
- [ ] `common/phantom_material.py`: add bone; emit **per-region** composition + μ.
- [ ] `tools/snapshot_scenario.py`: per-region `phantom_material_*` files.
- [ ] `common/SCHEMA.md` + `latex/04_source_reference.tex`: document the
      multi-region medium; lift the "homogeneous-only" caveat for this variant.
- [ ] `analysis_transport/plot_geometry.py`, `validate_transport.py`: render /
      handle the head geometry (they currently assume a cylinder).
- Acceptance: a snapshot carries per-region medium; docs consistent.

### Phase 3 — branch `SOBP`

Drive a realistic spread-out field through the heterogeneous head, normalize to
1 Gy in a brain target, and freeze the scenario — using the methods real proton
planning uses rather than ad-hoc knobs: **water-equivalent path length (WEPL)**
for the bone heterogeneity, a **central-axis dose** profile for the shape, and
**R80 + plateau uniformity** for acceptance. Pencil *and* SOBP both selectable.

> **Pencil vs SOBP is already free.** `BeamConfig::SobpEnabled()` drives the run
> tag, so a head pencil run and a head SOBP run auto-separate into
> `data/runs/mird_head_pencil_*` and `…_sobp_*` with no clobber. Phase 3 only
> adds the SOBP macro(s); the existing pencil macros stay as the shape test.

#### Step 0 — fix the head depth reference (correctness; must land first)
- **What's wrong.** For the head, `fHalfZ` — the z-origin that the target box and
  the *reported* depths are measured from — is set only in `BuildCylinder` (to
  80 mm) and left there for the head, which actually spans ±72 mm along the beam
  (the scalp semi-axis) with the beam entering the scalp at z=−72. So
  `TargetProxZ = −80 + 55 = −25` places the target box **8 mm too shallow**, and
  `run_meta` *reports* 55–105 mm while the box physically sits at 47–97 mm from
  the scalp.
- **Why first.** Every "is the plateau flat over the target" judgement compares
  the realized dose to where the target box actually is. An 8 mm registration
  error silently invalidates the acceptance metric, so nothing downstream is
  trustworthy until this is fixed.
- **How.** In `BuildMirdHead`/`BuildUniformHead` set `fHalfZ = fBeamHalfExtent`
  (`= kScalpAxMM`). Then the entrance face (−`fHalfZ`), `TargetProxZ/DistZ`, the
  `depth_dose` registration, and the `run_meta` depths are all mutually
  consistent. Also fix `TargetMass()`: it uses "the phantom material", which for
  the head is the mother LV (soft-tissue **scalp**), but the box is in **brain** —
  look up the material at the box centre via the region list and use its density,
  so `dose = edep/mass` is dose-to-brain.
- **Files:** `DetectorConstruction.{cc,hh}`.
- **Done when:** a head run reports target depths that match the physical box, and
  target mass uses brain density.
- **Done (this branch).** Plus a regression net so this class of bug can't recur
  silently: `analysis_transport/check_run.py` (independent Python invariant
  checker — box-in-medium, target-mass density, dose ≥ whole-phantom, Np
  consistency, depths/regions in bounds; would have failed on both halves of the
  bug), `analysis_transport/plot_phantom.py` (regions + target box + beam in two
  cross-sections, wired into `make_figures`), and a C++ run-time guard warning
  when the box centre is in air.

#### Step 1 — central-axis dose profile (the curve we actually grade)
- **What's wrong.** `depth_dose.csv` tallies edep per z-bin integrated over the
  *whole* transverse plane. Through an ellipsoid that conflates the SOBP shape
  with (a) the cross-section shrinking toward the poles and (b) dose spikes where
  the axis crosses the dense bone shells — so it is **not** dose(z).
- **Why.** The standard central-axis depth dose is **dose-to-medium along a thin
  column on the beam axis**. Over the brain (≈constant density) a thin
  fixed-radius core's edep(z) is already ∝ dose, removing the cross-section
  contamination; this is the curve R80 and uniformity are read from.
- **How.** Add a parallel tally in `StageARun` for steps within a small radius
  `r_core` (≈5 mm) of the axis, binned in z over the head extent; write an
  `edep_core_MeV` column (and optionally `dose_core_Gy = edep / (ρ_bin·V_bin)`,
  `ρ_bin` from the region at the bin centre, so the bone bins are comparable too).
  `r_core` a `StageAConfig` constant.
- **Files:** `StageARun.{cc,hh}`, `RunAction.cc` (new column), `StageAConfig.hh`,
  `common/SCHEMA.md`.
- **Done when:** the core profile over the brain is a clean plateau/peak without
  the ellipsoid taper.

#### Step 2 — design the field in WEPL (the principled bone correction)
- **What's wrong.** `sobp.py` designs a flat SOBP in a *single water-equivalent*
  medium. The head path is soft tissue → ~0.8 cm bone → brain → bone → soft
  tissue; cortical bone's relative stopping power (RSP ≈ 1.6) makes 0.8 cm of
  skull ≈ 1.3 cm water-equivalent, so a water design lands the whole stack
  proximal.
- **Why.** Real planning never designs in geometric depth through bone — it
  converts to **water-equivalent path length** (`WEPL = ∫ RSP · dx`) and designs
  there, where the Bragg curve is universal. This absorbs the bone offset *by
  construction*; no `exp(µR)` fudge factor.
- **How.**
  1. **RSP per material** — relative stopping power vs water from each material's
     density + mean excitation energy (already in `phantom_material.py`) via the
     Bethe stopping-power ratio at a representative energy (~150 MeV; RSP is only
     weakly energy-dependent). Brain ≈ 1.03, cortical bone ≈ 1.6, soft tissue ≈ 1.0.
  2. **Ray-trace the central axis** (+z) through the priority-ordered ellipsoids
     in `phantom_regions.csv` to get each material's geometric thickness vs depth,
     and integrate RSP → `WEPL(geometric depth)`.
  3. **Map the geometric target window** (55–105 mm) to its WEPL window and hand
     that *radiological* window to `sobp.py`'s existing Bortfeld/Abel design. The
     layer energies then place Bragg peaks at the right WEPL — i.e. the right
     geometric depth in the head.
- **Files:** `field_design/sobp.py` (a heterogeneous/WEPL mode taking
  `phantom_regions.csv` + RSP; the water/cylinder path stays the default),
  `common/phantom_material.py` (expose RSP), `latex/02_beam_design.tex` (document
  the WEPL design).
- **Done when:** layers generated for the head's WEPL window; the design's nominal
  peak placement matches the target geometrically.

#### Step 3 — run + verify with R80 and plateau uniformity (the standard metrics)
- **Why these.** Range is specified as **R80** — the distal depth at 80 % of the
  plateau dose — because it is nearly independent of energy spread, hence robust
  and reproducible (not "where the peak is"). SOBP quality is **uniformity** over
  the modulation width. These are the field's acceptance.
- **How.** Run `mird_head_sobp` at 1e6 (design iteration) then 1e7 (freeze); same
  for `uniform_head_sobp` — the **no-bone control**: same field, isolates the
  skull's effect on the plateau and on R80. Extend `plot_sobp` to read the core
  profile, mark R80, shade the target, and print uniformity = spread/mean over the
  target window.
- **Acceptance.** R80 within a few mm of the target distal edge; uniformity ≲ ±5 %
  over the target (clinical is ±2–3 %; relaxed here because this is a
  detector-independent *source* and the activity map — the real deliverable —
  falls out of the MC regardless of plateau cosmetics).
- **Files:** `field_design/plot_sobp.py` (R80 + uniformity on the core profile),
  `macros/mird_head_sobp.mac`, `macros/uniform_head_sobp.mac`.
- **Done when:** both heads run, R80 + uniformity reported, plateau acceptable.

#### Step 4 — NNLS weight optimization (fallback; only if Step 3 isn't flat enough)
- **What.** If the analytic WEPL design still ripples (bone is not a perfectly
  constant offset; the distal shell, range straggling, and multiple scattering all
  smear it), fit the layer weights to the *simulated* response instead of the
  analytic one — a 1-D version of the inverse planning a real TPS does.
- **How.** Run each candidate energy as a single-energy disk through the head,
  tally its central-axis core depth-dose `D_i(z)` → response matrix; solve
  `min ‖Σ w_i D_i(z) − target‖²` s.t. `w ≥ 0` (`scipy.optimize.nnls`) over the
  target window; re-run with the optimized weights. Bounded: one response matrix +
  a solve.
- **Files:** `field_design/sobp.py` (or a `sobp_opt.py`) + a driver to collect the
  per-layer responses.
- **Done when:** optimized weights give acceptable uniformity.

#### Step 5 — normalize, freeze, document
- **Normalization.** Scale to **1 Gy in the (now correctly placed) brain target
  box**; the existing `P_j(D) = count_j · D / target_dose` machinery, target mass =
  brain.
- **Freeze.** `mird_head_sobp_1e7` (and optionally `uniform_head_sobp_1e7`) into
  `ptcrysp-scenarios`: `make_figures` (head plots + core SOBP plateau), build PDFs,
  `snapshot_scenario.py <run_tag>`.
- **Docs.** `latex/02_beam_design.tex` (WEPL design + R80), `04_source_reference`
  / `SCHEMA.md` (new `edep_core`/`dose_core` columns; the head field is
  WEPL-designed), this plan's status.
- **Done when:** scenario frozen + pushed; pencil and SOBP both selectable via
  macros; docs consistent.

**Acceptance (phase):** correctly registered brain target; clean central-axis
dose profile; R80 within a few mm of the target distal edge and uniformity ≲ ±5 %
through the head path; `mird_head_sobp_1e7` frozen + pushed; pencil and SOBP both
runnable.

> **Why this replaces the old A/B/C-with-µ framing.** Dropped the `exp(µR)`
> re-tune (no clinical analog). The bone correction is now **WEPL** (what planners
> actually do), the profile is **central-axis dose-to-medium** (the standard
> scorer, not transverse-integrated edep), acceptance is **R80 + uniformity** (the
> standard metrics), and **NNLS weight optimization** is the rigorous fallback
> (1-D inverse planning) — not an exotic escalation.

## Cross-cutting caveats

- **Heterogeneous breaks the homogeneous invariant** baked into the pipeline
  (SCHEMA File 3, `phantom_material.py`, `04_source_reference`). Phase 2 is where
  that's resolved — for the analytic 3-region head, not full voxel CT.
- **Lose the clean Parodi Table-2 integral cross-check** (different geometry);
  the bone-driven O15/C11 change becomes the interesting output instead.
- **Biological washout** is still not modeled (out of scope unless added later).
- **MIRD head is stylized**, not a patient CT. A voxelized-CT head
  (`extended/medical/DICOM`) is the faithful-but-much-larger future option, not
  part of this plan.
- **WEPL/RSP is analytic, not CT-calibrated.** Phase 3 computes relative stopping
  powers from the analytic materials (Bethe ratio), and ray-traces WEPL through
  the analytic ellipsoids — the stylized-phantom analog of a TPS's HU→RSP
  calibration + CT ray-cast, good enough to place the field, not a clinical dose
  calc.

## Status

| phase | branch | status |
|---|---|---|
| 1 | `geometry-proof` | **done, merged to main** (builds, no overlaps; ~1.9k emitters, all in-head; skull O15/C11=1.37 vs brain 2.17) |
| 2 | `heterogeneous-medium` | **done, merged to main** (phantom_regions.csv + per-region μ; bone added; explicit `geometry` label; **3 cases**: cylinder / uniform_head / mird_head, same head envelope; docs + data-driven plot) |
| — | `data-reorg` | **done, merged to main** (per-run `data/runs/<tag>/` dirs; geometry-aware `make_figures.py`; snapshot tied to a run tag — the infra Phase 3 freezes through) |
| 3 | `SOBP` | **scoped** (WEPL design + central-axis dose + R80/uniformity; pencil & SOBP selectable). Awaiting OK to implement. |

Update this table and tick the boxes as each phase lands.
