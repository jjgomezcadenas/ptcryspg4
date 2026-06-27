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
Realistic field through the heterogeneous head, then freeze the scenario.
- [ ] SOBP for the bone+brain path: `sobp.py` assumes one water-equivalent
      density, so re-tune (`mu`) or move to NLLS weights against the simulated
      heterogeneous depth-dose.
- [ ] Define the target at the skull base; dose normalization for the head.
- [ ] Verify plateau across the target through the head path.
- [ ] Freeze `head_mird_1e7` into `ptcrysp-scenarios` (build PDFs, snapshot).
- Acceptance: acceptable plateau; scenario frozen + pushed.

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

## Status

| phase | branch | status |
|---|---|---|
| 1 | `geometry-proof` | **done, merged to main** (builds, no overlaps; ~1.9k emitters, all in-head; skull O15/C11=1.37 vs brain 2.17) |
| 2 | `heterogeneous-medium` | done — pending review (phantom_regions.csv + per-region μ; bone added; explicit `geometry` label; **3 cases**: cylinder / uniform_head / mird_head, same head envelope; docs + data-driven plot) |
| 3 | `SOBP` | not started |

Update this table and tick the boxes as each phase lands.
