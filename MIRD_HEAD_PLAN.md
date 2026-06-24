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
| Head (outer) | elliptical tube ∪ ellipsoid cap | 7 × 10 × 7.75 (tube) | `G4_TISSUE_SOFT_ICRP` |
| Skull | ellipsoid(6.8×9.8×8.3) − ellipsoid(6×9×6.5) | ~0.8 cm cranium shell | `G4_BONE_CORTICAL_ICRP` |
| Brain | ellipsoid | 6 × 9 × 6.5 | `G4_BRAIN_ICRP` |

Full head ≈ 14 (L-R) × 20 (A-P) × 16 (S-I) cm. Lateral path ≈ 14 cm.

## Phased plan

### Phase 1 — branch `geometry-proof`
Geometry only; prove it works and is more realistic; **no schema changes**.
- [ ] Replicate the 3-region MIRD head (head/skull/brain) in
      `DetectorConstruction`, re-centered at origin, oriented for a lateral beam.
- [ ] NIST materials (brain/skull/scalp as above); confirm the bone choice.
- [ ] Make geometry selectable: `/stageA/phantom/geometry cylinder|mird_head`
      (cylinder default). New messenger command.
- [ ] Single pencil beam, lateral entry.
- [ ] Run; verify activity sits in brain + skull, the bone changes the local mix,
      yields are sane, no overlaps/crashes.
- Files: `DetectorConstruction.{cc,hh}`, `DetectorMessenger.{cc,hh}`,
  `StageAConfig.hh` (head dims/material constants), `PrimaryGeneratorAction`
  (entry point for the head extent).
- Keep `run_meta.csv` columns unchanged (label the geometry; full per-region
  medium is Phase 2). Do **not** touch `phantom_material.py` / `SCHEMA.md`.
- Acceptance: builds, runs, `emitters.csv` shows production inside the head;
  TrackingAction (Z,A capture) is geometry-independent so should just work.

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
| 1 | `geometry-proof` | not started |
| 2 | `heterogeneous-medium` | not started |
| 3 | `SOBP` | not started |

Update this table and tick the boxes as each phase lands.
