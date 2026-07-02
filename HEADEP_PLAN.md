# HEADEP_PLAN.md — a MIRD head with an ependymoma, treated by a posterior beam

Plan for a new phantom, `headep`: the MIRD head plus a posterior-fossa tumour
(ependymoma, from `~/Papers/CryspBrain/figs/ependymoma_mri.png`), irradiated by a
straight posterior beam aimed at the tumour. The existing `mird_head` (skull +
lateral beam, brain-centred target) stays unchanged.

Design choices, already fixed with the user:
- **(A) baked**: `headep` *is* "tumour + posterior beam aimed at it". Beam
  direction is not a separate knob; it is intrinsic to the geometry. The run tag
  `headep_sobp_<N>` is therefore unique.
- **tumour-as-target**: the dose target (and the SOBP the field is designed for)
  is the tumour, not the generic 5 cm box.
- **straight posterior**: an axis-aligned occiput→anterior field (no obliquity).

## The mechanism this rests on

The gun always fires along **world +z**, centred on the axis (disk/pencil around
x=y=0), starting just upstream of the head. The target box, the WEPL ray-trace,
and the dose scoring are all defined along +z. So the beam direction lives entirely
in **how the head is rotated and translated** relative to that fixed +z axis:

- **orientation** — which head axis lies along +z: `lateral` = L-R (current,
  `RotateY 90°`); `posterior` = A-P (`RotateX 90°`, beam enters the occiput
  travelling anterior).
- **on-axis translation** — what sits on the beam axis at the target depth: the
  base case centres the brain; `headep` centres the tumour.

Consequence: to irradiate an off-axis tumour, the gun, the target box, and
`sobp.py` are not touched — the head is placed so the tumour sits on the beam axis,
and the whole +z pipeline hits it.

## Geometry and placement (concrete numbers)

Head-local frame: x = L-R, y = A-P (+y anterior), z = S-I (+z superior); origin at
skull centre; brain centre at (0, 0, +10). Existing MIRD semi-axes (mm): scalp
(72, 102, 87), skull-out (68, 98, 83), brain (60, 90, 65).

**Tumour** (posterior fossa / 4th ventricle; midline, posterior, inferior):
- head-local centre **(0, −25, −35)**, semi-axes **(18, 20, 18) mm** (≈ Ø 3.8 cm),
  material **`G4_TISSUE_SOFT_ICRP`**. Sits inside the brain ellipsoid.

**Posterior placement** (`RotateX 90°`, world = (x, −z, y) in head-local):
- translation **T = (0, −35, 0)** drops the tumour onto the beam axis.
- beam line ends up at head **midline (L-R = 0), inferior (S-I = −35)**, straight
  along A-P: a suboccipital midline field through occiput → cerebellum → tumour →
  brainstem.
- tumour centre → world (0, 0, −25); **depth ≈ 77 mm** from the posterior entrance
  (A-P scalp semi 102). `BeamAxisHalfExtent` = A-P scalp semi = 102 mm.

**Dose target = the tumour**: box radius ≈ 22 mm, depth window ≈ **57–97 mm**
(brackets the tumour along the beam); beam disk radius ≈ 22 mm.

**Region → world mapping** for the posterior orientation: for a head-local
ellipsoid with semi-axes (sx, sy, sz) at centre (cx, cy, cz),
- world semi-axes = (sx, sz, sy)
- world centre = (cx, −cz − 35, cy)

`HeadRegion` becomes orientation-aware (90° rotations keep regions axis-aligned;
only permute/sign the axes and centre).

## Parameters to confirm (all `StageAConfig` constants, easy to tune)

1. **Tumour material** — `G4_TISSUE_SOFT_ICRP` (default). [confirm]
2. **Tumour size/centre** — (0, −25, −35), semi-axes (18, 20, 18) mm. [confirm]
3. **Target window / disk** — ≈ 57–97 mm depth, ≈ 22 mm radius. [confirm]

## Step-by-step plan

Each step: code → test → plot → (docs) → commit, the established pattern.

### Step 0 — parameters & config
`StageAConfig.hh`: tumour ellipsoid constants (centre, semi-axes, material), the
`headep` geometry string, an orientation constant. No behaviour yet.

### Step 1 — geometry (`BuildHeadEP`)
`DetectorConstruction`:
- generalise head placement to take an **orientation** (rotation) and the on-axis
  **translation**; make `HeadRegion` orientation-aware.
- build scalp/skull/brain **+ tumour** (highest priority, carves the brain); place
  posterior with T = (0, −35, 0); set `BeamAxisHalfExtent` = A-P scalp semi; set
  the target window to bracket the tumour.
- `DetectorMessenger`: register `headep`.

**Test:** extend `MaterialAt` checks — points inside the tumour return tumour
material; the on-axis point at 77 mm depth is tumour. `check_run.py` carve-order
assertion becomes `tumour < brain < skull < scalp`. Build + a 1e4 smoke run to
emit `phantom_regions.csv` + `run_meta.csv`.

### Step 2 — draw the head and the beam
Two paths, both delivered:
- **Python (the committed figure).** Extend `plot_phantom.py` / `plot_mird_head.py`
  to colour the tumour region distinctly and render the posterior case. The head is
  placed so the beam is +z and the tumour is on-axis, so the existing z–x / z–y
  cross-section plotter renders it with no structural change: occiput-to-front head
  section, tumour on the axis, beam column. Output `figures/phantom.png`, checkable
  **against the MRI**. Wire into `make_figures.py`.
- **G4 machinery (interactive sanity check).** A `headep_vis.mac` for the Qt viewer:
  draws the actual volumes in 3D and overlays a few proton trajectories.

Python is the reproducible figure for the docs/paper; G4 Qt is for eyeballing the
placement.

### Step 3 — SOBP field for HeadEP
Run `field_design/sobp.py --from-run` on a `headep` smoke run. It reads
`phantom_regions.csv` + the target window and ray-traces WEPL along +z through
occiput → bone → cerebellum → tumour — expected to need **no code change**; verify
the ray-trace picks up the posterior path and the WEPL window is sane. Add
`headep_sobp.mac`.

### Step 4 — production run + verification + figures
Run `headep_sobp` (1e6, then 1e7). Gate with `check_run.py` (R80 on the tumour
distal edge, plateau uniformity). `make_figures.py`; the `dose_activity` overlay
showing A(x)→0 on the tumour's distal edge.

### Step 5 — docs
Add the HeadEP case (tumour, posterior beam, the drawing) to `ptcrysp_guide.tex` /
`ptcrysp_physics.tex` and `CLAUDE.md`.

## Status

| Step | What | State |
|---|---|---|
| 0 | parameters & config | done |
| 1 | geometry `BuildHeadEP` + messenger | done (smoke run, check_run 10/10) |
| 2 | draw head + beam (Python + G4) | done (phantom.png + headep_vis.mac) |
| 3 | SOBP field (`sobp.py --from-run`) | not started |
| 4 | production run + verification + figures | not started |
| 5 | docs | not started |
