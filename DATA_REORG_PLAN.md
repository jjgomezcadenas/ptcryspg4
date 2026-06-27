# Plan — data organization rework (branch `data-reorg`)

Off `main` after Phase 2 merged. Kills the shared-`data/` coupling: each run gets
its own directory, and a snapshot is tied to a *specific run* (identity, file
set, and figures derived from that run's `run_meta`) instead of a hand-typed name
plus whatever happens to be in `data/`.

This is a living plan: written before compaction so no scope is lost. Implement
with the confirmed defaults below; do not start until the user gives the OK.

## The problem (today)

A snapshot is defined by **(a hand-typed name) + (the transient contents of the
shared `data/` scratch dir)** — nothing binds the two. Consequences:
- **Name decoupled from run.** `data/` is overwritten by every run; the name is
  the only thing distinguishing cases and nothing checks it (you can freeze a
  `uniform_head` run as `head_mird_1e7`).
- **Hardcoded file/figure assumptions.** `snapshot_scenario.py` copies a fixed
  file list + a fixed (cylinder-oriented) figure set, so a head *pencil* run with
  no `sampling_budget`/`sobp_layers` silently bundles wrong/old figures.
- **No coexistence.** All runs funnel through one `data/`; cases can't sit side
  by side; reproducing one means re-running.

## Confirmed decisions (defaults)

1. **Run tag = `<geometry>_<beam>_<N>`, AUTO-DERIVED** (e.g. `mird_head_pencil_1e5`,
   `cylinder_sobp_1e7`). `RunAction` builds it from the actual run — `geometry`
   (`fDet`), `beam` (pencil/sobp — pass `BeamConfig` into `RunAction`), `N`
   (n_protons). Robust: the tag can't disagree with what ran.
2. **Frozen scenarios in `ptcrysp-scenarios` stay unchanged** — this only changes
   the transient `data/` side.
3. **Add `make_figures.py`** as the geometry-aware figure dispatcher.

## New layout

```
data/runs/<run_tag>/        one self-contained run
    emitters.csv  run_meta.csv  phantom_regions.csv  depth_dose.csv
    sampling_budget_*.csv (+ _meta)      (added by budget.py for this run)
    figures/                              (geometry-appropriate control plots)
```
Runs no longer clobber each other; same config re-run overwrites only itself.
`sobp_layers.csv` (a beam *input* made by `field_design/sobp.py`) is read by the
gun and copied into the run dir if used.

## Components

- **Per-run output dir.** Macros set the base (`/stageA/output/dir data/runs`);
  `RunAction` appends `<run_tag>` and writes all Stage-A outputs there. The tag is
  computed at end-of-run (n_protons known). Plumb `BeamConfig` (or a pencil/sobp
  label) into `RunAction` via `ActionInitialization` for the beam part of the tag.
- **Snapshot.** `snapshot_scenario.py <run_tag>` (scenario name defaults to the
  tag; `--name` to override): read `data/runs/<run_tag>/run_meta.csv` for
  identity; copy the files that **actually exist** (glob, not a hardcoded list);
  copy `figures/` from the run dir; keep the global bits (doc PDFs, `isotopes.csv`,
  `SCHEMA.md`, generated README).
- **Geometry-aware figures.** `make_figures.py <run_dir>` reads
  `run_meta.geometry` and calls the right plotters into `<run_dir>/figures/`:
  cylinder → `validate_transport` + `plot_sobp`; head → `plot_mird_head`. (Plot
  scripts already take a data dir, so this is mostly dispatch.)
- **Budget / analysis take a `<run_dir>`.** `budget.py`, `budget_gen.py`,
  `activity_plot.py` operate on one run.

## Files touched

- `stageA_transport`: `RunAction.{hh,cc}` (tag + subdir creation),
  `ActionInitialization.{hh,cc}` (pass beam to RunAction), `macros/*.mac` (base
  output dir).
- `tools/snapshot_scenario.py` (run-dir arg + glob + figures from run dir).
- `decay_sampling/budget.py`, `budget_gen.py`, `activity_plot.py` (run-dir arg).
- `analysis_transport/`: new `make_figures.py`; plot scripts already take a dir.
- `.gitignore` (`data/runs/`), `CLAUDE.md` (run/snapshot instructions).

## Acceptance

- Two different cases run back-to-back both present under `data/runs/`, neither
  clobbered.
- `snapshot_scenario.py mird_head_pencil_1e5` freezes that exact run with its head
  figure and only the files it has.
- Re-running the same config is idempotent.

## Status

Not started — scope written, awaiting the OK to implement (after compaction).
