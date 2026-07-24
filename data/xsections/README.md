# Proton-production cross sections

This directory stores the public cross-section inputs used to quantify the
positron-emitter production uncertainty. It has three data layers:

- `raw/`: immutable files as distributed by EXFOR, TENDL and JENDL;
- `normalized/`: one-to-one conversions into the common point schema;
- `models/`: documented interpolation-ready curves used by the reweighting
  calculation.

Covariance tables derived from an evaluated source live in `covariance/`.
`sources.csv` identifies every upstream distribution or file, its retrieval
location, revision and checksum. `SCHEMA.md` defines the normalized and model
tables.

## Current raw snapshot

- EXFOR Plot File commit `57052199247e2f64862a320e3930d174318f5e6c`, dated
  2026-06-29: 67 CX4 experimental series.
  - 30 series for C-12(p,x)C-11;
  - 12 series for O-16(p,x)C-11;
  - 12 series for O-16(p,x)O-15;
  - 13 series for N-14(p,x)N-13.
- Incident-proton files served by the TENDL-2023 distribution for C-12, N-14
  and O-16. Their internal headers identify the adopted LANL ENDF/B-VII.1
  evaluations dated 1996--1997.
- JENDL-4.0/HE residual-production tables for proton interactions on C-12,
  N-14 and O-16.

The Geant4 reference curves are generated locally with the thin-target runs and
will enter `models/` after validation. They are not imported source files.
