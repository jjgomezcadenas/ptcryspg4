# File formats

Units: position mm, energy keV, time ns, dose Gy.
Isotope codes are in isotopes.csv.

## Coordinate frame

All positions (emitters.csv, the phantom, the target box) share one frame:
the phantom is a cylinder **centred at the origin** with its axis along **+z**,
the direction the proton beam travels.

- Phantom spans `z` in [-L/2, +L/2] and radius `r = sqrt(x^2+y^2)` in [0, R],
  with L = `phantom_length_mm` and R = `phantom_diameter_mm`/2 (run_meta.csv).
- The beam enters at the `z = -L/2` face. Depth into the phantom is
  `d = z + L/2`, i.e. `z = d - L/2` (so the target box depths in run_meta.csv,
  measured from the entrance face, map to `z = depth - L/2`).
- emitters.csv prod/anh coordinates are in this same frame — place the phantom
  for gamma transport with the origin at its centre and you are co-registered.

## emitters.csv — one row per positron emitter
    event_id                       proton that produced it
    isotope_id                     isotope code (see isotopes.csv)
    prod_x_mm prod_y_mm prod_z_mm  where the isotope was produced
    anh_x_mm  anh_y_mm  anh_z_mm   where the positron annihilated — the source seen by the detector

The positron range per emitter is `anh - prod`. Positrons that would leave the
phantom into air are stopped at the surface, so a small fraction (~0.6%) have
their `anh` pinned to the phantom boundary.

## run_meta.csv — one row of run-level numbers
    n_protons              protons simulated
    beam_energy_MeV        nominal proton energy; for a SOBP run this is the
                           fallback single energy — the real spectrum is in
                           sobp_layers.csv
    beam_sigma_mm          Gaussian beam sigma at entrance (pencil mode)
    phantom_material       Geant4 NIST material name (see phantom_material_*.csv)
    phantom_diameter_mm    phantom cylinder diameter
    phantom_length_mm      phantom cylinder length (along z)
    phantom_mass_g         phantom mass
    edep_total_MeV         total energy deposited in the phantom
    dose_total_Gy          whole-phantom dose for n_protons
    target_dose_Gy         dose in the target box for n_protons
    target_mass_g          target-box mass
    target_radius_mm       target box (cylinder) radius
    target_prox_depth_mm   near face of the target box, depth from entrance
    target_dist_depth_mm   far face of the target box, depth from entrance
    Np_per_Gy              protons for 1 Gy (= n_protons / target_dose_Gy)
    geant4_version, physics_list, random_seed

## sampling_budget_<scenario>.csv — measured decays per isotope
    isotope_id    isotope code
    N_expected    expected decays measured for 1 Gy at this scenario's timing

Each *_meta.csv next to it gives the dose and timing (t_irr, t_del, t_meas) used.
The scenarios differ only in the acquisition timing (see the README legend);
the spatial source (emitters.csv) is shared by all of them.

## isotopes.csv — isotope code table
    isotope_id     code used in emitters.csv and the budgets
    name           e.g. O15
    half_life_s    physical half-life
    endpoint_MeV   beta+ endpoint energy (governs positron range)
    prompt_gamma   1 if the isotope emits a prompt de-excitation gamma in
                   coincidence with the decay (relevant to PET backgrounds)

## sobp_layers.csv — the beam definition
The proton energy layers (energy and weight) that make up the spread-out Bragg
peak. Input to the run; kept so the field is fully documented.

## phantom_material_<name>.csv — the medium for gamma transport
The source is detector-independent but not medium-independent: a PET sim must
propagate the 511 keV annihilation photons through the phantom and needs the
attenuation coefficient for reconstruction. This file is the phantom's elemental
composition (one row per element: element, Z, A_g_mol, mass_fraction).
    phantom_material_<name>_meta.csv  density_g_cm3, mean_excitation_eV,
        mu_rho_cm2_g and mu_cm_inv at energy_keV=511 (annihilation), mean_free_path_cm
<name> is the run_meta.csv phantom_material, lower-cased (e.g. g4_brain_icrp).
Valid only for this homogeneous phantom; a heterogeneous phantom would instead
need a voxel attenuation map.
