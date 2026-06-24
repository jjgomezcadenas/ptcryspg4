# File formats

Units: position mm, energy keV, time ns, dose Gy.
Isotope codes are in isotopes.csv.

## emitters.csv — one row per positron emitter
    event_id                       proton that produced it
    isotope_id                     isotope code (see isotopes.csv)
    prod_x_mm prod_y_mm prod_z_mm  where the isotope was produced
    anh_x_mm  anh_y_mm  anh_z_mm   where the positron annihilated — the source seen by the detector

## run_meta.csv — one row of run-level numbers
    n_protons          protons simulated
    beam_energy_MeV    proton energy
    phantom_material   phantom material
    target_dose_Gy     dose in the target box for n_protons
    Np_per_Gy          protons for 1 Gy (= n_protons / target_dose_Gy)
    geant4_version, physics_list, random_seed

## sampling_budget_<scenario>.csv — measured decays per isotope
    isotope_id    isotope code
    N_expected    expected decays measured for 1 Gy at this scenario's timing

Each *_meta.csv next to it gives the dose and timing (t_irr, t_del, t_meas) used.

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
