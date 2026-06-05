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
