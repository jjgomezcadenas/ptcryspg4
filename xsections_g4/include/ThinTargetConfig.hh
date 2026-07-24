#pragma once

#include "globals.hh"

struct ThinTargetConfig {
  G4String target = "O16";
  G4double energy_MeV = 100.0;
  G4double areal_mg_cm2 = 5.0;
  G4long protons = 10000;
  G4long seed = 12345;
  G4int threads = 1;
  G4String physics_list = "QGSP_BIC_HP";
  G4String output = "thin_target_counts.csv";
};
