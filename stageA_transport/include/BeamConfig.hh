#ifndef STAGEA_BEAMCONFIG_HH
#define STAGEA_BEAMCONFIG_HH

#include "globals.hh"
#include <vector>

class BeamMessenger;

// Single shared beam configuration. One instance is created on the master (in
// main) and passed to every worker gun, so the SOBP energy-layer table is
// loaded once (by a UI command on the master) and read concurrently by the
// worker guns. Read-only after loading, so the concurrent reads are safe.
// Until a layer file is loaded the gun uses a single fixed energy.
class BeamConfig {
 public:
  BeamConfig();
  ~BeamConfig();

  // Load an SOBP layer table (CSV: energy_MeV,weight) and enable SOBP mode.
  void LoadLayers(const G4String& path);

  bool SobpEnabled() const { return !fEnergyMeV.empty(); }

  // Sample one layer's energy [MeV] in proportion to its weight, using the
  // calling worker thread's RNG (stochastic delivery, see docs/sobp.tex).
  G4double SampleEnergyMeV() const;

 private:
  std::vector<G4double> fEnergyMeV;
  std::vector<G4double> fCumWeight;  // cumulative weights, normalized to 1
  BeamMessenger* fMessenger = nullptr;
};

#endif  // STAGEA_BEAMCONFIG_HH
