#include "BeamConfig.hh"

#include "BeamMessenger.hh"

#include "Randomize.hh"
#include "G4ios.hh"

#include <algorithm>
#include <fstream>
#include <sstream>

BeamConfig::BeamConfig() { fMessenger = new BeamMessenger(this); }

BeamConfig::~BeamConfig() { delete fMessenger; }

void BeamConfig::LoadLayers(const G4String& path) {
  std::ifstream f(path);
  if (!f) {
    G4cerr << "[Beam] ERROR: cannot open SOBP layer file " << path << G4endl;
    return;
  }
  fEnergyMeV.clear();
  fCumWeight.clear();

  std::vector<G4double> w;
  std::string line;
  std::getline(f, line);  // skip header "energy_MeV,weight"
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string e, wt;
    std::getline(ss, e, ',');
    std::getline(ss, wt, ',');
    fEnergyMeV.push_back(std::stod(e));
    w.push_back(std::stod(wt));
  }

  // Build the normalized cumulative distribution for sampling.
  G4double sum = 0.;
  for (G4double x : w) sum += x;
  G4double cum = 0.;
  for (G4double x : w) {
    cum += x / sum;
    fCumWeight.push_back(cum);
  }

  G4cout << "[Beam] loaded " << fEnergyMeV.size() << " SOBP layers from " << path
         << "  (" << fEnergyMeV.front() << "-" << fEnergyMeV.back() << " MeV)"
         << G4endl;
}

G4double BeamConfig::SampleEnergyMeV() const {
  // Inverse-CDF sampling: draw u, take the first layer whose cumulative weight
  // reaches it.
  const G4double u = G4UniformRand();
  auto it = std::lower_bound(fCumWeight.begin(), fCumWeight.end(), u);
  const std::size_t i = (it == fCumWeight.end())
                            ? fCumWeight.size() - 1
                            : static_cast<std::size_t>(it - fCumWeight.begin());
  return fEnergyMeV[i];
}
