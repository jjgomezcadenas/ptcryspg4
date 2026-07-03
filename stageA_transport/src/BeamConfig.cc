#include "BeamConfig.hh"

#include "BeamMessenger.hh"

#include "Randomize.hh"
#include "G4ios.hh"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

namespace {
// A bad layer table must stop the program at load time: sampling a malformed
// or zero-weight table would deliver the wrong beam while the run is tagged
// sobp (a NaN cumulative distribution sends every primary to the last layer).
void FailLayers(const G4String& path, const std::string& why) {
  G4Exception("BeamConfig::LoadLayers", "StageA.BadLayerTable", FatalException,
              (path + ": " + why).c_str());
}
}  // namespace

BeamConfig::BeamConfig() { fMessenger = new BeamMessenger(this); }

BeamConfig::~BeamConfig() { delete fMessenger; }

void BeamConfig::LoadLayers(const G4String& path) {
  std::ifstream f(path);
  if (!f) {
    FailLayers(path, "cannot open SOBP layer file");
    return;
  }
  fLayerPath = path;  // remembered so the run can copy it into its own dir
  fEnergyMeV.clear();
  fCumWeight.clear();

  std::vector<G4double> w;
  std::string line;
  std::getline(f, line);  // skip header "energy_MeV,weight"
  int lineNo = 1;
  while (std::getline(f, line)) {
    ++lineNo;
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string e, wt;
    std::getline(ss, e, ',');
    std::getline(ss, wt, ',');
    G4double energy = 0., weight = 0.;
    try {
      energy = std::stod(e);
      weight = std::stod(wt);
    } catch (const std::exception&) {
      FailLayers(path, "line " + std::to_string(lineNo) +
                           " is not numeric: '" + line + "'");
    }
    if (energy <= 0.)
      FailLayers(path, "non-positive energy on line " + std::to_string(lineNo));
    if (weight < 0.)
      FailLayers(path, "negative weight on line " + std::to_string(lineNo));
    fEnergyMeV.push_back(energy);
    w.push_back(weight);
  }
  if (fEnergyMeV.empty())
    FailLayers(path, "no layers (empty or header-only file)");

  // Build the normalized cumulative distribution for sampling.
  G4double sum = 0.;
  for (G4double x : w) sum += x;
  if (sum <= 0.) FailLayers(path, "layer weights sum to zero");
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
