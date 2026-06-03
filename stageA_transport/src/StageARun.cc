#include "StageARun.hh"

#include "StageAConfig.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"

#include <algorithm>
#include <cmath>
#include <string>

namespace {
// Depth binning along z (phantom centred at origin, axis along z).
constexpr G4double kZMin = -0.5 * stageA::kPhantomLengthMM;  // mm
constexpr G4double kBinW = stageA::kPhantomLengthMM / StageARun::kNZBins;  // mm
}  // namespace

StageARun::StageARun() = default;

void StageARun::RecordEvent(const G4Event* event) {
  // Sum the phantom energy-deposit scorer for this event.
  if (fCollID < 0) {
    const std::string name =
        std::string(stageA::kScorerMFD) + "/" + stageA::kScorerEdep;
    fCollID = G4SDManager::GetSDMpointer()->GetCollectionID(name);
  }

  auto* hce = event->GetHCofThisEvent();
  if (hce && fCollID >= 0) {
    auto* hits = static_cast<G4THitsMap<G4double>*>(hce->GetHC(fCollID));
    if (hits) {
      for (const auto& entry : *hits->GetMap()) {
        fEdep += *(entry.second);
      }
    }
  }

  G4Run::RecordEvent(event);  // keep G4's event counter correct
}

void StageARun::AddEdepAlongStep(G4double z1, G4double z2, G4double edep,
                                 bool primary) {
  const G4double zlo = std::min(z1, z2);
  const G4double zhi = std::max(z1, z2);
  const G4double len = zhi - zlo;

  // Degenerate (zero-length) step: assign to its single bin.
  if (len < 1e-9) {
    const int bin = static_cast<int>((zlo - kZMin) / kBinW);
    if (bin >= 0 && bin < kNZBins) {
      fEdepZTotal[bin] += edep;
      if (primary) fEdepZPrimary[bin] += edep;
    }
    return;
  }

  // Spread edep over the spanned bins, in proportion to the overlap length.
  int binLo = static_cast<int>(std::floor((zlo - kZMin) / kBinW));
  int binHi = static_cast<int>(std::floor((zhi - kZMin) / kBinW));
  binLo = std::max(binLo, 0);
  binHi = std::min(binHi, kNZBins - 1);
  for (int b = binLo; b <= binHi; ++b) {
    const G4double bz0 = kZMin + b * kBinW;
    const G4double overlap = std::min(zhi, bz0 + kBinW) - std::max(zlo, bz0);
    if (overlap <= 0.) continue;
    const G4double e = edep * (overlap / len);
    fEdepZTotal[b] += e;
    if (primary) fEdepZPrimary[b] += e;
  }
}

void StageARun::Merge(const G4Run* aRun) {
  const auto* local = static_cast<const StageARun*>(aRun);
  fEmitters.append(local->fEmitters);
  fEdep += local->fEdep;
  for (int i = 0; i < kNZBins; ++i) {
    fEdepZTotal[i] += local->fEdepZTotal[i];
    fEdepZPrimary[i] += local->fEdepZPrimary[i];
  }
  G4Run::Merge(aRun);
}
