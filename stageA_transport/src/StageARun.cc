#include "StageARun.hh"

#include "StageAConfig.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"

#include <algorithm>
#include <cmath>
#include <string>

// Depth binning along z (phantom centred at origin, axis along z): the window
// is the phantom's beam-axis extent, passed in by RunAction::GenerateRun.
StageARun::StageARun(G4double halfZ_mm)
    : fZMinMM(-halfZ_mm), fBinWMM(2. * halfZ_mm / kNZBins) {}

// Called by the run manager for every event (after EndOfEventAction). We use it
// to pull the per-event scorer result; emitter rows and depth-dose are filled
// separately by EventAction / SteppingAction.
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
                                 bool primary, bool inCore) {
  const G4double zlo = std::min(z1, z2);
  const G4double zhi = std::max(z1, z2);
  const G4double len = zhi - zlo;

  // Degenerate (zero-length) step: assign to its single bin.
  if (len < 1e-9) {
    const int bin = static_cast<int>((zlo - fZMinMM) / fBinWMM);
    if (bin >= 0 && bin < kNZBins) {
      fEdepZTotal[bin] += edep;
      if (primary) fEdepZPrimary[bin] += edep;
      if (inCore) fEdepZCore[bin] += edep;
    }
    return;
  }

  // Spread edep over the spanned bins, in proportion to the overlap length.
  int binLo = static_cast<int>(std::floor((zlo - fZMinMM) / fBinWMM));
  int binHi = static_cast<int>(std::floor((zhi - fZMinMM) / fBinWMM));
  binLo = std::max(binLo, 0);
  binHi = std::min(binHi, kNZBins - 1);
  for (int b = binLo; b <= binHi; ++b) {
    const G4double bz0 = fZMinMM + b * fBinWMM;
    const G4double overlap = std::min(zhi, bz0 + fBinWMM) - std::max(zlo, bz0);
    if (overlap <= 0.) continue;
    const G4double e = edep * (overlap / len);
    fEdepZTotal[b] += e;
    if (primary) fEdepZPrimary[b] += e;
    if (inCore) fEdepZCore[b] += e;
  }
}

// Called on the master once per worker run at end of run: fold that worker's
// accumulators into the master's. This is the *only* cross-thread combine, so
// no locking is needed anywhere in the per-event hot path.
void StageARun::Merge(const G4Run* aRun) {
  const auto* local = static_cast<const StageARun*>(aRun);
  fEmitters.append(local->fEmitters);
  fEdep += local->fEdep;
  fTargetEdep += local->fTargetEdep;
  for (int i = 0; i < kNZBins; ++i) {
    fEdepZTotal[i] += local->fEdepZTotal[i];
    fEdepZPrimary[i] += local->fEdepZPrimary[i];
    fEdepZCore[i] += local->fEdepZCore[i];
  }
  G4Run::Merge(aRun);
}
