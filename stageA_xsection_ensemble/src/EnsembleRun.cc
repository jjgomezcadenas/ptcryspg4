#include "EnsembleRun.hh"

EnsembleRun::EnsembleRun(double halfZ_mm, double ebin_width_MeV,
                         double emax_MeV, double zbin_width_mm)
    : fHalfZmm(halfZ_mm) {
  fGrid.Init(ebin_width_MeV, emax_MeV, zbin_width_mm, 2. * halfZ_mm);
  fEdepZTotal.assign(fGrid.DepthBins(), 0.);
  fEdepZPrimary.assign(fGrid.DepthBins(), 0.);
  fEdepZCore.assign(fGrid.DepthBins(), 0.);
}

void EnsembleRun::AddDepthEdep(int bin, double edep_MeV, bool primary,
                               bool core) {
  if (bin < 0 || bin >= DepthBins()) return;
  fEdepZTotal[bin] += edep_MeV;
  if (primary) fEdepZPrimary[bin] += edep_MeV;
  if (core) fEdepZCore[bin] += edep_MeV;
}

void EnsembleRun::Merge(const G4Run* other_run) {
  const auto* other = static_cast<const EnsembleRun*>(other_run);
  fGrid.MergeFrom(other->fGrid);
  for (const auto& [key, count] : other->fNative) fNative[key] += count;
  fSampled.insert(fSampled.end(), other->fSampled.begin(),
                  other->fSampled.end());
  fBank.insert(fBank.end(), other->fBank.begin(), other->fBank.end());
  for (const auto& [event, xyz] : other->fAnnihilations)
    fAnnihilations[event] = xyz;
  fEdepTotal += other->fEdepTotal;
  fTargetEdep += other->fTargetEdep;
  for (int i = 0; i < DepthBins(); ++i) {
    fEdepZTotal[i] += other->fEdepZTotal[i];
    fEdepZPrimary[i] += other->fEdepZPrimary[i];
    fEdepZCore[i] += other->fEdepZCore[i];
  }
  G4Run::Merge(other_run);
}
