#include "ThinTargetRun.hh"

void ThinTargetRun::CountInelastic(G4double energy) {
  ++inelastic_;
  inelastic_energy_sum_ += energy;
}

void ThinTargetRun::CountResidual(G4int z, G4int a, G4double energy) {
  if (z == 6 && a == 11) { ++c11_; c11_energy_sum_ += energy; }
  if (z == 7 && a == 13) { ++n13_; n13_energy_sum_ += energy; }
  if (z == 8 && a == 15) { ++o15_; o15_energy_sum_ += energy; }
}

void ThinTargetRun::Merge(const G4Run* other_run) {
  const auto* other = static_cast<const ThinTargetRun*>(other_run);
  incident_ += other->incident_;
  inelastic_ += other->inelastic_;
  c11_ += other->c11_;
  n13_ += other->n13_;
  o15_ += other->o15_;
  inelastic_energy_sum_ += other->inelastic_energy_sum_;
  c11_energy_sum_ += other->c11_energy_sum_;
  n13_energy_sum_ += other->n13_energy_sum_;
  o15_energy_sum_ += other->o15_energy_sum_;
  continuous_loss_ += other->continuous_loss_;
  G4Run::Merge(other_run);
}
