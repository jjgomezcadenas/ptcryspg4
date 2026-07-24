#pragma once

#include "G4Run.hh"
#include "globals.hh"

class ThinTargetRun : public G4Run {
 public:
  void CountIncident() { ++incident_; }
  void CountInelastic(G4double energy);
  void CountResidual(G4int z, G4int a, G4double energy);
  void AddContinuousLoss(G4double loss) { continuous_loss_ += loss; }
  void Merge(const G4Run* other) override;

  G4long Incident() const { return incident_; }
  G4long Inelastic() const { return inelastic_; }
  G4long C11() const { return c11_; }
  G4long N13() const { return n13_; }
  G4long O15() const { return o15_; }
  G4double InelasticEnergySum() const { return inelastic_energy_sum_; }
  G4double C11EnergySum() const { return c11_energy_sum_; }
  G4double N13EnergySum() const { return n13_energy_sum_; }
  G4double O15EnergySum() const { return o15_energy_sum_; }
  G4double ContinuousLoss() const { return continuous_loss_; }

 private:
  G4long incident_ = 0;
  G4long inelastic_ = 0;
  G4long c11_ = 0;
  G4long n13_ = 0;
  G4long o15_ = 0;
  G4double inelastic_energy_sum_ = 0.0;
  G4double c11_energy_sum_ = 0.0;
  G4double n13_energy_sum_ = 0.0;
  G4double o15_energy_sum_ = 0.0;
  G4double continuous_loss_ = 0.0;
};
