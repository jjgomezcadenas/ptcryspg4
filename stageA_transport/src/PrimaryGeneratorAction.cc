#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "BeamConfig.hh"
#include "StageAConfig.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "Randomize.hh"

// Configure the gun once: a single proton along +z. The energy is set per event
// in GeneratePrimaries (it may vary, in SOBP mode).
PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction* det,
                                               const BeamConfig* beam)
    : fDet(det), fBeam(beam) {
  fEnergy = stageA::kBeamEnergyMeV * MeV;
  fSigmaXY = stageA::kBeamSigmaMM * mm;

  fGun = new G4ParticleGun(1);  // "1" = one primary particle per vertex
  auto* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGun->SetParticleDefinition(proton);
  fGun->SetParticleMomentumDirection({0., 0., 1.});
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fGun; }

// Called by the run manager once at the start of every event: pick this
// proton's energy (SOBP layer sample or the single value) and entry point.
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // Energy: sample an SOBP layer if a table is loaded, else the fixed value.
  const G4double energy =
      (fBeam && fBeam->SobpEnabled()) ? fBeam->SampleEnergyMeV() * MeV : fEnergy;
  fGun->SetParticleEnergy(energy);

  // Start just upstream of the phantom entrance face, on the -z side.
  const G4double z0 = -fDet->PhantomHalfLength() - 1. * mm;
  const G4double x = G4RandGauss::shoot(0., fSigmaXY);
  const G4double y = G4RandGauss::shoot(0., fSigmaXY);
  fGun->SetParticlePosition({x, y, z0});
  fGun->GeneratePrimaryVertex(event);
}
