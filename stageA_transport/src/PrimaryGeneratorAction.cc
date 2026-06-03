#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "StageAConfig.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "Randomize.hh"

// Configure the gun once: a single proton of fixed energy along +z. Only the
// per-event transverse position changes (set in GeneratePrimaries).
PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction* det)
    : fDet(det) {
  fEnergy = stageA::kBeamEnergyMeV * MeV;
  fSigmaXY = stageA::kBeamSigmaMM * mm;

  fGun = new G4ParticleGun(1);  // "1" = one primary particle per vertex
  auto* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGun->SetParticleDefinition(proton);
  fGun->SetParticleEnergy(fEnergy);
  fGun->SetParticleMomentumDirection({0., 0., 1.});
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fGun; }

// Called by the run manager once at the start of every event: sample this
// proton's entry point and hand the vertex to the event.
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // Start just upstream of the phantom entrance face, on the -z side.
  const G4double z0 = -fDet->PhantomHalfLength() - 1. * mm;
  const G4double x = G4RandGauss::shoot(0., fSigmaXY);
  const G4double y = G4RandGauss::shoot(0., fSigmaXY);
  fGun->SetParticlePosition({x, y, z0});
  fGun->GeneratePrimaryVertex(event);
}
