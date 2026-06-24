#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "BeamConfig.hh"
#include "StageAConfig.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include "Randomize.hh"

#include <cmath>

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

  // Transverse start point: a uniform disk over the target cross-section (the
  // realistic field, r = sqrt(u)·R for uniform area sampling), or a Gaussian
  // pencil if no disk radius is set.
  G4double x, y;
  if (fBeam && fBeam->DiskMode()) {
    const G4double r = fBeam->DiskRadius() * std::sqrt(G4UniformRand());
    const G4double phi = twopi * G4UniformRand();
    x = r * std::cos(phi);
    y = r * std::sin(phi);
  } else {
    x = G4RandGauss::shoot(0., fSigmaXY);
    y = G4RandGauss::shoot(0., fSigmaXY);
  }

  // Start just upstream of the phantom entrance face, on the -z side. The
  // half-extent is geometry-aware (cylinder half-length, or the head's L-R
  // semi-axis), so the pencil always begins in air ahead of the phantom.
  const G4double z0 = -fDet->BeamAxisHalfExtent() - 1. * mm;
  fGun->SetParticlePosition({x, y, z0});
  fGun->GeneratePrimaryVertex(event);
}
