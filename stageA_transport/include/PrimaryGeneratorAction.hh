#ifndef STAGEA_PRIMARYGENERATORACTION_HH
#define STAGEA_PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class BeamConfig;

// Proton pencil beam (spec sec 2.1), directed along +z into the phantom face.
// Energy: a single fixed value by default, or — if an SOBP layer table has been
// loaded into BeamConfig — sampled per primary from the layers (Spread-Out
// Bragg Peak, see latex/02_beam_design.tex). Transverse profile is a Gaussian pencil
// (lateral spreading over the target is a later step).
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction(const DetectorConstruction* det, const BeamConfig* beam);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ParticleGun* fGun = nullptr;
  const DetectorConstruction* fDet = nullptr;
  const BeamConfig* fBeam = nullptr;  // shared; SOBP layers if loaded

  G4double fEnergy = 100. * 1.0;   // MeV, single-energy fallback (set in ctor)
  G4double fSigmaXY = 3. * 1.0;    // mm, transverse Gaussian sigma
};

#endif  // STAGEA_PRIMARYGENERATORACTION_HH
