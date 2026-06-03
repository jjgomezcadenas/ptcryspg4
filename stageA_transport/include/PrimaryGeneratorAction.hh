#ifndef STAGEA_PRIMARYGENERATORACTION_HH
#define STAGEA_PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

// Proton pencil beam (spec sec 2.1): mono-energetic, baseline 100 MeV, Gaussian
// transverse profile, directed along +z, entering the phantom front face.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  explicit PrimaryGeneratorAction(const DetectorConstruction* det);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ParticleGun* fGun = nullptr;
  const DetectorConstruction* fDet = nullptr;

  G4double fEnergy = 100. * 1.0;   // MeV, set in ctor
  G4double fSigmaXY = 3. * 1.0;    // mm, transverse Gaussian sigma
};

#endif  // STAGEA_PRIMARYGENERATORACTION_HH
