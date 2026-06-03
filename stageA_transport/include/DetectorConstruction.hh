#ifndef STAGEA_DETECTORCONSTRUCTION_HH
#define STAGEA_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;

// Stage A geometry (spec sec 2.2): a homogeneous right-circular PMMA cylinder,
// axis along +z, acting both as the production target and as the passive
// 511 keV attenuator. Dimensions come from StageAConfig.
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction() = default;
  ~DetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;  // dose scorer on the phantom

  // Geometry accessors (used by the primary generator to place the beam).
  G4double PhantomHalfLength() const { return fHalfZ; }
  G4double PhantomRadius() const { return fRadius; }
  // Phantom mass (for the dose normalization); valid after Construct().
  G4double PhantomMass() const;

 private:
  G4double fRadius = 0.;  // set in Construct()
  G4double fHalfZ = 0.;
  G4LogicalVolume* fPhantomLV = nullptr;
};

#endif  // STAGEA_DETECTORCONSTRUCTION_HH
