#ifndef STAGEA_DETECTORCONSTRUCTION_HH
#define STAGEA_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class DetectorMessenger;

// Stage A geometry (spec sec 2.2): a homogeneous right-circular cylinder, axis
// along +z, that serves both as the β⁺ production target and as the passive
// 511 keV attenuator. Geant4 calls the two virtuals below automatically: it is
// the mandatory "what does the apparatus look like" user class.
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction();             // builds the messenger; sets default material
  ~DetectorConstruction() override;   // deletes the messenger

  // Build the volume tree (world + phantom) and return the world physical
  // volume. Called once by the run manager at /run/initialize.
  G4VPhysicalVolume* Construct() override;

  // Attach the energy-deposit scorer to the phantom. Called after Construct(),
  // once per worker thread in MT, so each thread owns its sensitive detector.
  void ConstructSDandField() override;

  // Phantom material (NIST name, e.g. G4_TISSUE_SOFT_ICRP). Settable from a
  // macro before /run/initialize via DetectorMessenger; default from StageAConfig.
  void SetMaterialName(const G4String& name) { fMaterialName = name; }
  const G4String& MaterialName() const { return fMaterialName; }

  // Beam placement helpers, read by the primary generator.
  G4double PhantomHalfLength() const { return fHalfZ; }
  G4double PhantomRadius() const { return fRadius; }
  // Phantom mass (ρ·V), for the dose normalization. Valid only after Construct().
  G4double PhantomMass() const;

 private:
  G4double fRadius = 0.;   // cylinder radius, set in Construct()
  G4double fHalfZ = 0.;    // cylinder half-length (G4Tubs uses half-z)
  G4LogicalVolume* fPhantomLV = nullptr;  // kept for the scorer + mass query
  G4String fMaterialName;                 // NIST name; default set in the ctor
  DetectorMessenger* fMessenger = nullptr;
};

#endif  // STAGEA_DETECTORCONSTRUCTION_HH
