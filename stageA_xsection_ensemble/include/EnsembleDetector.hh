#ifndef ENSEMBLE_DETECTOR_HH
#define ENSEMBLE_DETECTOR_HH

#include "EnsembleConfig.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;

/// Phantom geometries for the exposure calculation: the homogeneous cylinder
/// and the homogeneous head placements (lateral uniform_head, posterior
/// uniform_headep), replicating the stageA_transport shapes and constants.
class EnsembleDetector : public G4VUserDetectorConstruction {
 public:
  explicit EnsembleDetector(const ensemble::EnsembleCli& cli);
  G4VPhysicalVolume* Construct() override;

  G4LogicalVolume* PhantomLogical() const { return fPhantomLV; }
  /// Beam-axis half-extent of the phantom (depth window = 2x this).
  G4double BeamAxisHalfExtent() const { return fHalfZ; }
  G4double PhantomMass() const;
  G4double TargetMass() const;
  /// Target-box window in world z (from the configured entrance-face depths).
  G4double TargetProxZ() const { return -fHalfZ + fTargetProx; }
  G4double TargetDistZ() const { return -fHalfZ + fTargetDist; }
  G4double TargetRadius() const { return fTargetRadius; }
  /// True when the on-axis point (0,0,z) lies inside the phantom medium.
  bool ContainsOnAxis(G4double z) const;

 private:
  const ensemble::EnsembleCli& fCli;
  G4LogicalVolume* fPhantomLV = nullptr;
  G4double fHalfZ = 0.;
  G4double fTargetRadius = 0., fTargetProx = 0., fTargetDist = 0.;
};

#endif  // ENSEMBLE_DETECTOR_HH
