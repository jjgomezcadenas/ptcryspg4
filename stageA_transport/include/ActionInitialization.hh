#ifndef STAGEA_ACTIONINITIALIZATION_HH
#define STAGEA_ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;
class BeamConfig;

// Factory for the user action classes. Geant4 calls BuildForMaster() once on the
// master and Build() on every worker thread, so each thread gets its own action
// objects. Holds the shared (non-owned) DetectorConstruction and BeamConfig to
// hand to those actions.
class ActionInitialization : public G4VUserActionInitialization {
 public:
  ActionInitialization(const DetectorConstruction* det, const BeamConfig* beam)
      : fDet(det), fBeam(beam) {}
  ~ActionInitialization() override = default;

  void Build() const override;           // per-worker thread
  void BuildForMaster() const override;  // master thread (MT)

 private:
  const DetectorConstruction* fDet = nullptr;
  const BeamConfig* fBeam = nullptr;
};

#endif  // STAGEA_ACTIONINITIALIZATION_HH
