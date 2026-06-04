#ifndef STAGEA_STEPPINGACTION_HH
#define STAGEA_STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"

class DetectorConstruction;

// Per-step bookkeeping: bins each step's energy deposit by depth z (the
// depth-dose), accumulates the part that falls inside the target box (the
// transparent dose-scoring region, bounds from DetectorConstruction), and kills
// positrons that escape the phantom into air.
class SteppingAction : public G4UserSteppingAction {
 public:
  explicit SteppingAction(const DetectorConstruction* det) : fDet(det) {}
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  const DetectorConstruction* fDet = nullptr;
};

#endif  // STAGEA_STEPPINGACTION_HH
