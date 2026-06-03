#ifndef STAGEA_STEPPINGACTION_HH
#define STAGEA_STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"

// Bins each step's energy deposit by depth z into the thread's StageARun,
// separating the primary proton's contribution (parentID==0) from everything
// else, so the depth-dose curve and the secondary contribution can be compared.
class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction() = default;
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;
};

#endif  // STAGEA_STEPPINGACTION_HH
