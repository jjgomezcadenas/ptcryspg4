#include "SteppingAction.hh"

#include "StageARun.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

void SteppingAction::UserSteppingAction(const G4Step* step) {
  G4Track* track = step->GetTrack();

  // --- depth-binned dose (only steps that deposit energy) -------------------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) {
    // Distribute the step's edep across the z-bins it spans (pre -> post).
    const G4double z1 = step->GetPreStepPoint()->GetPosition().z();
    const G4double z2 = step->GetPostStepPoint()->GetPosition().z();
    const bool primary = (track->GetParentID() == 0);

    auto* run = static_cast<StageARun*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEdepAlongStep(z1 / mm, z2 / mm, edep / MeV, primary);
  }

  // --- kill positrons that escape the phantom into air ----------------------
  // A positron born in the first few mm at the entrance face can travel back
  // out into air, where its range is huge, and annihilate hundreds of mm away
  // (0.64% of emitters). Killing it the moment it leaves the phantom pins the
  // capture: the TrackingAction records ANH at this boundary point instead
  // (keep-at-surface). The annihilation photons are not generated, but Stage A
  // does not use them -- Stage B regenerates the pair from anh_xyz.
  if (track->GetDefinition() == G4Positron::Definition()) {
    const auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    if (postVol == nullptr || postVol->GetName() != "Phantom") {
      track->SetTrackStatus(fStopAndKill);
    }
  }
}
