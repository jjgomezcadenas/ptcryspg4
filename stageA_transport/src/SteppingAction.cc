#include "SteppingAction.hh"

#include "StageARun.hh"
#include "StageAConfig.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

void SteppingAction::UserSteppingAction(const G4Step* step) {
  G4Track* track = step->GetTrack();

  // --- depth-binned dose + target-box dose (steps that deposit energy) ------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) {
    const auto& pre = step->GetPreStepPoint()->GetPosition();
    const auto& post = step->GetPostStepPoint()->GetPosition();
    const bool primary = (track->GetParentID() == 0);

    auto* run = static_cast<StageARun*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());

    // Transverse radius sqrt(x^2+y^2), shared by the on-axis core and the target.
    const G4double rmid = 0.5 * (pre.perp() + post.perp());

    // Depth-dose: distribute edep across the z-bins the step spans; steps within
    // the thin on-axis core also feed the central-axis profile.
    const bool inCore = rmid <= stageA::kCoreRadiusMM * mm;
    run->AddEdepAlongStep(pre.z() / mm, post.z() / mm, edep / MeV, primary,
                          inCore);

    // Target-box dose: add edep if the step midpoint is inside the target
    // (cylinder: target z-range and radius from DetectorConstruction).
    const G4double zmid = 0.5 * (pre.z() + post.z());
    if (zmid >= fDet->TargetProxZ() && zmid <= fDet->TargetDistZ() &&
        rmid <= fDet->TargetRadius()) {
      run->AddTargetEdep(edep);
    }
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
