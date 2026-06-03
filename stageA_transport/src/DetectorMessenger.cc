#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det) : fDet(det) {
  // "/stageA/" is the shared command root. DetectorConstruction is built before
  // the run actions, so this messenger creates the root; OutputMessenger then
  // only adds "/stageA/output/".
  fDirStageA = new G4UIdirectory("/stageA/");
  fDirStageA->SetGuidance("Stage-A proton-transport controls.");

  fDirPhantom = new G4UIdirectory("/stageA/phantom/");
  fDirPhantom->SetGuidance("Phantom controls.");

  fMatCmd = new G4UIcmdWithAString("/stageA/phantom/material", this);
  fMatCmd->SetGuidance("Phantom material by NIST name, "
                       "e.g. G4_TISSUE_SOFT_ICRP, G4_PLEXIGLASS, G4_WATER.");
  fMatCmd->SetParameterName("name", false);
  // Geometry is constructed at /run/initialize, so only allow this beforehand.
  fMatCmd->AvailableForStates(G4State_PreInit);
}

DetectorMessenger::~DetectorMessenger() {
  delete fMatCmd;
  delete fDirPhantom;
  delete fDirStageA;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value) {
  if (command == fMatCmd) fDet->SetMaterialName(value);
}
