#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

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

  fGeomCmd = new G4UIcmdWithAString("/stageA/phantom/geometry", this);
  fGeomCmd->SetGuidance("Phantom geometry: 'cylinder' (default, homogeneous) or "
                        "'mird_head' (heterogeneous scalp/skull/brain).");
  fGeomCmd->SetParameterName("name", false);
  fGeomCmd->SetCandidates("cylinder mird_head");
  fGeomCmd->AvailableForStates(G4State_PreInit);

  // --- target box (dose-normalization scoring region) -----------------------
  fDirTarget = new G4UIdirectory("/stageA/target/");
  fDirTarget->SetGuidance("Target (tumour) box for dose normalization.");

  fTRadiusCmd = new G4UIcmdWithADoubleAndUnit("/stageA/target/radius", this);
  fTRadiusCmd->SetGuidance("Target radius.");
  fTRadiusCmd->SetParameterName("r", false);
  fTRadiusCmd->SetUnitCategory("Length");
  fTRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTProxCmd = new G4UIcmdWithADoubleAndUnit("/stageA/target/proximal", this);
  fTProxCmd->SetGuidance("Target proximal depth (from the entrance face).");
  fTProxCmd->SetParameterName("d", false);
  fTProxCmd->SetUnitCategory("Length");
  fTProxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTDistCmd = new G4UIcmdWithADoubleAndUnit("/stageA/target/distal", this);
  fTDistCmd->SetGuidance("Target distal depth (from the entrance face).");
  fTDistCmd->SetParameterName("d", false);
  fTDistCmd->SetUnitCategory("Length");
  fTDistCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {
  delete fTDistCmd;
  delete fTProxCmd;
  delete fTRadiusCmd;
  delete fDirTarget;
  delete fGeomCmd;
  delete fMatCmd;
  delete fDirPhantom;
  delete fDirStageA;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value) {
  if (command == fMatCmd) {
    fDet->SetMaterialName(value);
  } else if (command == fGeomCmd) {
    fDet->SetGeometry(value);
  } else if (command == fTRadiusCmd) {
    fDet->SetTargetRadius(fTRadiusCmd->GetNewDoubleValue(value));
  } else if (command == fTProxCmd) {
    fDet->SetTargetProxDepth(fTProxCmd->GetNewDoubleValue(value));
  } else if (command == fTDistCmd) {
    fDet->SetTargetDistDepth(fTDistCmd->GetNewDoubleValue(value));
  }
}
