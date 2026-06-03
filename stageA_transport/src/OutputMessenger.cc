#include "OutputMessenger.hh"

#include "RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

OutputMessenger::OutputMessenger(RunAction* runAction) : fRunAction(runAction) {
  fDirStageA = new G4UIdirectory("/stageA/");
  fDirStageA->SetGuidance("Stage-A proton-transport controls.");

  fDirOutput = new G4UIdirectory("/stageA/output/");
  fDirOutput->SetGuidance("Output file controls.");

  fDirCmd = new G4UIcmdWithAString("/stageA/output/dir", this);
  fDirCmd->SetGuidance("Directory for emitters.csv and run_meta.csv.");
  fDirCmd->SetParameterName("path", false);
  fDirCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

OutputMessenger::~OutputMessenger() {
  delete fDirCmd;
  delete fDirOutput;
  delete fDirStageA;
}

void OutputMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fDirCmd) {
    fRunAction->SetOutputDir(newValue);
  }
}
