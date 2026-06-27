#include "OutputMessenger.hh"

#include "RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

OutputMessenger::OutputMessenger(RunAction* runAction) : fRunAction(runAction) {
  // "/stageA/" already exists (created by DetectorMessenger, which is built
  // first); here we only add the output sub-directory and its command.
  fDirOutput = new G4UIdirectory("/stageA/output/");
  fDirOutput->SetGuidance("Output file controls.");

  fDirCmd = new G4UIcmdWithAString("/stageA/output/dir", this);
  fDirCmd->SetGuidance("Base output dir; Stage A appends a per-run subdir "
                       "<geometry>_<beam>_<N> and writes the CSVs there.");
  fDirCmd->SetParameterName("path", false);
  fDirCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

OutputMessenger::~OutputMessenger() {
  delete fDirCmd;
  delete fDirOutput;
}

void OutputMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fDirCmd) fRunAction->SetOutputDir(newValue);
}
