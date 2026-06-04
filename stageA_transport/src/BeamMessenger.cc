#include "BeamMessenger.hh"

#include "BeamConfig.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

BeamMessenger::BeamMessenger(BeamConfig* beam) : fBeam(beam) {
  // "/stageA/" already exists (DetectorMessenger, built first); add "beam/".
  fDirBeam = new G4UIdirectory("/stageA/beam/");
  fDirBeam->SetGuidance("Beam / gun controls.");

  fLayersCmd = new G4UIcmdWithAString("/stageA/beam/layers", this);
  fLayersCmd->SetGuidance("Load an SOBP energy-layer table (energy_MeV,weight CSV).");
  fLayersCmd->SetParameterName("file", false);
  fLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDiskCmd = new G4UIcmdWithADoubleAndUnit("/stageA/beam/disk", this);
  fDiskCmd->SetGuidance("Uniform-disk lateral profile of this radius "
                        "(0 = Gaussian pencil).");
  fDiskCmd->SetParameterName("radius", false);
  fDiskCmd->SetUnitCategory("Length");
  fDiskCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

BeamMessenger::~BeamMessenger() {
  delete fDiskCmd;
  delete fLayersCmd;
  delete fDirBeam;
}

void BeamMessenger::SetNewValue(G4UIcommand* command, G4String value) {
  if (command == fLayersCmd) {
    fBeam->LoadLayers(value);
  } else if (command == fDiskCmd) {
    fBeam->SetDiskRadius(fDiskCmd->GetNewDoubleValue(value));
  }
}
