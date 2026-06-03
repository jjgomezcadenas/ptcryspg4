// Stage A: proton transport in a PMMA phantom, producing the detector-
// independent beta+ emitter source (spec sec 2). Milestone 1: geometry, gun,
// physics, MT transport, and interactive Qt visualization. PROD/ANH capture
// and the prod_anh.h5 writer come in later milestones.
//
// Usage:
//   ./proton_transport              -> interactive Qt session (executes vis.mac)
//   ./proton_transport run.mac      -> batch, runs the given macro

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "QGSP_BIC_HP.hh"  // hadronic list suited to proton isotope production

int main(int argc, char** argv) {
  // Detect interactive mode (no macro argument) and open the UI session first,
  // so G4UIExecutive can pick the Qt session this build provides.
  G4UIExecutive* uiExec = nullptr;
  if (argc == 1) {
    uiExec = new G4UIExecutive(argc, argv);
  }

  // Tasking run manager (G4's default MT engine; scales well on Apple Silicon).
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  auto* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(new QGSP_BIC_HP);
  runManager->SetUserInitialization(new ActionInitialization(det));

  // Visualization manager. This build has Qt + ToolsSG drivers; in batch it is
  // simply created and left idle.
  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  auto* ui = G4UImanager::GetUIpointer();

  if (uiExec) {
    // Interactive: set up the viewer and shoot a few protons to inspect.
    // Macros are copied next to the binary by CMake, so reference them by name.
    ui->ApplyCommand("/control/execute vis.mac");
    uiExec->SessionStart();
    delete uiExec;
  } else {
    // Batch: execute the supplied macro.
    ui->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }

  delete visManager;
  delete runManager;
  return 0;
}
