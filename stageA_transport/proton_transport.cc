// Stage A: proton transport in a PMMA phantom, producing the detector-
// independent beta+ emitter source (spec sec 2). Milestone 1: geometry, gun,
// physics, MT transport, and interactive Qt visualization. PROD/ANH capture
// and the prod_anh.h5 writer come in later milestones.
//
// Usage:
//   ./proton_transport                 -> interactive Qt session (executes vis.mac)
//   ./proton_transport -i [macro]      -> interactive Qt session running [macro]
//                                         (default vis.mac) — for headep_vis.mac etc.
//   ./proton_transport run.mac         -> batch, runs the given macro

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "BeamConfig.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "QGSP_BIC_HP.hh"  // hadronic list suited to proton isotope production

int main(int argc, char** argv) {
  // Invocation modes:
  //   (no args)        interactive Qt, runs vis.mac
  //   -i [macro]       interactive Qt, runs [macro] (default vis.mac)
  //   macro            batch, runs macro
  const bool interactive = (argc == 1) || (argc >= 2 && G4String(argv[1]) == "-i");
  const G4String macro =
      interactive ? (argc >= 3 ? G4String(argv[2]) : G4String("vis.mac"))
                  : G4String(argv[1]);

  // Open the UI session first (with only argv[0], so the extra flag/macro do not
  // confuse the session-type autodetection), so G4UIExecutive can pick the Qt
  // session this build provides before any graphics system is created.
  G4UIExecutive* uiExec = interactive ? new G4UIExecutive(1, argv) : nullptr;

  // Tasking run manager (G4's default MT engine; scales well on Apple Silicon).
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // The three mandatory user initializations the run manager takes ownership of:
  // geometry (what), physics (how particles interact), actions (what we record).
  auto* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(new QGSP_BIC_HP);

  // Shared beam config: owns the /stageA/beam/ messenger and the SOBP layer
  // table (loaded by macro). Created after det so "/stageA/" already exists.
  auto* beam = new BeamConfig;
  runManager->SetUserInitialization(new ActionInitialization(det, beam));

  // Visualization manager. This build has Qt + ToolsSG drivers; in batch it is
  // simply created and left idle.
  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  auto* ui = G4UImanager::GetUIpointer();

  ui->ApplyCommand(G4String("/control/execute ") + macro);
  if (uiExec) {
    // Interactive: keep the Qt window open for inspection after the macro runs.
    uiExec->SessionStart();
    delete uiExec;
  }

  delete visManager;
  delete runManager;
  delete beam;
  return 0;
}
