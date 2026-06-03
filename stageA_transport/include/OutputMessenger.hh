#ifndef STAGEA_OUTPUTMESSENGER_HH
#define STAGEA_OUTPUTMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

// Messenger exposing
//     /stageA/output/dir <path>
// so a macro can send the CSV outputs (emitters / run_meta / depth_dose) to the
// repo data/ dir. The shared "/stageA/" root is created by DetectorMessenger.
class OutputMessenger : public G4UImessenger {
 public:
  explicit OutputMessenger(RunAction* runAction);
  ~OutputMessenger() override;

  // The UI manager calls this on command entry; we route the path to RunAction.
  void SetNewValue(G4UIcommand* command, G4String newValue) override;

 private:
  RunAction* fRunAction = nullptr;
  G4UIdirectory* fDirOutput = nullptr;
  G4UIcmdWithAString* fDirCmd = nullptr;
};

#endif  // STAGEA_OUTPUTMESSENGER_HH
