#ifndef STAGEA_OUTPUTMESSENGER_HH
#define STAGEA_OUTPUTMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

// Exposes /stageA/output/dir <path> so the macro can send the CSVs to data/.
class OutputMessenger : public G4UImessenger {
 public:
  explicit OutputMessenger(RunAction* runAction);
  ~OutputMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String newValue) override;

 private:
  RunAction* fRunAction = nullptr;
  G4UIdirectory* fDirStageA = nullptr;
  G4UIdirectory* fDirOutput = nullptr;
  G4UIcmdWithAString* fDirCmd = nullptr;
};

#endif  // STAGEA_OUTPUTMESSENGER_HH
