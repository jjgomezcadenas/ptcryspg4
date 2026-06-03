#ifndef STAGEA_RUNACTION_HH
#define STAGEA_RUNACTION_HH

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <string>

class G4Run;
class DetectorConstruction;
class OutputMessenger;
class StageARun;

// Creates the custom StageARun on every thread. On the master at end of run it
// writes data/emitters.csv + data/run_meta.csv and prints a short summary.
class RunAction : public G4UserRunAction {
 public:
  RunAction(G4bool onMaster, const DetectorConstruction* det);
  ~RunAction() override;

  G4Run* GenerateRun() override;
  void EndOfRunAction(const G4Run* run) override;

  void SetOutputDir(const G4String& dir) { fOutputDir = dir; }

 private:
  void WriteEmittersCsv(const StageARun* run) const;
  void WriteMetaCsv(const StageARun* run) const;
  void WriteDepthDoseCsv(const StageARun* run) const;
  void PrintSummary(const StageARun* run) const;

  const DetectorConstruction* fDet = nullptr;
  OutputMessenger* fMessenger = nullptr;  // master only
  std::string fOutputDir = ".";
};

#endif  // STAGEA_RUNACTION_HH
