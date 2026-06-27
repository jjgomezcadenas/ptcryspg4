#ifndef STAGEA_RUNACTION_HH
#define STAGEA_RUNACTION_HH

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <string>

class G4Run;
class DetectorConstruction;
class BeamConfig;
class OutputMessenger;
class StageARun;

// Creates the custom StageARun on every thread. On the master at end of run it
// derives a per-run tag (<geometry>_<beam>_<N>) from what actually ran, makes a
// self-contained run directory <base>/<tag>/, and writes the CSV outputs there.
class RunAction : public G4UserRunAction {
 public:
  RunAction(G4bool onMaster, const DetectorConstruction* det,
            const BeamConfig* beam);
  ~RunAction() override;

  G4Run* GenerateRun() override;
  void EndOfRunAction(const G4Run* run) override;

  // The base directory under which a per-run subdir is created (set by macro).
  void SetOutputDir(const G4String& dir) { fBaseDir = dir; }

 private:
  // <geometry>_<beam>_<N>, auto-derived from the run (cannot disagree with it).
  std::string RunTag(const StageARun* run) const;

  void WriteEmittersCsv(const StageARun* run) const;
  void WriteMetaCsv(const StageARun* run) const;
  void WriteRegionsCsv() const;
  void WriteDepthDoseCsv(const StageARun* run) const;
  void PrintSummary(const StageARun* run) const;

  const DetectorConstruction* fDet = nullptr;
  const BeamConfig* fBeam = nullptr;
  OutputMessenger* fMessenger = nullptr;  // master only
  std::string fBaseDir = ".";             // base; run subdir appended at write
  std::string fOutputDir = ".";           // resolved <base>/<tag> for this run
};

#endif  // STAGEA_RUNACTION_HH
