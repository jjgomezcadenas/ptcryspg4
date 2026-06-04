#ifndef STAGEA_BEAMMESSENGER_HH
#define STAGEA_BEAMMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class BeamConfig;
class G4UIdirectory;
class G4UIcmdWithAString;

// Exposes
//     /stageA/beam/layers <file>
// to load an SOBP energy-layer table (the field_design/sobp.py output) into the
// shared BeamConfig. Without it the gun fires a single fixed energy. Owned by
// the single BeamConfig, so the command is registered once on the master.
class BeamMessenger : public G4UImessenger {
 public:
  explicit BeamMessenger(BeamConfig* beam);
  ~BeamMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String value) override;

 private:
  BeamConfig* fBeam = nullptr;
  G4UIdirectory* fDirBeam = nullptr;
  G4UIcmdWithAString* fLayersCmd = nullptr;
};

#endif  // STAGEA_BEAMMESSENGER_HH
