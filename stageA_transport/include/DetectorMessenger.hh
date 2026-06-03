#ifndef STAGEA_DETECTORMESSENGER_HH
#define STAGEA_DETECTORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

// A "messenger" is Geant4's bridge between typed UI/macro commands and C++
// objects. This one exposes
//     /stageA/phantom/material <NIST name>
// so the phantom material (G4_TISSUE_SOFT_ICRP / G4_PLEXIGLASS / G4_WATER …)
// can be chosen from a macro without recompiling. Must be issued before
// /run/initialize, since the geometry is built then.
class DetectorMessenger : public G4UImessenger {
 public:
  explicit DetectorMessenger(DetectorConstruction* det);
  ~DetectorMessenger() override;

  // The UI manager calls this when a registered command is entered; we route
  // the typed value to the DetectorConstruction.
  void SetNewValue(G4UIcommand* command, G4String value) override;

 private:
  DetectorConstruction* fDet = nullptr;
  G4UIdirectory* fDirStageA = nullptr;    // "/stageA/" (shared root)
  G4UIdirectory* fDirPhantom = nullptr;   // "/stageA/phantom/"
  G4UIcmdWithAString* fMatCmd = nullptr;  // "/stageA/phantom/material"
};

#endif  // STAGEA_DETECTORMESSENGER_HH
