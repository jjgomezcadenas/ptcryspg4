#ifndef STAGEA_DETECTORMESSENGER_HH
#define STAGEA_DETECTORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

// A "messenger" is Geant4's bridge between typed UI/macro commands and C++
// objects. This one exposes the phantom material
//     /stageA/phantom/material <NIST name>     (before /run/initialize)
// and the target box for dose normalization
//     /stageA/target/radius|proximal|distal <length>
// so geometry/target parameters can be set from a macro without recompiling.
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
  G4UIdirectory* fDirTarget = nullptr;    // "/stageA/target/"
  G4UIcmdWithAString* fMatCmd = nullptr;  // "/stageA/phantom/material"
  G4UIcmdWithAString* fGeomCmd = nullptr;  // "/stageA/phantom/geometry"
  G4UIcmdWithADoubleAndUnit* fTRadiusCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fTProxCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fTDistCmd = nullptr;
};

#endif  // STAGEA_DETECTORMESSENGER_HH
