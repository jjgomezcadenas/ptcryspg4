#ifndef STAGEA_ACTIONINITIALIZATION_HH
#define STAGEA_ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization {
 public:
  explicit ActionInitialization(const DetectorConstruction* det) : fDet(det) {}
  ~ActionInitialization() override = default;

  void Build() const override;        // per-worker thread
  void BuildForMaster() const override;  // master thread (MT)

 private:
  const DetectorConstruction* fDet = nullptr;
};

#endif  // STAGEA_ACTIONINITIALIZATION_HH
