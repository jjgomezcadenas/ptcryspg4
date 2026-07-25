#ifndef ENSEMBLE_ACTIONS_HH
#define ENSEMBLE_ACTIONS_HH

#include "EnsembleConfig.hh"

#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <unordered_map>
#include <vector>

class EnsembleDetector;
class EnsembleRun;
class G4Material;
class G4ParticleGun;

/// Proton gun along +z: SOBP layer sampling with a uniform fluence disk when a
/// layer table is loaded, otherwise a Gaussian pencil at the fixed energy.
class EnsemblePrimaryGenerator : public G4VUserPrimaryGeneratorAction {
 public:
  EnsemblePrimaryGenerator(const ensemble::EnsembleCli& cli,
                           const EnsembleDetector& det);
  ~EnsemblePrimaryGenerator() override;
  void GeneratePrimaries(G4Event* event) override;

 private:
  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
  G4ParticleGun* fGun = nullptr;
  std::vector<G4double> fLayerEnergyMeV;
  std::vector<G4double> fLayerCumWeight;
};

/// Master-side output: run directory, exposure table, native counters,
/// depth dose, raw metadata.
class EnsembleRunAction : public G4UserRunAction {
 public:
  EnsembleRunAction(const ensemble::EnsembleCli& cli,
                    const EnsembleDetector& det);
  G4Run* GenerateRun() override;
  void EndOfRunAction(const G4Run* run) override;

 private:
  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
};

/// Per-step tallies: proton exposure, native beta-plus residual production,
/// dose and depth dose.
class EnsembleSteppingAction : public G4UserSteppingAction {
 public:
  EnsembleSteppingAction(const ensemble::EnsembleCli& cli,
                         const EnsembleDetector& det);
  void UserSteppingAction(const G4Step* step) override;

 private:
  const std::array<G4double, ensemble::kNTargets>& Densities(
      const G4Material* material);

  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
  std::unordered_map<const G4Material*,
                     std::array<G4double, ensemble::kNTargets>>
      fDensityCache;  ///< target nuclei per cm3, per material
};

class EnsembleActionInitialization : public G4VUserActionInitialization {
 public:
  EnsembleActionInitialization(const ensemble::EnsembleCli& cli,
                               const EnsembleDetector& det)
      : fCli(cli), fDet(det) {}
  void BuildForMaster() const override;
  void Build() const override;

 private:
  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
};

#endif  // ENSEMBLE_ACTIONS_HH
