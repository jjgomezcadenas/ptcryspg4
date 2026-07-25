#ifndef ENSEMBLE_ACTIONS_HH
#define ENSEMBLE_ACTIONS_HH

#include "EnsembleConfig.hh"
#include "SamplingCurves.hh"

#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <random>
#include <unordered_map>
#include <vector>

class EnsembleDetector;
class EnsembleRun;
class G4Material;
class G4ParticleGun;

namespace ensemble {
/// One row of sampled_productions.csv, read back by the emitter-transport
/// mode and by the emitters.csv writer.
struct EmitterSeed {
  int event_id = 0;
  int isotope_index = 0;  ///< index into kBetaPlusResiduals
  double x_mm = 0., y_mm = 0., z_mm = 0.;
};
std::vector<EmitterSeed> ReadEmitterSeeds(const std::string& path);
}  // namespace ensemble

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

/// Emitter-transport mode: one event per sampled production, the residual
/// ion created at rest at its production point; Geant4's radioactive decay
/// and positron transport then set the annihilation point.
class EmitterPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
 public:
  explicit EmitterPrimaryGenerator(const ensemble::EnsembleCli& cli);
  ~EmitterPrimaryGenerator() override;
  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ParticleGun* fGun = nullptr;
  std::vector<ensemble::EmitterSeed> fSeeds;
};

/// Emitter-transport mode: associate the beta-plus positron of the primary
/// ion with its event and record its end point as the annihilation.
class EnsembleTrackingAction : public G4UserTrackingAction {
 public:
  void PreUserTrackingAction(const G4Track* track) override;
  void PostUserTrackingAction(const G4Track* track) override;

 private:
  G4int fPositronTrackID = -1;
};

/// Master-side output: run directory, exposure table, native counters,
/// depth dose, raw metadata; sampled productions when the sampler is on;
/// emitters.csv in emitter-transport mode.
class EnsembleRunAction : public G4UserRunAction {
 public:
  EnsembleRunAction(const ensemble::EnsembleCli& cli,
                    const EnsembleDetector& det,
                    const ensemble::SamplingCurves* curves);
  G4Run* GenerateRun() override;
  void EndOfRunAction(const G4Run* run) override;

 private:
  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
  const ensemble::SamplingCurves* fCurves;
};

/// Per-step tallies: proton exposure, native beta-plus residual production,
/// dose and depth dose; in-flight production sampling when curves are
/// loaded (from an engine separate from the transport stream).
class EnsembleSteppingAction : public G4UserSteppingAction {
 public:
  EnsembleSteppingAction(const ensemble::EnsembleCli& cli,
                         const EnsembleDetector& det,
                         const ensemble::SamplingCurves* curves);
  void UserSteppingAction(const G4Step* step) override;

 private:
  const std::array<G4double, ensemble::kNTargets>& Densities(
      const G4Material* material);

  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
  const ensemble::SamplingCurves* fCurves;
  std::mt19937_64 fSampleEngine;
  std::uniform_real_distribution<double> fUniform{0.0, 1.0};
  std::unordered_map<const G4Material*,
                     std::array<G4double, ensemble::kNTargets>>
      fDensityCache;  ///< target nuclei per cm3, per material
};

class EnsembleActionInitialization : public G4VUserActionInitialization {
 public:
  EnsembleActionInitialization(const ensemble::EnsembleCli& cli,
                               const EnsembleDetector& det,
                               const ensemble::SamplingCurves* curves)
      : fCli(cli), fDet(det), fCurves(curves) {}
  void BuildForMaster() const override;
  void Build() const override;

 private:
  const ensemble::EnsembleCli& fCli;
  const EnsembleDetector& fDet;
  const ensemble::SamplingCurves* fCurves;
};

#endif  // ENSEMBLE_ACTIONS_HH
