#pragma once

#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "ThinTargetConfig.hh"

class G4ParticleGun;
class G4Run;
class G4Step;
class ThinTargetDetector;

class ThinTargetPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
 public:
  ThinTargetPrimaryGenerator(const ThinTargetConfig& config,
                             const ThinTargetDetector& detector);
  ~ThinTargetPrimaryGenerator() override;
  void GeneratePrimaries(G4Event* event) override;
 private:
  G4ParticleGun* gun_;
  const ThinTargetDetector& detector_;
};

class ThinTargetRunAction : public G4UserRunAction {
 public:
  ThinTargetRunAction(const ThinTargetConfig& config,
                      const ThinTargetDetector& detector);
  G4Run* GenerateRun() override;
  void EndOfRunAction(const G4Run* run) override;
 private:
  ThinTargetConfig config_;
  const ThinTargetDetector& detector_;
};

class ThinTargetSteppingAction : public G4UserSteppingAction {
 public:
  explicit ThinTargetSteppingAction(const ThinTargetDetector& detector)
      : detector_(detector) {}
  void UserSteppingAction(const G4Step* step) override;
 private:
  const ThinTargetDetector& detector_;
};

class ThinTargetActionInitialization : public G4VUserActionInitialization {
 public:
  ThinTargetActionInitialization(const ThinTargetConfig& config,
                                 const ThinTargetDetector& detector)
      : config_(config), detector_(detector) {}
  void BuildForMaster() const override;
  void Build() const override;
 private:
  ThinTargetConfig config_;
  const ThinTargetDetector& detector_;
};
