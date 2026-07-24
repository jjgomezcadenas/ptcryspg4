#include "ThinTargetActions.hh"

#include <filesystem>
#include <fstream>
#include <iomanip>

#include "G4Event.hh"
#include "G4GenericIon.hh"
#include "G4HadronicProcessType.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4Version.hh"
#include "ThinTargetDetector.hh"
#include "ThinTargetRun.hh"

ThinTargetPrimaryGenerator::ThinTargetPrimaryGenerator(
    const ThinTargetConfig& config, const ThinTargetDetector& detector)
    : gun_(new G4ParticleGun(1)), detector_(detector) {
  gun_->SetParticleDefinition(G4Proton::Definition());
  gun_->SetParticleEnergy(config.energy_MeV * MeV);
  gun_->SetParticleMomentumDirection({0.0, 0.0, 1.0});
}

ThinTargetPrimaryGenerator::~ThinTargetPrimaryGenerator() { delete gun_; }

void ThinTargetPrimaryGenerator::GeneratePrimaries(G4Event* event) {
  gun_->SetParticlePosition({0.0, 0.0, -(0.5 * detector_.ThicknessMm() + 0.001) * mm});
  gun_->GeneratePrimaryVertex(event);
  auto* run = static_cast<ThinTargetRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountIncident();
}

ThinTargetRunAction::ThinTargetRunAction(const ThinTargetConfig& config,
                                         const ThinTargetDetector& detector)
    : config_(config), detector_(detector) {}

G4Run* ThinTargetRunAction::GenerateRun() { return new ThinTargetRun; }

namespace {
G4double mean_or_zero(G4double sum, G4long count) { return count ? sum / count : 0.0; }
}

void ThinTargetRunAction::EndOfRunAction(const G4Run* base_run) {
  if (!IsMaster()) return;
  const auto* run = static_cast<const ThinTargetRun*>(base_run);
  const G4double atomic_mass = detector_.TargetA();
  const G4double nuclear_areal_density =
      config_.areal_mg_cm2 * 1.0e-3 / atomic_mass * 6.02214076e23;
  std::filesystem::path output(std::string(config_.output));
  if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
  const bool write_header = !std::filesystem::exists(output) || std::filesystem::file_size(output) == 0;
  std::ofstream stream(output, std::ios::app);
  if (write_header) {
    stream << "target,target_z,target_a,energy_MeV,areal_mg_cm2,density_g_cm3,"
              "thickness_mm,nuclear_areal_density_cm2,n_protons,n_inelastic,n_c11,n_n13,n_o15,"
              "mean_continuous_loss_MeV,mean_inelastic_energy_MeV,mean_c11_energy_MeV,"
              "mean_n13_energy_MeV,mean_o15_energy_MeV,physics_list,geant4_version,seed,threads\n";
  }
  G4String version = G4Version;
  for (auto& character : version) if (character == ',') character = ' ';
  stream << std::setprecision(12)
         << config_.target << ',' << detector_.TargetZ() << ',' << detector_.TargetA() << ','
         << config_.energy_MeV << ',' << config_.areal_mg_cm2 << ',' << detector_.DensityGcm3() << ','
         << detector_.ThicknessMm() << ',' << nuclear_areal_density << ',' << run->Incident() << ','
         << run->Inelastic() << ',' << run->C11() << ',' << run->N13() << ',' << run->O15() << ','
         << mean_or_zero(run->ContinuousLoss(), run->Incident()) / MeV << ','
         << mean_or_zero(run->InelasticEnergySum(), run->Inelastic()) / MeV << ','
         << mean_or_zero(run->C11EnergySum(), run->C11()) / MeV << ','
         << mean_or_zero(run->N13EnergySum(), run->N13()) / MeV << ','
         << mean_or_zero(run->O15EnergySum(), run->O15()) / MeV << ','
         << config_.physics_list << ',' << version << ',' << config_.seed << ',' << config_.threads << '\n';
}

void ThinTargetSteppingAction::UserSteppingAction(const G4Step* step) {
  const auto* track = step->GetTrack();
  if (track->GetParticleDefinition() != G4Proton::Definition()) return;
  if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()
      != detector_.TargetLogical()) return;
  const auto* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  const bool inelastic = process && process->GetProcessSubType() == fHadronInelastic;
  auto* run = static_cast<ThinTargetRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  if (track->GetParentID() == 0 && !inelastic) {
    const G4double loss = step->GetPreStepPoint()->GetKineticEnergy()
                          - step->GetPostStepPoint()->GetKineticEnergy();
    if (loss > 0.0) run->AddContinuousLoss(loss);
  }
  if (!inelastic) return;
  const G4double interaction_energy = step->GetPreStepPoint()->GetKineticEnergy();
  run->CountInelastic(interaction_energy);
  for (const auto* secondary : *step->GetSecondaryInCurrentStep()) {
    const auto* definition = secondary->GetParticleDefinition();
    if (definition->GetParticleType() != "nucleus") continue;
    const G4int z = definition->GetAtomicNumber();
    const G4int a = definition->GetAtomicMass();
    run->CountResidual(z, a, interaction_energy);
  }
}

void ThinTargetActionInitialization::BuildForMaster() const {
  SetUserAction(new ThinTargetRunAction(config_, detector_));
}

void ThinTargetActionInitialization::Build() const {
  SetUserAction(new ThinTargetPrimaryGenerator(config_, detector_));
  SetUserAction(new ThinTargetRunAction(config_, detector_));
  SetUserAction(new ThinTargetSteppingAction(detector_));
}
