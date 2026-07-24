#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include "G4CrossSectionDataStore.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicProcessType.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysListFactory.hh"
#include "G4ProcessManager.hh"
#include "G4Proton.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4VParticleChange.hh"
#include "G4Version.hh"
#include "Randomize.hh"
#include "ThinTargetConfig.hh"
#include "ThinTargetDetector.hh"

namespace {

struct SamplerConfig {
  G4String target = "O16";
  G4double energy_MeV = 100.0;
  G4long interactions = 10000;
  G4long seed = 12345;
  G4String physics_list = "QGSP_BIC_HP";
  G4String output = "denominator_counts.csv";
};

void usage() {
  std::cout << "Usage: denominator_sampler --target C12|N14|O16 "
               "--energy-mev E --interactions N --seed S "
               "--physics-list NAME --output FILE\n";
}

SamplerConfig parse(int argc, char** argv) {
  SamplerConfig config;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--help") {
      usage();
      std::exit(0);
    }
    if (i + 1 >= argc) {
      throw std::runtime_error("Missing value for " + argument);
    }
    const std::string value = argv[++i];
    if (argument == "--target") config.target = value;
    else if (argument == "--energy-mev") config.energy_MeV = std::stod(value);
    else if (argument == "--interactions") config.interactions = std::stol(value);
    else if (argument == "--seed") config.seed = std::stol(value);
    else if (argument == "--physics-list") config.physics_list = value;
    else if (argument == "--output") config.output = value;
    else throw std::runtime_error("Unknown option: " + argument);
  }
  if (config.energy_MeV <= 0.0 || config.interactions <= 0) {
    throw std::runtime_error("Energy and interaction count must be positive");
  }
  return config;
}

G4HadronicProcess* proton_inelastic_process() {
  auto* manager = G4Proton::Definition()->GetProcessManager();
  auto* processes = manager->GetProcessList();
  for (G4int index = 0; index < manager->GetProcessListLength(); ++index) {
    auto* process = (*processes)[index];
    if (process->GetProcessSubType() == fHadronInelastic) {
      auto* hadronic = dynamic_cast<G4HadronicProcess*>(process);
      if (hadronic) return hadronic;
    }
  }
  throw std::runtime_error("The selected physics list has no proton inelastic process");
}

std::string data_set_names(G4HadronicProcess* process) {
  std::ostringstream names;
  const auto& data_sets = process->GetCrossSectionDataStore()->GetDataSetList();
  for (std::size_t index = 0; index < data_sets.size(); ++index) {
    if (index) names << '|';
    names << data_sets[index]->GetName();
  }
  return names.str();
}

std::string serialized_counts(const std::map<std::string, G4long>& counts) {
  std::ostringstream output;
  bool first = true;
  for (const auto& [name, count] : counts) {
    if (!first) output << '|';
    first = false;
    output << name << ':' << count;
  }
  return output.str();
}

void write_result(const SamplerConfig& config, const ThinTargetDetector& detector,
                  G4double sigma_inelastic, G4long c11, G4long n13, G4long o15,
                  G4long secondary_count, G4long nucleus_count,
                  const std::map<std::string, G4long>& residual_counts,
                  const std::string& data_sets,
                  const std::map<std::string, G4long>& model_counts) {
  std::filesystem::path output(std::string(config.output));
  if (output.has_parent_path()) {
    std::filesystem::create_directories(output.parent_path());
  }
  const bool header = (!std::filesystem::exists(output)
                       || std::filesystem::file_size(output) == 0);
  std::ofstream stream(output, std::ios::app);
  if (header) {
    stream << "target,target_z,target_a,energy_MeV,n_interactions,n_c11,n_n13,n_o15,"
              "sigma_inelastic_mb,n_secondaries,n_nuclei,residual_counts,"
              "cross_section_data_sets,model_counts,physics_list,"
              "geant4_version,seed\n";
  }
  G4String version = G4Version;
  for (auto& character : version) if (character == ',') character = ' ';
  stream << std::setprecision(12)
         << config.target << ',' << detector.TargetZ() << ',' << detector.TargetA() << ','
         << config.energy_MeV << ',' << config.interactions << ','
         << c11 << ',' << n13 << ',' << o15 << ',' << sigma_inelastic / millibarn << ','
         << secondary_count << ',' << nucleus_count << ','
         << serialized_counts(residual_counts) << ',' << data_sets << ','
         << serialized_counts(model_counts) << ','
         << config.physics_list << ',' << version << ',' << config.seed << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto config = parse(argc, argv);
    G4Random::setTheSeed(config.seed);

    ThinTargetConfig detector_config;
    detector_config.target = config.target;
    detector_config.areal_mg_cm2 = 5.0;
    auto* run_manager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    auto* detector = new ThinTargetDetector(detector_config);
    run_manager->SetUserInitialization(detector);
    G4PhysListFactory factory;
    auto* physics = factory.GetReferencePhysList(config.physics_list);
    if (!physics) {
      throw std::runtime_error(
          "Unknown Geant4 physics list: " + std::string(config.physics_list));
    }
    run_manager->SetUserInitialization(physics);
    run_manager->Initialize();

    auto* process = proton_inelastic_process();
    auto* store = process->GetCrossSectionDataStore();
    G4DynamicParticle query(
        G4Proton::Definition(), G4ThreeVector(0.0, 0.0, 1.0),
        config.energy_MeV * MeV);
    const G4double sigma_inelastic = store->GetCrossSection(
        &query, detector->TargetZ(), detector->TargetA(),
        detector->TargetIsotope(), detector->TargetElement(),
        detector->TargetMaterial());

    G4long c11 = 0;
    G4long n13 = 0;
    G4long o15 = 0;
    G4long secondary_count = 0;
    G4long nucleus_count = 0;
    std::map<std::string, G4long> residual_counts;
    std::map<std::string, G4long> model_counts;
    auto* couple = detector->TargetLogical()->GetMaterialCutsCouple();
    G4MaterialCutsCouple fallback_couple(detector->TargetMaterial());
    if (!couple) {
      // No event loop is started, so the logical-volume couple may not yet be
      // assigned.  The hadronic process needs its material through the step.
      couple = &fallback_couple;
    }
    for (G4long interaction = 0; interaction < config.interactions; ++interaction) {
      auto* particle = new G4DynamicParticle(
          G4Proton::Definition(), G4ThreeVector(0.0, 0.0, 1.0),
          config.energy_MeV * MeV);
      G4Track track(particle, 0.0, G4ThreeVector());
      G4Step step;
      step.SetTrack(&track);
      step.GetPreStepPoint()->SetMaterial(detector->TargetMaterial());
      step.GetPostStepPoint()->SetMaterial(detector->TargetMaterial());
      step.GetPreStepPoint()->SetMaterialCutsCouple(couple);
      step.GetPostStepPoint()->SetMaterialCutsCouple(couple);
      track.SetStep(&step);
      store->ComputeCrossSection(track.GetDynamicParticle(), detector->TargetMaterial());
      auto* change = process->PostStepDoIt(track, step);
      const auto* interaction_model = process->GetHadronicInteraction();
      const std::string model = interaction_model
          ? std::string(interaction_model->GetModelName()) : "unknown";
      ++model_counts[model];
      const G4int secondaries = change->GetNumberOfSecondaries();
      secondary_count += secondaries;
      for (G4int index = 0; index < secondaries; ++index) {
        auto* secondary = change->GetSecondary(index);
        const auto* definition = secondary->GetParticleDefinition();
        if (definition->GetParticleType() == "nucleus") {
          ++nucleus_count;
          const G4int z = definition->GetAtomicNumber();
          const G4int a = definition->GetAtomicMass();
          ++residual_counts[std::to_string(z) + "-" + std::to_string(a)];
          if (z == 6 && a == 11) ++c11;
          if (z == 7 && a == 13) ++n13;
          if (z == 8 && a == 15) ++o15;
        }
        delete secondary;
      }
      change->Clear();
    }

    write_result(config, *detector, sigma_inelastic, c11, n13, o15,
                 secondary_count, nucleus_count, residual_counts,
                 data_set_names(process), model_counts);
    delete run_manager;
  } catch (const std::exception& error) {
    std::cerr << "denominator_sampler: " << error.what() << '\n';
    return 2;
  }
  return 0;
}
