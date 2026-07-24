#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Threading.hh"
#include "Randomize.hh"
#include "ThinTargetActions.hh"
#include "ThinTargetConfig.hh"
#include "ThinTargetDetector.hh"

namespace {
void usage() {
  std::cout << "Usage: thin_target --target C12|N14|O16 --energy-mev E "
               "--areal-mg-cm2 T --protons N --seed S --threads N "
               "--physics-list NAME --output FILE\n";
}

ThinTargetConfig parse(int argc, char** argv) {
  ThinTargetConfig config;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--help") { usage(); std::exit(0); }
    if (i + 1 >= argc) throw std::runtime_error("Missing value for " + argument);
    const std::string value = argv[++i];
    if (argument == "--target") config.target = value;
    else if (argument == "--energy-mev") config.energy_MeV = std::stod(value);
    else if (argument == "--areal-mg-cm2") config.areal_mg_cm2 = std::stod(value);
    else if (argument == "--protons") config.protons = std::stol(value);
    else if (argument == "--seed") config.seed = std::stol(value);
    else if (argument == "--threads") config.threads = std::stoi(value);
    else if (argument == "--physics-list") config.physics_list = value;
    else if (argument == "--output") config.output = value;
    else throw std::runtime_error("Unknown option: " + argument);
  }
  if (config.energy_MeV <= 0 || config.areal_mg_cm2 <= 0 || config.protons <= 0 || config.threads <= 0)
    throw std::runtime_error("Energy, thickness, protons and threads must be positive");
  return config;
}
}

int main(int argc, char** argv) {
  try {
    const auto config = parse(argc, argv);
    G4Random::setTheSeed(config.seed);
    auto* run_manager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
    run_manager->SetNumberOfThreads(config.threads);
    auto* detector = new ThinTargetDetector(config);
    run_manager->SetUserInitialization(detector);
    G4PhysListFactory factory;
    auto* physics = factory.GetReferencePhysList(config.physics_list);
    if (!physics) throw std::runtime_error("Unknown Geant4 physics list: " + std::string(config.physics_list));
    run_manager->SetUserInitialization(physics);
    run_manager->SetUserInitialization(new ThinTargetActionInitialization(config, *detector));
    run_manager->Initialize();
    run_manager->BeamOn(config.protons);
    delete run_manager;
  } catch (const std::exception& error) {
    std::cerr << "thin_target: " << error.what() << '\n';
    return 2;
  }
  return 0;
}
