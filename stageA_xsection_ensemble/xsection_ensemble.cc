/// Exposure application of the cross-section ensemble calculation
/// (docs/xsections_plan.tex; workshop/xsections_phases.md phase 1a).
///
/// Transport mode: protons through the phantom, writing the target-resolved
/// proton exposure table, native-route counters, depth dose and raw metadata.
///
/// Step-table mode (--steps-csv): accumulate a synthetic proton-step CSV into
/// an exposure table with the same code path as the transport scorer, for the
/// C++/Python equality test against accumulate_step_exposure.

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "EnsembleActions.hh"
#include "EnsembleConfig.hh"
#include "EnsembleDetector.hh"
#include "ExposureGrid.hh"
#include "G4PhysListFactory.hh"
#include "G4RunManagerFactory.hh"
#include "Randomize.hh"

namespace {

void usage() {
  std::cout
      << "Transport: xsection_ensemble --protons N [--geometry "
         "uniform_headep|uniform_head|cylinder]\n"
         "  [--layers sobp.csv] [--disk-radius-mm R] [--energy-mev E]\n"
         "  [--target-radius-mm R --target-prox-mm P --target-dist-mm D]\n"
         "  [--material NAME --phantom-diameter-mm D --phantom-length-mm L]\n"
         "  [--seed S] [--threads N] [--physics-list NAME]\n"
         "  [--ebin-width-mev W] [--emax-mev E] [--zbin-width-mm W]\n"
         "  [--output-base DIR]\n"
         "Step table: xsection_ensemble --steps-csv IN --steps-out OUT "
         "--zmax-mm Z\n"
         "  [--ebin-width-mev W] [--emax-mev E] [--zbin-width-mm W]\n";
}

ensemble::EnsembleCli parse(int argc, char** argv) {
  ensemble::EnsembleCli cli;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help") { usage(); std::exit(0); }
    if (i + 1 >= argc) throw std::runtime_error("missing value for " + arg);
    const std::string value = argv[++i];
    if (arg == "--geometry") cli.geometry = value;
    else if (arg == "--layers") cli.layers = value;
    else if (arg == "--disk-radius-mm") cli.disk_radius_mm = std::stod(value);
    else if (arg == "--energy-mev") cli.energy_MeV = std::stod(value);
    else if (arg == "--beam-sigma-mm") cli.beam_sigma_mm = std::stod(value);
    else if (arg == "--target-radius-mm") cli.target_radius_mm = std::stod(value);
    else if (arg == "--target-prox-mm") cli.target_prox_mm = std::stod(value);
    else if (arg == "--target-dist-mm") cli.target_dist_mm = std::stod(value);
    else if (arg == "--material") cli.material = value;
    else if (arg == "--phantom-diameter-mm") cli.phantom_diameter_mm = std::stod(value);
    else if (arg == "--phantom-length-mm") cli.phantom_length_mm = std::stod(value);
    else if (arg == "--protons") cli.protons = std::stol(value);
    else if (arg == "--seed") cli.seed = std::stol(value);
    else if (arg == "--threads") cli.threads = std::stoi(value);
    else if (arg == "--physics-list") cli.physics_list = value;
    else if (arg == "--ebin-width-mev") cli.ebin_width_MeV = std::stod(value);
    else if (arg == "--emax-mev") cli.emax_MeV = std::stod(value);
    else if (arg == "--zbin-width-mm") cli.zbin_width_mm = std::stod(value);
    else if (arg == "--output-base") cli.output_base = value;
    else if (arg == "--steps-csv") cli.steps_csv = value;
    else if (arg == "--steps-out") cli.steps_out = value;
    else if (arg == "--zmax-mm") cli.zmax_mm = std::stod(value);
    else throw std::runtime_error("unknown option: " + arg);
  }
  return cli;
}

/// Offline accumulation of a synthetic step table
/// (columns of accumulate_step_exposure: target, proton_weight,
/// target_number_density_cm3, step_length_cm, energy_MeV, depth_mm).
int RunStepTableMode(const ensemble::EnsembleCli& cli) {
  if (cli.steps_out.empty() || cli.zmax_mm <= 0.)
    throw std::runtime_error("--steps-csv requires --steps-out and --zmax-mm");
  ensemble::ExposureGrid grid;
  grid.Init(cli.ebin_width_MeV, cli.emax_MeV, cli.zbin_width_mm, cli.zmax_mm);

  std::ifstream f(cli.steps_csv);
  if (!f) throw std::runtime_error("cannot open " + cli.steps_csv);
  std::string line;
  std::getline(f, line);  // header
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string target, weight, density, length, energy, depth;
    std::getline(ss, target, ',');
    std::getline(ss, weight, ',');
    std::getline(ss, density, ',');
    std::getline(ss, length, ',');
    std::getline(ss, energy, ',');
    std::getline(ss, depth, ',');
    int t = -1;
    for (int i = 0; i < ensemble::kNTargets; ++i)
      if (target == ensemble::kTargetNames[i]) t = i;
    if (t < 0) throw std::runtime_error("unknown target: " + target);
    grid.Add(t, std::stod(energy), std::stod(depth),
             std::stod(weight) * std::stod(density) * std::stod(length));
  }
  if (grid.Overflow() > 0)
    throw std::runtime_error("step-table overflow: " +
                             std::to_string(grid.Overflow()));
  grid.WriteCsv(cli.steps_out);
  std::cout << "accumulated step table -> " << cli.steps_out << "\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto cli = parse(argc, argv);
    if (!cli.steps_csv.empty()) return RunStepTableMode(cli);
    if (cli.protons <= 0)
      throw std::runtime_error("--protons must be positive (or use --steps-csv)");

    G4Random::setTheSeed(cli.seed);
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
    runManager->SetNumberOfThreads(cli.threads);
    auto* detector = new EnsembleDetector(cli);
    runManager->SetUserInitialization(detector);
    G4PhysListFactory factory;
    auto* physics = factory.GetReferencePhysList(cli.physics_list);
    if (!physics)
      throw std::runtime_error("unknown Geant4 physics list: " +
                               cli.physics_list);
    runManager->SetUserInitialization(physics);
    runManager->SetUserInitialization(
        new EnsembleActionInitialization(cli, *detector));
    runManager->Initialize();
    runManager->BeamOn(static_cast<G4int>(cli.protons));
    delete runManager;
  } catch (const std::exception& error) {
    std::cerr << "xsection_ensemble: " << error.what() << "\n";
    return 2;
  }
  return 0;
}
