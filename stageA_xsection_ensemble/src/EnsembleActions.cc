#include "EnsembleActions.hh"

#include "EnsembleDetector.hh"
#include "EnsembleRun.hh"

#include "G4Event.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicProcessType.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4Version.hh"
#include "Randomize.hh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {

// Compact proton-count rendering for the run tag (stageA convention).
std::string FormatProtonCount(long n) {
  if (n <= 0) return std::to_string(n);
  long m = n;
  int k = 0;
  while (m % 10 == 0) { m /= 10; ++k; }
  if (m < 10 && k > 0) return std::to_string(m) + "e" + std::to_string(k);
  return std::to_string(n);
}

// Element symbols for target labels of the native-route counters.
std::string NucleusLabel(int Z, int A) {
  static const char* symbols[] = {"?",  "H",  "He", "Li", "Be", "B",  "C",
                                  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al",
                                  "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca"};
  const std::string symbol =
      (Z >= 1 && Z <= 20) ? symbols[Z] : ("Z" + std::to_string(Z));
  return symbol + std::to_string(A);
}

}  // namespace

namespace ensemble {

std::vector<EmitterSeed> ReadEmitterSeeds(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  std::vector<EmitterSeed> seeds;
  std::string line;
  std::getline(f, line);  // header
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string field;
    EmitterSeed seed;
    std::getline(ss, field, ',');  // event_id
    seed.event_id = std::stoi(field);
    std::getline(ss, field, ',');  // channel_id
    std::getline(ss, field, ',');  // target
    std::getline(ss, field, ',');  // residual
    std::getline(ss, field, ',');  // isotope_id
    seed.isotope_index = std::stoi(field);
    std::getline(ss, field, ',');
    seed.x_mm = std::stod(field);
    std::getline(ss, field, ',');
    seed.y_mm = std::stod(field);
    std::getline(ss, field, ',');
    seed.z_mm = std::stod(field);
    seeds.push_back(seed);
  }
  return seeds;
}

}  // namespace ensemble

// ---------------------------------------------------------------------------
// Emitter-transport mode

EmitterPrimaryGenerator::EmitterPrimaryGenerator(
    const ensemble::EnsembleCli& cli)
    : fGun(new G4ParticleGun(1)),
      fSeeds(ensemble::ReadEmitterSeeds(cli.emitters_in)) {}

EmitterPrimaryGenerator::~EmitterPrimaryGenerator() { delete fGun; }

void EmitterPrimaryGenerator::GeneratePrimaries(G4Event* event) {
  const auto& seed = fSeeds.at(event->GetEventID());
  const auto& residual = ensemble::kBetaPlusResiduals[seed.isotope_index];
  auto* ion = G4IonTable::GetIonTable()->GetIon(residual.Z, residual.A, 0.);
  fGun->SetParticleDefinition(ion);
  fGun->SetParticleCharge(0.);
  fGun->SetParticleEnergy(0.);
  fGun->SetParticlePosition(
      {seed.x_mm * mm, seed.y_mm * mm, seed.z_mm * mm});
  fGun->GeneratePrimaryVertex(event);
}

void EnsembleTrackingAction::PreUserTrackingAction(const G4Track* track) {
  // The beta-plus positron is the positron child of the primary ion.
  if (track->GetParticleDefinition() == G4Positron::Definition() &&
      track->GetParentID() == 1) {
    fPositronTrackID = track->GetTrackID();
  }
}

void EnsembleTrackingAction::PostUserTrackingAction(const G4Track* track) {
  if (track->GetTrackID() != fPositronTrackID) return;
  fPositronTrackID = -1;
  auto* run = static_cast<EnsembleRun*>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  const auto* event =
      G4RunManager::GetRunManager()->GetCurrentEvent();
  const auto& position = track->GetPosition();
  run->RecordAnnihilation(event->GetEventID(), position.x() / mm,
                          position.y() / mm, position.z() / mm);
}

// ---------------------------------------------------------------------------
// Primary generator

EnsemblePrimaryGenerator::EnsemblePrimaryGenerator(
    const ensemble::EnsembleCli& cli, const EnsembleDetector& det)
    : fCli(cli), fDet(det), fGun(new G4ParticleGun(1)) {
  fGun->SetParticleDefinition(G4Proton::Definition());
  fGun->SetParticleMomentumDirection({0., 0., 1.});

  if (!fCli.layers.empty()) {
    std::ifstream f(fCli.layers);
    if (!f)
      G4Exception("EnsemblePrimaryGenerator", "BadLayerTable", FatalException,
                  ("cannot open SOBP layer file " + fCli.layers).c_str());
    std::string line;
    std::getline(f, line);  // header energy_MeV,weight
    std::vector<G4double> weights;
    while (std::getline(f, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line);
      std::string e, w;
      std::getline(ss, e, ',');
      std::getline(ss, w, ',');
      const G4double energy = std::stod(e);
      const G4double weight = std::stod(w);
      if (energy <= 0. || weight < 0.)
        G4Exception("EnsemblePrimaryGenerator", "BadLayerTable", FatalException,
                    ("invalid layer line: " + line).c_str());
      fLayerEnergyMeV.push_back(energy);
      weights.push_back(weight);
    }
    G4double sum = 0.;
    for (G4double w : weights) sum += w;
    if (fLayerEnergyMeV.empty() || sum <= 0.)
      G4Exception("EnsemblePrimaryGenerator", "BadLayerTable", FatalException,
                  "empty layer table or zero total weight");
    G4double cum = 0.;
    for (G4double w : weights) {
      cum += w / sum;
      fLayerCumWeight.push_back(cum);
    }
  }
}

EnsemblePrimaryGenerator::~EnsemblePrimaryGenerator() { delete fGun; }

void EnsemblePrimaryGenerator::GeneratePrimaries(G4Event* event) {
  G4double energy = fCli.energy_MeV * MeV;
  if (!fLayerEnergyMeV.empty()) {
    const G4double u = G4UniformRand();
    std::size_t i = 0;
    while (i + 1 < fLayerCumWeight.size() && fLayerCumWeight[i] < u) ++i;
    energy = fLayerEnergyMeV[i] * MeV;
  }
  fGun->SetParticleEnergy(energy);

  G4double x, y;
  if (!fLayerEnergyMeV.empty() && fCli.disk_radius_mm > 0.) {
    const G4double r = fCli.disk_radius_mm * mm * std::sqrt(G4UniformRand());
    const G4double phi = twopi * G4UniformRand();
    x = r * std::cos(phi);
    y = r * std::sin(phi);
  } else {
    x = G4RandGauss::shoot(0., fCli.beam_sigma_mm * mm);
    y = G4RandGauss::shoot(0., fCli.beam_sigma_mm * mm);
  }
  fGun->SetParticlePosition({x, y, -fDet.BeamAxisHalfExtent() - 1. * mm});
  fGun->GeneratePrimaryVertex(event);
}

// ---------------------------------------------------------------------------
// Stepping action

EnsembleSteppingAction::EnsembleSteppingAction(
    const ensemble::EnsembleCli& cli, const EnsembleDetector& det,
    const ensemble::SamplingCurves* curves)
    : fCli(cli), fDet(det), fCurves(curves),
      // Production sampling draws from its own engine, so the Geant4
      // transport stream is identical with the sampler on or off.
      fSampleEngine(static_cast<std::uint64_t>(cli.seed) * 1000003ULL +
                    static_cast<std::uint64_t>(G4Threading::G4GetThreadId() + 2)) {}

const std::array<G4double, ensemble::kNTargets>&
EnsembleSteppingAction::Densities(const G4Material* material) {
  auto found = fDensityCache.find(material);
  if (found != fDensityCache.end()) return found->second;
  std::array<G4double, ensemble::kNTargets> densities{};
  const auto* elements = material->GetElementVector();
  const G4double* atoms = material->GetVecNbOfAtomsPerVolume();
  for (std::size_t i = 0; i < material->GetNumberOfElements(); ++i) {
    const int Z = (*elements)[i]->GetZasInt();
    for (int t = 0; t < ensemble::kNTargets; ++t) {
      if (Z == ensemble::kTargetZ[t])
        densities[t] = atoms[i] * ensemble::kTargetAbundance[t] / (1. / cm3);
    }
  }
  return fDensityCache.emplace(material, densities).first->second;
}

void EnsembleSteppingAction::UserSteppingAction(const G4Step* step) {
  const G4Track* track = step->GetTrack();
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto* preVol = pre->GetPhysicalVolume();
  const bool inPhantom =
      preVol && preVol->GetLogicalVolume() == fDet.PhantomLogical();
  auto* run = static_cast<EnsembleRun*>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  const G4double halfZ = fDet.BeamAxisHalfExtent();

  // Steps can be several mm long at high energy, so booking a whole step at
  // its midpoint aliases the depth profile with the step-length period and
  // lumps the step's energy interval into one energy bin. Every tally
  // therefore subdivides the step into segments shorter than half a bin (in
  // depth, and for the exposure also in energy), interpolating position and
  // kinetic energy linearly along the path.
  const G4double preZ = pre->GetPosition().z(), postZ = post->GetPosition().z();
  const G4double preR = pre->GetPosition().perp(),
                 postR = post->GetPosition().perp();
  const G4double preE = pre->GetKineticEnergy(),
                 postE = post->GetKineticEnergy();
  const G4double stepLenMm = step->GetStepLength() / mm;

  auto subdivisions = [&](bool with_energy) {
    double n = std::ceil(stepLenMm / (0.5 * fCli.zbin_width_mm));
    if (with_energy)
      n = std::max(n, std::ceil(std::abs(preE - postE) / MeV /
                                (0.5 * fCli.ebin_width_MeV)));
    return static_cast<int>(std::min(std::max(n, 1.0), 256.0));
  };

  // --- dose tallies -------------------------------------------------------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) {
    if (inPhantom) run->AddEdep(edep);
    const bool primary = track->GetParentID() == 0;
    const int n = subdivisions(false);
    for (int k = 0; k < n; ++k) {
      const G4double f = (k + 0.5) / n;
      const G4double z = preZ + f * (postZ - preZ);
      const G4double r = preR + f * (postR - preR);
      const int bin = static_cast<int>((z + halfZ) / mm / fCli.zbin_width_mm);
      run->AddDepthEdep(bin, edep / MeV / n, primary,
                        r <= ensemble::kCoreRadiusMM * mm);
      if (z >= fDet.TargetProxZ() && z <= fDet.TargetDistZ() &&
          r <= fDet.TargetRadius()) {
        run->AddTargetEdep(edep / n);
      }
    }
  }

  if (!inPhantom) return;

  // --- proton exposure and in-flight production sampling ------------------
  if (track->GetParticleDefinition() == G4Proton::Definition()) {
    const auto& densities = Densities(pre->GetMaterial());
    const G4double weight = track->GetWeight();
    const int n = subdivisions(true);
    const G4double stepLenCm = step->GetStepLength() / cm;
    const auto& prePos = pre->GetPosition();
    const auto& postPos = post->GetPosition();
    for (int k = 0; k < n; ++k) {
      const G4double f = (k + 0.5) / n;
      const G4double energyMeV = (preE + f * (postE - preE)) / MeV;
      const G4double depthMm = (preZ + f * (postZ - preZ) + halfZ) / mm;
      for (int t = 0; t < ensemble::kNTargets; ++t) {
        if (densities[t] > 0.)
          run->Grid().Add(t, energyMeV, depthMm,
                          weight * densities[t] * stepLenCm / n);
      }
      if (!fCurves) continue;
      for (std::size_t c = 0; c < fCurves->Channels().size(); ++c) {
        const auto& channel = fCurves->Channels()[c];
        const G4double density = densities[channel.target_index];
        if (density <= 0.) continue;
        const G4double sigma = channel.SigmaMb(energyMeV);
        if (sigma <= 0.) continue;
        const G4double probability =
            weight * density * (stepLenCm / n) * sigma * 1.0e-27;
        if (fUniform(fSampleEngine) >= probability) continue;
        const G4double u = fUniform(fSampleEngine);
        const G4double g = (k + u) / n;
        const auto position = prePos + g * (postPos - prePos);
        run->AddSampledProduction(
            {G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(),
             static_cast<int>(c), position.x() / mm, position.y() / mm,
             position.z() / mm, energyMeV});
      }
    }
  }

  // --- keep escaping positrons at the phantom surface ---------------------
  // stageA convention: a positron leaving into air is stopped at the
  // boundary, so its end point is the recorded annihilation position.
  if (track->GetParticleDefinition() == G4Positron::Definition()) {
    const auto* postVol = post->GetPhysicalVolume();
    if (postVol == nullptr || postVol->GetMotherLogical() == nullptr) {
      step->GetTrack()->SetTrackStatus(fStopAndKill);
    }
  }

  // --- native beta-plus residual production -------------------------------
  const auto* process = post->GetProcessDefinedStep();
  if (process && process->GetProcessSubType() == fHadronInelastic) {
    const auto* hadronic = dynamic_cast<const G4HadronicProcess*>(process);
    std::string target = "unknown";
    if (hadronic) {
      auto* nucleus = const_cast<G4HadronicProcess*>(hadronic)->GetTargetNucleus();
      if (nucleus) target = NucleusLabel(nucleus->GetZ_asInt(), nucleus->GetA_asInt());
    }
    const int zbin = static_cast<int>(
        (post->GetPosition().z() + halfZ) / mm / fCli.zbin_width_mm);
    for (const auto* secondary : *step->GetSecondaryInCurrentStep()) {
      const auto* definition = secondary->GetParticleDefinition();
      if (definition->GetParticleType() != "nucleus") continue;
      const int Z = definition->GetAtomicNumber();
      const int A = definition->GetAtomicMass();
      for (std::size_t r = 0; r < ensemble::kBetaPlusResiduals.size(); ++r) {
        if (Z == ensemble::kBetaPlusResiduals[r].Z &&
            A == ensemble::kBetaPlusResiduals[r].A) {
          run->CountNative({track->GetParticleDefinition()->GetParticleName(),
                            target, static_cast<int>(r), zbin});
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Run action

EnsembleRunAction::EnsembleRunAction(const ensemble::EnsembleCli& cli,
                                     const EnsembleDetector& det,
                                     const ensemble::SamplingCurves* curves)
    : fCli(cli), fDet(det), fCurves(curves) {}

G4Run* EnsembleRunAction::GenerateRun() {
  return new EnsembleRun(fDet.BeamAxisHalfExtent() / mm, fCli.ebin_width_MeV,
                         fCli.emax_MeV, fCli.zbin_width_mm);
}

void EnsembleRunAction::EndOfRunAction(const G4Run* base_run) {
  if (!IsMaster()) return;
  const auto* run = static_cast<const EnsembleRun*>(base_run);

  // Emitter-transport mode: join the input seeds with the recorded
  // annihilations and write the core emitters.csv next to the input file.
  if (!fCli.emitters_in.empty()) {
    const auto seeds = ensemble::ReadEmitterSeeds(fCli.emitters_in);
    const auto& annihilations = run->Annihilations();
    const std::string dir =
        std::filesystem::path(fCli.emitters_in).parent_path().string();
    std::ofstream f(dir + "/emitters.csv");
    f << "event_id,isotope_id,prod_x_mm,prod_y_mm,prod_z_mm,"
         "anh_x_mm,anh_y_mm,anh_z_mm\n";
    f << std::setprecision(7);
    long missing = 0;
    for (std::size_t i = 0; i < seeds.size(); ++i) {
      const auto found = annihilations.find(static_cast<int>(i));
      if (found == annihilations.end()) { ++missing; continue; }
      const auto& seed = seeds[i];
      f << seed.event_id << ',' << seed.isotope_index << ',' << seed.x_mm
        << ',' << seed.y_mm << ',' << seed.z_mm << ','
        << found->second[0] << ',' << found->second[1] << ','
        << found->second[2] << '\n';
    }
    G4cout << "[ensemble] wrote " << seeds.size() - missing << " emitters ("
           << missing << " without a captured positron) -> " << dir
           << "/emitters.csv" << G4endl;
    return;
  }

  if (run->Grid().Overflow() > 0) {
    G4Exception("EnsembleRunAction", "ExposureOverflow", FatalException,
                ("exposure grid overflow: " +
                 std::to_string(run->Grid().Overflow()) +
                 " step contributions fell outside the energy or depth window")
                    .c_str());
  }

  const std::string beam = fCli.layers.empty() ? "pencil" : "sobp";
  const std::string tag = fCli.geometry + "_" + beam + "_" + fCli.physics_list +
                          "_" + FormatProtonCount(run->GetNumberOfEvent());
  const std::string dir = fCli.output_base + "/" + tag;
  std::filesystem::create_directories(dir);
  G4cout << "[ensemble] run directory -> " << dir << G4endl;

  // Exposure table.
  run->Grid().WriteCsv(dir + "/proton_exposure.csv");

  // Native-route counters.
  {
    std::ofstream f(dir + "/native_route_counts.csv");
    f << "projectile,target,residual,depth_low_mm,depth_high_mm,"
         "production_count\n";
    f << std::setprecision(12);
    for (const auto& [key, count] : run->Native()) {
      const auto& [projectile, target, residual, zbin] = key;
      f << projectile << ',' << target << ','
        << ensemble::kBetaPlusResiduals[residual].name << ','
        << zbin * fCli.zbin_width_mm << ',' << (zbin + 1) * fCli.zbin_width_mm
        << ',' << count << '\n';
    }
  }

  // Depth dose (stageA 5-column format; z is the world coordinate).
  {
    std::ofstream f(dir + "/depth_dose.csv");
    f << "z_mm,edep_total_MeV,edep_primary_MeV,edep_core_MeV,dose_core_Gy\n";
    f << std::setprecision(9);
    const double binw = fCli.zbin_width_mm;
    const double rho =
        fDet.PhantomLogical()->GetMaterial()->GetDensity();
    const double rcore = ensemble::kCoreRadiusMM * mm;
    const double vbin = pi * rcore * rcore * (binw * mm);
    for (int i = 0; i < run->DepthBins(); ++i) {
      const double zc = -run->HalfZmm() + (i + 0.5) * binw;
      const double doseGy =
          ((run->EdepZCore()[i] * MeV) / (rho * vbin)) / gray;
      f << zc << ',' << run->EdepZTotal()[i] << ',' << run->EdepZPrimary()[i]
        << ',' << run->EdepZCore()[i] << ',' << doseGy << '\n';
    }
  }

  // Raw metadata; the Python driver adds the digest and software revision and
  // writes the validated exposure_meta.json.
  {
    const long nProtons = run->GetNumberOfEvent();
    const double targetMass = fDet.TargetMass();
    const double targetDoseGy =
        (targetMass > 0.) ? (run->TargetEdep() / targetMass) / gray : 0.;
    const int major = G4VERSION_NUMBER / 100;
    const int minor = (G4VERSION_NUMBER / 10) % 10;
    const int patch = G4VERSION_NUMBER % 10;

    std::ofstream f(dir + "/run_meta_raw.json");
    f << std::setprecision(17);
    f << "{\n";
    f << "  \"schema_version\": 1,\n";
    f << "  \"run_id\": \"" << tag << "\",\n";
    f << "  \"exposure_file\": \"proton_exposure.csv\",\n";
    f << "  \"n_protons\": " << nProtons << ",\n";
    f << "  \"target_dose_Gy\": " << targetDoseGy << ",\n";
    f << "  \"physics_list\": \"" << fCli.physics_list << "\",\n";
    f << "  \"geant4_version\": \"" << major << '.' << minor << '.' << patch
      << "\",\n";
    f << "  \"random_seed\": " << fCli.seed << ",\n";
    f << "  \"beam_axis\": \"z\",\n";
    f << "  \"depth_origin\": \"phantom entrance face; depth = z_world + "
      << run->HalfZmm() << " mm\",\n";
    f << "  \"depth_unit\": \"mm\",\n";
    f << "  \"geometry\": \"" << fCli.geometry << "\",\n";
    f << "  \"beam\": \"" << beam << "\",\n";
    f << "  \"sobp_layers\": \"" << fCli.layers << "\",\n";
    f << "  \"disk_radius_mm\": " << fCli.disk_radius_mm << ",\n";
    f << "  \"target_radius_mm\": " << fCli.target_radius_mm << ",\n";
    f << "  \"target_prox_mm\": " << fCli.target_prox_mm << ",\n";
    f << "  \"target_dist_mm\": " << fCli.target_dist_mm << ",\n";
    f << "  \"edep_total_MeV\": " << run->EdepTotal() / MeV << ",\n";
    f << "  \"phantom_mass_g\": " << fDet.PhantomMass() / g << ",\n";
    f << "  \"target_mass_g\": " << targetMass / g << ",\n";
    f << "  \"total_exposure_cm2_inv\": " << run->Grid().TotalExposure()
      << ",\n";
    f << "  \"energy_edges_MeV\": [";
    for (int i = 0; i <= run->Grid().EnergyBins(); ++i)
      f << (i ? ", " : "") << run->Grid().EnergyEdge(i);
    f << "],\n";
    f << "  \"depth_edges_mm\": [";
    for (int i = 0; i <= run->Grid().DepthBins(); ++i)
      f << (i ? ", " : "") << run->Grid().DepthEdge(i);
    f << "]\n";
    f << "}\n";

    G4cout << "[ensemble] protons " << nProtons << ", target dose "
           << targetDoseGy << " Gy, total exposure "
           << run->Grid().TotalExposure() << " /cm2, native routes "
           << run->Native().size() << G4endl;
  }

  // Productions sampled in flight from the fitted curves.
  if (fCurves) {
    std::ofstream f(dir + "/sampled_productions.csv");
    f << "event_id,channel_id,target,residual,isotope_id,prod_x_mm,"
         "prod_y_mm,prod_z_mm,proton_energy_MeV\n";
    f << std::setprecision(9);
    for (const auto& production : run->Sampled()) {
      const auto& channel = fCurves->Channels()[production.channel_index];
      f << production.event_id << ',' << channel.channel_id << ','
        << ensemble::kTargetNames[channel.target_index] << ','
        << ensemble::kBetaPlusResiduals[channel.residual_index].name << ','
        << channel.residual_index << ',' << production.x_mm << ','
        << production.y_mm << ',' << production.z_mm << ','
        << production.proton_energy_MeV << '\n';
    }
    std::ofstream meta(dir + "/sampling_meta.json");
    meta << "{\n  \"curves_file\": \"" << fCurves->Path() << "\",\n"
         << "  \"sampled_productions\": " << run->Sampled().size() << "\n}\n";
    G4cout << "[ensemble] sampled " << run->Sampled().size()
           << " productions from " << fCurves->Path() << G4endl;
  }

  // Copy the SOBP layer table (and its provenance meta) into the run dir.
  if (!fCli.layers.empty()) {
    std::error_code ec;
    std::filesystem::copy_file(fCli.layers, dir + "/sobp_layers.csv",
                               std::filesystem::copy_options::overwrite_existing,
                               ec);
    const std::string meta =
        fCli.layers.substr(0, fCli.layers.size() - 4) + "_meta.csv";
    if (std::filesystem::exists(meta, ec))
      std::filesystem::copy_file(meta, dir + "/sobp_layers_meta.csv",
                                 std::filesystem::copy_options::overwrite_existing,
                                 ec);
  }
}

// ---------------------------------------------------------------------------

void EnsembleActionInitialization::BuildForMaster() const {
  SetUserAction(new EnsembleRunAction(fCli, fDet, fCurves));
}

void EnsembleActionInitialization::Build() const {
  if (fCli.emitters_in.empty()) {
    SetUserAction(new EnsemblePrimaryGenerator(fCli, fDet));
  } else {
    SetUserAction(new EmitterPrimaryGenerator(fCli));
    SetUserAction(new EnsembleTrackingAction());
  }
  SetUserAction(new EnsembleRunAction(fCli, fDet, fCurves));
  SetUserAction(new EnsembleSteppingAction(fCli, fDet, fCurves));
}
