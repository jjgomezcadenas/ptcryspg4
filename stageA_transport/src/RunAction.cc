#include "RunAction.hh"

#include "StageARun.hh"
#include "StageAConfig.hh"
#include "OutputMessenger.hh"
#include "DetectorConstruction.hh"
#include "Isotopes.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <array>
#include <fstream>
#include <iomanip>

RunAction::RunAction(G4bool onMaster, const DetectorConstruction* det)
    : fDet(det) {
  // Only the master owns the messenger, to avoid duplicate UI-command
  // registration from the worker threads.
  if (onMaster) fMessenger = new OutputMessenger(this);
}

RunAction::~RunAction() { delete fMessenger; }

// Called by the run manager (on master and on each worker) to create the run
// object. Returning our subclass is how we attach custom per-run accumulators.
G4Run* RunAction::GenerateRun() { return new StageARun; }

// Called when a run finishes. In MT each worker's EndOfRunAction fires first,
// then G4 merges the worker runs into the master's run and calls this on the
// master — so only the master holds the complete data and writes the output.
void RunAction::EndOfRunAction(const G4Run* run) {
  if (!IsMaster()) return;  // only the merged master run has all the data
  const auto* stRun = static_cast<const StageARun*>(run);

  WriteEmittersCsv(stRun);
  WriteMetaCsv(stRun);
  WriteDepthDoseCsv(stRun);
  PrintSummary(stRun);
}

// One row per captured β⁺ emitter: event_id, isotope, prod point, anh point.
void RunAction::WriteEmittersCsv(const StageARun* run) const {
  const auto& e = run->Emitters();
  const std::string path = fOutputDir + "/emitters.csv";
  std::ofstream f(path);
  if (!f) {
    G4cerr << "[Stage A] ERROR: cannot open " << path << " for writing."
           << G4endl;
    return;
  }
  f << "event_id,isotope_id,prod_x_mm,prod_y_mm,prod_z_mm,"
       "anh_x_mm,anh_y_mm,anh_z_mm\n";
  f << std::setprecision(7);
  for (std::size_t i = 0; i < e.size(); ++i) {
    f << e.event_id[i] << ',' << static_cast<int>(e.isotope_id[i]) << ','
      << e.prod_xyz[i][0] << ',' << e.prod_xyz[i][1] << ',' << e.prod_xyz[i][2]
      << ',' << e.anh_xyz[i][0] << ',' << e.anh_xyz[i][1] << ','
      << e.anh_xyz[i][2] << '\n';
  }
  G4cout << "[Stage A] wrote " << e.size() << " emitters -> " << path << G4endl;
}

// One wide row of run-level metadata (the would-be HDF5 root attributes):
// beam, phantom, dose, provenance. Dose = total edep / phantom mass.
void RunAction::WriteMetaCsv(const StageARun* run) const {
  const G4long nProtons = run->GetNumberOfEvent();
  const G4double edepMeV = run->EdepTotal() / MeV;
  const G4double mass = fDet->PhantomMass();
  const G4double doseGy = (mass > 0.) ? (run->EdepTotal() / mass) / gray : 0.;

  // Target-box dose (the normalization): D_sim = target edep / target mass.
  const G4double tMass = fDet->TargetMass();
  const G4double tDoseGy = (tMass > 0.) ? (run->TargetEdep() / tMass) / gray : 0.;
  const G4double npPerGy = (tDoseGy > 0.) ? nProtons / tDoseGy : 0.;
  const G4double tProxDepth = fDet->TargetProxZ() + fDet->PhantomHalfLength();
  const G4double tDistDepth = fDet->TargetDistZ() + fDet->PhantomHalfLength();

  const int major = G4VERSION_NUMBER / 100;
  const int minor = (G4VERSION_NUMBER / 10) % 10;
  const int patch = G4VERSION_NUMBER % 10;

  const std::string path = fOutputDir + "/run_meta.csv";
  std::ofstream f(path);
  if (!f) {
    G4cerr << "[Stage A] ERROR: cannot open " << path << " for writing."
           << G4endl;
    return;
  }
  f << "n_protons,beam_energy_MeV,beam_sigma_mm,phantom_material,"
       "phantom_diameter_mm,phantom_length_mm,phantom_mass_g,edep_total_MeV,"
       "dose_total_Gy,target_dose_Gy,target_mass_g,target_radius_mm,"
       "target_prox_depth_mm,target_dist_depth_mm,Np_per_Gy,"
       "geant4_version,physics_list,random_seed\n";
  f << std::setprecision(7) << nProtons << ',' << stageA::kBeamEnergyMeV << ','
    << stageA::kBeamSigmaMM << ',' << fDet->PhantomLabel() << ','
    << stageA::kPhantomDiameterMM << ',' << stageA::kPhantomLengthMM << ','
    << (mass / g) << ',' << edepMeV << ',' << doseGy << ',' << tDoseGy << ','
    << (tMass / g) << ',' << (fDet->TargetRadius() / mm) << ','
    << (tProxDepth / mm) << ',' << (tDistDepth / mm) << ',' << npPerGy << ','
    << major << '.' << minor << '.' << patch << ',' << stageA::kPhysicsList
    << ',' << G4Random::getTheSeed() << '\n';
  G4cout << "[Stage A] wrote run metadata -> " << path << G4endl;
}

// The Bragg profile: energy deposited per 1 mm z-bin, total and primary-only,
// with the z column as bin centres.
void RunAction::WriteDepthDoseCsv(const StageARun* run) const {
  const auto& tot = run->EdepZTotal();
  const auto& prim = run->EdepZPrimary();
  const double zmin = -0.5 * stageA::kPhantomLengthMM;
  const double binw = stageA::kPhantomLengthMM / StageARun::kNZBins;

  const std::string path = fOutputDir + "/depth_dose.csv";
  std::ofstream f(path);
  if (!f) {
    G4cerr << "[Stage A] ERROR: cannot open " << path << " for writing."
           << G4endl;
    return;
  }
  f << "z_mm,edep_total_MeV,edep_primary_MeV\n";
  f << std::setprecision(7);
  for (int i = 0; i < StageARun::kNZBins; ++i) {
    const double zc = zmin + (i + 0.5) * binw;  // bin centre
    f << zc << ',' << tot[i] << ',' << prim[i] << '\n';
  }
  G4cout << "[Stage A] wrote depth-dose profile -> " << path << G4endl;
}

// Human-readable stdout digest for an at-a-glance sanity check (yields per
// isotope, dose, primary/secondary split); the real analysis is the Python.
void RunAction::PrintSummary(const StageARun* run) const {
  const auto& e = run->Emitters();
  const G4long nProtons = run->GetNumberOfEvent();
  const G4double edepMeV = run->EdepTotal() / MeV;
  const G4double mass = fDet->PhantomMass();
  const G4double doseGy = (mass > 0.) ? (run->EdepTotal() / mass) / gray : 0.;

  std::array<long, ptcrysp::kNEmitters> counts{};
  for (auto id : e.isotope_id) {
    if (id >= 0 && id < ptcrysp::kNEmitters) ++counts[id];
  }

  G4cout << "\n========== Stage A summary ==========\n"
         << "  protons simulated : " << nProtons << "\n"
         << "  emitters captured : " << e.size();
  if (nProtons > 0) {
    G4cout << "   (" << static_cast<double>(e.size()) / nProtons
           << " /proton)";
  }
  G4cout << "\n  per isotope:\n";
  for (int id = 0; id < ptcrysp::kNEmitters; ++id) {
    const double perP = nProtons > 0 ? static_cast<double>(counts[id]) / nProtons
                                     : 0.;
    G4cout << "    " << std::setw(4) << ptcrysp::EmitterName(id) << " : "
           << std::setw(8) << counts[id] << "   (" << perP << " /proton)\n";
  }
  // Primary vs secondary share of the deposited energy (from the depth tally).
  double sumTot = 0., sumPrim = 0.;
  for (auto v : run->EdepZTotal()) sumTot += v;
  for (auto v : run->EdepZPrimary()) sumPrim += v;
  const double primFrac = sumTot > 0. ? sumPrim / sumTot : 0.;

  G4cout << "  edep total        : " << edepMeV << " MeV";
  if (nProtons > 0) G4cout << "   (" << edepMeV / nProtons << " MeV/proton)";
  const G4double tMass = fDet->TargetMass();
  const G4double tDoseGy = (tMass > 0.) ? (run->TargetEdep() / tMass) / gray : 0.;
  const G4double npPerGy = (tDoseGy > 0.) ? nProtons / tDoseGy : 0.;

  G4cout << "\n  dose (whole phantom): " << doseGy << " Gy";
  if (nProtons > 0) G4cout << "   (" << doseGy / nProtons << " Gy/proton)";
  G4cout << "\n  dose (target box)   : " << tDoseGy << " Gy   -> N_p(1 Gy) = "
         << npPerGy;
  G4cout << "\n  primary proton edep : " << 100. * primFrac
         << " %   (secondaries " << 100. * (1. - primFrac) << " %)";
  G4cout << "\n=====================================\n" << G4endl;
}
