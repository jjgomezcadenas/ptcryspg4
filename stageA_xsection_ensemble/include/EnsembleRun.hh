#ifndef ENSEMBLE_RUN_HH
#define ENSEMBLE_RUN_HH

#include "EnsembleConfig.hh"
#include "ExposureGrid.hh"

#include "G4Run.hh"
#include "globals.hh"

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace ensemble {

/// Native-route key: projectile name, target label (element symbol + A),
/// beta-plus residual index into kBetaPlusResiduals, depth bin.
using NativeKey = std::tuple<std::string, std::string, int, int>;

/// One production sampled in flight from the fitted curves.
struct SampledProduction {
  int event_id = 0;
  int channel_index = 0;  ///< index into the loaded SamplingCurves channels
  double x_mm = 0., y_mm = 0., z_mm = 0.;
  double proton_energy_MeV = 0.;
};

}  // namespace ensemble

/// Per-run accumulators: the exposure grid, the native-route counters, and the
/// dose tallies (total, target box, depth profile). Workers fill their own
/// copy; Merge folds them into the master run.
class EnsembleRun : public G4Run {
 public:
  EnsembleRun(double halfZ_mm, double ebin_width_MeV, double emax_MeV,
              double zbin_width_mm);

  void Merge(const G4Run* other) override;

  ensemble::ExposureGrid& Grid() { return fGrid; }
  const ensemble::ExposureGrid& Grid() const { return fGrid; }

  void CountNative(const ensemble::NativeKey& key) { ++fNative[key]; }
  const std::map<ensemble::NativeKey, long long>& Native() const {
    return fNative;
  }

  void AddSampledProduction(const ensemble::SampledProduction& production) {
    fSampled.push_back(production);
  }
  const std::vector<ensemble::SampledProduction>& Sampled() const {
    return fSampled;
  }

  void RecordAnnihilation(int event_id, double x_mm, double y_mm, double z_mm) {
    fAnnihilations[event_id] = {x_mm, y_mm, z_mm};
  }
  const std::map<int, std::array<double, 3>>& Annihilations() const {
    return fAnnihilations;
  }

  void AddEdep(double edep) { fEdepTotal += edep; }
  void AddTargetEdep(double edep) { fTargetEdep += edep; }
  void AddDepthEdep(int bin, double edep_MeV, bool primary, bool core);

  double EdepTotal() const { return fEdepTotal; }
  double TargetEdep() const { return fTargetEdep; }
  double HalfZmm() const { return fHalfZmm; }
  int DepthBins() const { return static_cast<int>(fEdepZTotal.size()); }
  const std::vector<double>& EdepZTotal() const { return fEdepZTotal; }
  const std::vector<double>& EdepZPrimary() const { return fEdepZPrimary; }
  const std::vector<double>& EdepZCore() const { return fEdepZCore; }

 private:
  ensemble::ExposureGrid fGrid;
  std::map<ensemble::NativeKey, long long> fNative;
  std::vector<ensemble::SampledProduction> fSampled;
  std::map<int, std::array<double, 3>> fAnnihilations;
  double fEdepTotal = 0.;   // G4 units
  double fTargetEdep = 0.;  // G4 units
  double fHalfZmm = 0.;
  std::vector<double> fEdepZTotal, fEdepZPrimary, fEdepZCore;  // MeV per bin
};

#endif  // ENSEMBLE_RUN_HH
