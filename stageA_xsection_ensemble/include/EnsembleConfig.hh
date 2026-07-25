#ifndef ENSEMBLE_CONFIG_HH
#define ENSEMBLE_CONFIG_HH

/// \file
/// Configuration of the exposure application: the command-line record and the
/// geometry/beam constants. The constants replicate the stageA_transport
/// values; agreement of the two applications in dose and depth dose is an
/// explicit validation (workshop/xsections_phases.md, phase 1a).

#include <array>
#include <string>

namespace ensemble {

/// Geometry names (subset of stageA_transport; heterogeneous heads arrive with
/// a later phase).
inline constexpr const char* kGeometryCylinder = "cylinder";
inline constexpr const char* kGeometryUniformHead = "uniform_head";
inline constexpr const char* kGeometryUniformHeadEP = "uniform_headep";

/// Head envelope (MIRD scalp ellipsoid) and placement constants, as in
/// stageA_transport/include/StageAConfig.hh.
inline constexpr double kScalpAxMM = 72.0;
inline constexpr double kScalpByMM = 102.0;
inline constexpr double kScalpCzMM = 87.0;
inline constexpr double kBrainOffsetZMM = 10.0;
inline constexpr double kTumourPosXMM = 0.0;
inline constexpr double kTumourPosZMM = -30.0;
inline constexpr const char* kBrainMaterial = "G4_BRAIN_ICRP";

/// Central-axis depth-dose core radius (stageA convention).
inline constexpr double kCoreRadiusMM = 5.0;

/// Target isotopes of the folded channels, indexed 0..2.
inline constexpr int kNTargets = 3;
inline constexpr std::array<const char*, kNTargets> kTargetNames{
    "C12", "N14", "O16"};
/// Natural isotopic abundances of the target isotopes within their elements.
inline constexpr std::array<double, kNTargets> kTargetAbundance{
    0.9893, 0.99636, 0.99757};
/// Element Z of each target isotope.
inline constexpr std::array<int, kNTargets> kTargetZ{6, 7, 8};

/// Beta-plus residuals tallied by the native-route counters (Z, A, name),
/// matching common/isotopes.py.
struct BetaPlusResidual {
  int Z;
  int A;
  const char* name;
};
inline constexpr std::array<BetaPlusResidual, 5> kBetaPlusResiduals{{
    {8, 15, "O15"}, {6, 11, "C11"}, {7, 13, "N13"},
    {6, 10, "C10"}, {8, 14, "O14"},
}};

/// Command-line record. Defaults are the uniform_headep standard case.
struct EnsembleCli {
  std::string geometry = kGeometryUniformHeadEP;
  std::string layers;                    ///< SOBP layer CSV; empty = pencil
  double disk_radius_mm = 20.0;          ///< lateral fluence disk (SOBP mode)
  double energy_MeV = 100.0;             ///< pencil energy
  double beam_sigma_mm = 3.0;            ///< pencil transverse sigma
  double target_radius_mm = 20.0;
  double target_prox_mm = 57.0;          ///< depth from entrance face
  double target_dist_mm = 97.0;
  std::string material = kBrainMaterial; ///< cylinder medium
  double phantom_diameter_mm = 160.0;    ///< cylinder only
  double phantom_length_mm = 160.0;      ///< cylinder only
  long protons = 0;
  long seed = 12345;
  int threads = 18;
  std::string physics_list = "QGSP_BIC_HP";
  double ebin_width_MeV = 0.5;
  double emax_MeV = 130.0;
  double zbin_width_mm = 1.0;
  std::string output_base = "../../data/xsections/exposure";
  std::string steps_csv;                 ///< offline accumulation input
  std::string steps_out;                 ///< offline accumulation output
  double zmax_mm = 0.0;                  ///< offline depth-grid maximum
};

}  // namespace ensemble

#endif  // ENSEMBLE_CONFIG_HH
