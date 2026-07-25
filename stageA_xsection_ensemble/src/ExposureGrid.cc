#include "ExposureGrid.hh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace ensemble {

void ExposureGrid::Init(double ebin_width_MeV, double emax_MeV,
                        double zbin_width_mm, double zmax_mm) {
  if (ebin_width_MeV <= 0. || emax_MeV <= 0. || zbin_width_mm <= 0. ||
      zmax_mm <= 0.)
    throw std::runtime_error("ExposureGrid: bin widths and maxima must be positive");
  fEWidth = ebin_width_MeV;
  fZWidth = zbin_width_mm;
  fNE = static_cast<int>(std::ceil(emax_MeV / ebin_width_MeV - 1.0e-9));
  fNZ = static_cast<int>(std::ceil(zmax_mm / zbin_width_mm - 1.0e-9));
  for (int t = 0; t < kNTargets; ++t) {
    fExposure[t].assign(static_cast<std::size_t>(fNE) * fNZ, 0.0);
    fEMoment[t].assign(static_cast<std::size_t>(fNE) * fNZ, 0.0);
    fZMoment[t].assign(static_cast<std::size_t>(fNE) * fNZ, 0.0);
  }
  fOverflow = 0;
}

void ExposureGrid::Add(int t, double energy_MeV, double depth_mm,
                       double contribution) {
  if (contribution <= 0.) return;
  const int ie = static_cast<int>(energy_MeV / fEWidth);
  const int iz = static_cast<int>(depth_mm / fZWidth);
  if (energy_MeV < 0. || ie >= fNE || depth_mm < 0. || iz >= fNZ) {
    ++fOverflow;
    return;
  }
  const std::size_t k = static_cast<std::size_t>(ie) * fNZ + iz;
  fExposure[t][k] += contribution;
  fEMoment[t][k] += contribution * energy_MeV;
  fZMoment[t][k] += contribution * depth_mm;
}

void ExposureGrid::MergeFrom(const ExposureGrid& other) {
  for (int t = 0; t < kNTargets; ++t) {
    for (std::size_t k = 0; k < fExposure[t].size(); ++k) {
      fExposure[t][k] += other.fExposure[t][k];
      fEMoment[t][k] += other.fEMoment[t][k];
      fZMoment[t][k] += other.fZMoment[t][k];
    }
  }
  fOverflow += other.fOverflow;
}

double ExposureGrid::TotalExposure() const {
  double total = 0.;
  for (int t = 0; t < kNTargets; ++t)
    for (double v : fExposure[t]) total += v;
  return total;
}

void ExposureGrid::WriteCsv(const std::string& path) const {
  std::filesystem::path out(path);
  if (out.has_parent_path()) std::filesystem::create_directories(out.parent_path());
  std::ofstream f(out);
  if (!f) throw std::runtime_error("cannot open " + path + " for writing");
  f << "target,energy_low_MeV,energy_high_MeV,energy_mean_MeV,"
       "depth_low_mm,depth_high_mm,depth_mean_mm,target_exposure_cm2_inv\n";
  f << std::setprecision(12);
  for (int t = 0; t < kNTargets; ++t) {
    for (int ie = 0; ie < fNE; ++ie) {
      for (int iz = 0; iz < fNZ; ++iz) {
        const std::size_t k = static_cast<std::size_t>(ie) * fNZ + iz;
        const double exposure = fExposure[t][k];
        if (exposure <= 0.) continue;
        f << kTargetNames[t] << ',' << ie * fEWidth << ',' << (ie + 1) * fEWidth
          << ',' << fEMoment[t][k] / exposure << ',' << iz * fZWidth << ','
          << (iz + 1) * fZWidth << ',' << fZMoment[t][k] / exposure << ','
          << exposure << '\n';
      }
    }
  }
}

}  // namespace ensemble
