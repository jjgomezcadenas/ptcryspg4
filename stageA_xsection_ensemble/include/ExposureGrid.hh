#ifndef EXPOSURE_GRID_HH
#define EXPOSURE_GRID_HH

/// \file
/// The target-resolved proton exposure histogram: for every (target isotope,
/// energy bin, depth bin), the sum of weight x nuclei/cm3 x step length (cm),
/// with the exposure-weighted energy and depth moments that define the bin
/// means. One accumulation routine serves both the Geant4 stepping action and
/// the offline step-table mode, so the Python reference accumulator can be
/// compared against exactly the code that runs during transport.

#include "EnsembleConfig.hh"

#include <string>
#include <vector>

namespace ensemble {

class ExposureGrid {
 public:
  void Init(double ebin_width_MeV, double emax_MeV, double zbin_width_mm,
            double zmax_mm);

  /// Add one step contribution for target index t (0..2). Energy in MeV,
  /// depth in mm, contribution in 1/cm2. Out-of-range energies or depths are
  /// counted as overflow; the run refuses to write a table with overflow.
  void Add(int t, double energy_MeV, double depth_mm, double contribution);

  void MergeFrom(const ExposureGrid& other);
  void WriteCsv(const std::string& path) const;

  long long Overflow() const { return fOverflow; }
  int EnergyBins() const { return fNE; }
  int DepthBins() const { return fNZ; }
  double EnergyEdge(int i) const { return i * fEWidth; }
  double DepthEdge(int i) const { return i * fZWidth; }
  double TotalExposure() const;

 private:
  int fNE = 0, fNZ = 0;
  double fEWidth = 0., fZWidth = 0.;
  // Flat [target][energy*nZ + depth] arrays.
  std::array<std::vector<double>, kNTargets> fExposure, fEMoment, fZMoment;
  long long fOverflow = 0;
};

}  // namespace ensemble

#endif  // EXPOSURE_GRID_HH
