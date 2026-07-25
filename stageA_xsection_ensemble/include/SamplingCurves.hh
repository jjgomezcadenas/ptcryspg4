#ifndef SAMPLING_CURVES_HH
#define SAMPLING_CURVES_HH

/// \file
/// The fitted excitation functions used by the in-flight production sampler:
/// per channel, a threshold and a linearly interpolated curve exported by
/// analysis_transport/xsections/export_sampling_curves.py.

#include <string>
#include <vector>

namespace ensemble {

struct SamplingChannel {
  std::string channel_id;
  int target_index = -1;    ///< index into kTargetNames
  int residual_index = -1;  ///< index into kBetaPlusResiduals
  double threshold_MeV = 0.;
  std::vector<double> energy_MeV;
  std::vector<double> sigma_mb;
  std::vector<double> sigma_env_mb;  ///< bank sampling envelope

  /// Linear interpolation; zero at and below threshold and below the grid;
  /// throws above the fitted range.
  double SigmaMb(double energy) const;
  double SigmaEnvMb(double energy) const;
};

class SamplingCurves {
 public:
  /// Load sampling_curves.csv; throws on malformed input.
  static SamplingCurves Load(const std::string& path);

  const std::vector<SamplingChannel>& Channels() const { return fChannels; }
  const std::string& Path() const { return fPath; }

 private:
  std::vector<SamplingChannel> fChannels;
  std::string fPath;
};

}  // namespace ensemble

#endif  // SAMPLING_CURVES_HH
