#include "SamplingCurves.hh"

#include "EnsembleConfig.hh"

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace ensemble {

namespace {
double Interpolate(const std::vector<double>& grid,
                   const std::vector<double>& values, double threshold,
                   const std::string& label, double energy) {
  if (energy <= threshold || energy < grid.front()) return 0.;
  if (energy > grid.back())
    throw std::runtime_error(label + ": energy " + std::to_string(energy) +
                             " MeV exceeds the fitted range");
  auto upper = std::lower_bound(grid.begin(), grid.end(), energy);
  if (upper == grid.begin()) return values.front();
  const std::size_t hi = upper - grid.begin();
  const std::size_t lo = hi - 1;
  const double fraction = (energy - grid[lo]) / (grid[hi] - grid[lo]);
  return values[lo] + fraction * (values[hi] - values[lo]);
}
}  // namespace

double SamplingChannel::SigmaMb(double energy) const {
  return Interpolate(energy_MeV, sigma_mb, threshold_MeV, channel_id, energy);
}

double SamplingChannel::SigmaEnvMb(double energy) const {
  return Interpolate(energy_MeV, sigma_env_mb, threshold_MeV, channel_id,
                     energy);
}

SamplingCurves SamplingCurves::Load(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open sampling curves " + path);
  std::map<std::string, SamplingChannel> channels;
  std::string line;
  std::getline(f, line);  // header
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string id, target, residual, threshold, energy, sigma, envelope;
    std::getline(ss, id, ',');
    std::getline(ss, target, ',');
    std::getline(ss, residual, ',');
    std::getline(ss, threshold, ',');
    std::getline(ss, energy, ',');
    std::getline(ss, sigma, ',');
    std::getline(ss, envelope, ',');
    auto& channel = channels[id];
    if (channel.channel_id.empty()) {
      channel.channel_id = id;
      channel.threshold_MeV = std::stod(threshold);
      for (int t = 0; t < kNTargets; ++t)
        if (target == kTargetNames[t]) channel.target_index = t;
      for (std::size_t r = 0; r < kBetaPlusResiduals.size(); ++r)
        if (residual == kBetaPlusResiduals[r].name)
          channel.residual_index = static_cast<int>(r);
      if (channel.target_index < 0 || channel.residual_index < 0)
        throw std::runtime_error(id + ": unknown target or residual");
    }
    channel.energy_MeV.push_back(std::stod(energy));
    channel.sigma_mb.push_back(std::stod(sigma));
    channel.sigma_env_mb.push_back(
        envelope.empty() ? std::stod(sigma) : std::stod(envelope));
  }
  SamplingCurves curves;
  curves.fPath = path;
  for (auto& [id, channel] : channels) {
    if (channel.energy_MeV.size() < 2)
      throw std::runtime_error(id + ": curve has fewer than two points");
    for (std::size_t i = 1; i < channel.energy_MeV.size(); ++i)
      if (channel.energy_MeV[i] <= channel.energy_MeV[i - 1])
        throw std::runtime_error(id + ": energy grid is not increasing");
    curves.fChannels.push_back(std::move(channel));
  }
  if (curves.fChannels.empty())
    throw std::runtime_error(path + ": no channels loaded");
  return curves;
}

}  // namespace ensemble
