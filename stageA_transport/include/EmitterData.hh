#ifndef STAGEA_EMITTERDATA_HH
#define STAGEA_EMITTERDATA_HH

#include <array>
#include <cstdint>
#include <vector>

namespace ptcrysp {

// Struct-of-arrays holding the captured beta+ emitters. The column layout
// mirrors the future prod_anh.h5 datasets exactly (isotope_id int8,
// prod_xyz [N,3] float32, anh_xyz [N,3] float32), so milestone 3 swaps the CSV
// writer for HighFive with no change to this model. event_id is a milestone-2
// diagnostic (production multiplicity per primary). All positions in mm.
struct EmitterData {
  std::vector<std::int64_t> event_id;
  std::vector<std::int8_t> isotope_id;
  std::vector<std::array<float, 3>> prod_xyz;
  std::vector<std::array<float, 3>> anh_xyz;

  std::size_t size() const { return isotope_id.size(); }

  void clear() {
    event_id.clear();
    isotope_id.clear();
    prod_xyz.clear();
    anh_xyz.clear();
  }

  void add(std::int64_t ev, std::int8_t id, const std::array<float, 3>& prod,
           const std::array<float, 3>& anh) {
    event_id.push_back(ev);
    isotope_id.push_back(id);
    prod_xyz.push_back(prod);
    anh_xyz.push_back(anh);
  }

  // Concatenate another block (used by the MT run-merge on the master).
  void append(const EmitterData& o) {
    event_id.insert(event_id.end(), o.event_id.begin(), o.event_id.end());
    isotope_id.insert(isotope_id.end(), o.isotope_id.begin(), o.isotope_id.end());
    prod_xyz.insert(prod_xyz.end(), o.prod_xyz.begin(), o.prod_xyz.end());
    anh_xyz.insert(anh_xyz.end(), o.anh_xyz.begin(), o.anh_xyz.end());
  }
};

}  // namespace ptcrysp

#endif  // STAGEA_EMITTERDATA_HH
