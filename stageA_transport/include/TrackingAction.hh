#ifndef STAGEA_TRACKINGACTION_HH
#define STAGEA_TRACKINGACTION_HH

#include "G4UserTrackingAction.hh"
#include "globals.hh"

#include <array>
#include <cstdint>

class EventAction;

// Positron-centric capture (spec sec 2.4):
//  - PreUserTracking: tag listed emitter ions (Z,A -> isotope_id) in EventAction;
//    when a positron from such an ion starts, note PROD = its vertex.
//  - PostUserTracking: when that positron ends, ANH = its end point -> emit row.
class TrackingAction : public G4UserTrackingAction {
 public:
  explicit TrackingAction(EventAction* eventAction) : fEvent(eventAction) {}
  ~TrackingAction() override = default;

  void PreUserTrackingAction(const G4Track* track) override;
  void PostUserTrackingAction(const G4Track* track) override;

 private:
  EventAction* fEvent = nullptr;

  // State for the positron currently being tracked (tracking is depth-first, so
  // PreUserTracking and PostUserTracking of one track bracket it contiguously).
  bool fCapturing = false;
  G4int fPositronID = -1;
  std::int8_t fIsotopeID = -1;
  std::array<float, 3> fProd{};
};

#endif  // STAGEA_TRACKINGACTION_HH
