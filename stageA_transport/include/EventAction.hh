#ifndef STAGEA_EVENTACTION_HH
#define STAGEA_EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "EmitterData.hh"

#include <array>
#include <cstdint>
#include <unordered_map>

// Per-event bookkeeping shared with TrackingAction:
//  - fIonMap links an emitter-ion trackID to its isotope_id (filled when the
//    listed residual nucleus appears), so its positron can be labelled.
//  - fBuffer collects this event's captured {event_id, isotope, prod, anh} rows.
// At end of event the buffer is appended to the thread's StageARun.
class EventAction : public G4UserEventAction {
 public:
  EventAction() = default;
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  G4long EventID() const { return fEventID; }

  void RegisterEmitterIon(G4int trackID, std::int8_t isotopeID) {
    fIonMap[trackID] = isotopeID;
  }
  // If trackID is a registered emitter ion, return true and set isotopeID.
  bool LookupEmitterIon(G4int trackID, std::int8_t& isotopeID) const {
    auto it = fIonMap.find(trackID);
    if (it == fIonMap.end()) return false;
    isotopeID = it->second;
    return true;
  }

  void AddEmitter(std::int8_t isotopeID, const std::array<float, 3>& prod,
                  const std::array<float, 3>& anh) {
    fBuffer.add(fEventID, isotopeID, prod, anh);
  }

 private:
  G4long fEventID = 0;
  std::unordered_map<G4int, std::int8_t> fIonMap;
  ptcrysp::EmitterData fBuffer;
};

#endif  // STAGEA_EVENTACTION_HH
