#include "EventAction.hh"

#include "StageARun.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

// Called before each event's tracking: reset the per-event scratch state
// (the ion→isotope map and the captured-emitter buffer) and record the id.
void EventAction::BeginOfEventAction(const G4Event* event) {
  fEventID = event->GetEventID();
  fIonMap.clear();
  fBuffer.clear();
}

// Called after the event's tracking completes: hand this event's captured
// emitters to the thread's run object for the end-of-run merge.
void EventAction::EndOfEventAction(const G4Event*) {
  if (fBuffer.size() == 0) return;
  // Append to this thread's run (same thread, no lock needed).
  auto* run = static_cast<StageARun*>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->FillEmitters(fBuffer);
}
