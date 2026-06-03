#include "EventAction.hh"

#include "StageARun.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

void EventAction::BeginOfEventAction(const G4Event* event) {
  fEventID = event->GetEventID();
  fIonMap.clear();
  fBuffer.clear();
}

void EventAction::EndOfEventAction(const G4Event*) {
  if (fBuffer.size() == 0) return;
  // Append to this thread's run (same thread, no lock needed).
  auto* run = static_cast<StageARun*>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->FillEmitters(fBuffer);
}
