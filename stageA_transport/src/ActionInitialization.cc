#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"

void ActionInitialization::BuildForMaster() const {
  // The master receives the merged run and writes the CSVs; it owns the
  // output messenger and never generates primaries.
  SetUserAction(new RunAction(/*onMaster=*/true, fDet));
}

void ActionInitialization::Build() const {
  SetUserAction(new PrimaryGeneratorAction(fDet, fBeam));

  auto* eventAction = new EventAction;
  SetUserAction(eventAction);
  SetUserAction(new TrackingAction(eventAction));
  SetUserAction(new SteppingAction(fDet));
  SetUserAction(new RunAction(/*onMaster=*/false, fDet));
}
