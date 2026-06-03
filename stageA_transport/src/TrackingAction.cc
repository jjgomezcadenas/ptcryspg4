#include "TrackingAction.hh"

#include "EventAction.hh"
#include "Isotopes.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

namespace {
// Convert a Geant4 position (G4ThreeVector: double, internal units) to the
// stored form -- a float[3] in millimetres. Dividing by `mm` (a G4SystemOfUnits
// constant) extracts the value *expressed in* mm; this is the unit-safe Geant4
// idiom (never read a raw coordinate assuming its unit). The narrowing to float
// matches EmitterData's float32 columns, which mirror the future HDF5 datasets;
// float's ~7 significant digits give sub-micron precision over the 200 mm
// phantom -- far more than needed.
std::array<float, 3> ToMM(const G4ThreeVector& v) {
  return {static_cast<float>(v.x() / mm), static_cast<float>(v.y() / mm),
          static_cast<float>(v.z() / mm)};
}
}  // namespace

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
  const G4ParticleDefinition* def = track->GetParticleDefinition();

  // (1) Listed emitter ion? Tag trackID -> isotope_id. EmitterId only matches
  // our five beta+ emitters, so daughters/targets fall through harmlessly.
  const std::int8_t ionID =
      ptcrysp::EmitterId(def->GetAtomicNumber(), def->GetAtomicMass());
  if (ionID >= 0) {
    fEvent->RegisterEmitterIon(track->GetTrackID(), ionID);
    return;
  }

  // (2) Positron whose parent is a tagged emitter ion: PROD = its vertex.
  if (def == G4Positron::Definition()) {
    std::int8_t iso = -1;
    if (fEvent->LookupEmitterIon(track->GetParentID(), iso)) {
      fCapturing = true;
      fPositronID = track->GetTrackID();
      fIsotopeID = iso;
      fProd = ToMM(track->GetVertexPosition());
    }
  }
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
  // ANH = where this positron ends (annihilation point -> the 511 keV pair).
  if (fCapturing && track->GetTrackID() == fPositronID) {
    const std::array<float, 3> anh =
        ToMM(track->GetStep()->GetPostStepPoint()->GetPosition());
    fEvent->AddEmitter(fIsotopeID, fProd, anh);
    fCapturing = false;
  }
}
