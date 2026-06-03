#ifndef STAGEA_STAGEARUN_HH
#define STAGEA_STAGEARUN_HH

#include "G4Run.hh"
#include "globals.hh"

#include "EmitterData.hh"

#include <vector>

// Custom run object. Each worker thread fills its own instance (emitter rows via
// EventAction, dose via RecordEvent, depth-dose via SteppingAction); Geant4
// calls Merge() on the master to concatenate/sum all threads -- no manual
// locking. The master's merged instance is what RunAction writes to CSV.
class StageARun : public G4Run {
 public:
  static constexpr int kNZBins = 200;  // 1 mm bins over the 200 mm phantom

  StageARun();
  ~StageARun() override = default;

  void RecordEvent(const G4Event* event) override;  // accumulate total dose
  void Merge(const G4Run* aRun) override;            // concatenate/sum threads

  // Append one event's captured emitters (called by EventAction, same thread).
  void FillEmitters(const ptcrysp::EmitterData& buf) { fEmitters.append(buf); }

  // Bin a step's energy deposit by depth, distributing it across the z-bins the
  // step spans (weighted by path-length overlap) so long plateau steps give a
  // smooth profile rather than midpoint spikes. z in mm, edep in MeV.
  void AddEdepAlongStep(G4double z1_mm, G4double z2_mm, G4double edep_MeV,
                        bool primary);

  const ptcrysp::EmitterData& Emitters() const { return fEmitters; }
  G4double EdepTotal() const { return fEdep; }  // G4 internal energy units
  const std::vector<G4double>& EdepZTotal() const { return fEdepZTotal; }
  const std::vector<G4double>& EdepZPrimary() const { return fEdepZPrimary; }

 private:
  ptcrysp::EmitterData fEmitters;
  G4double fEdep = 0.;
  G4int fCollID = -1;  // cached scorer collection ID
  std::vector<G4double> fEdepZTotal{std::vector<G4double>(kNZBins, 0.)};
  std::vector<G4double> fEdepZPrimary{std::vector<G4double>(kNZBins, 0.)};
};

#endif  // STAGEA_STAGEARUN_HH
