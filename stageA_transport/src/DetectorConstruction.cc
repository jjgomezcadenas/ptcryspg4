#include "DetectorConstruction.hh"

#include "StageAConfig.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  // PMMA. NIST G4_PLEXIGLASS is C5H8O2 at rho = 1.19 g/cm^3; the spec quotes
  // 1.18 g/cm^3. We use the NIST material as-is for now (0.8% density
  // difference, negligible); revisit if absolute yields matter.
  G4Material* pmma = nist->FindOrBuildMaterial(stageA::kPhantomMaterial);
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  fRadius = 0.5 * stageA::kPhantomDiameterMM * mm;
  fHalfZ = 0.5 * stageA::kPhantomLengthMM * mm;

  // --- world: air box with a margin around the phantom -----------------------
  const G4double worldXY = fRadius + 20. * cm;
  const G4double worldHalfZ = fHalfZ + 20. * cm;
  auto* worldSolid = new G4Box("World", worldXY, worldXY, worldHalfZ);
  auto* worldLV = new G4LogicalVolume(worldSolid, air, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0, true);

  // --- phantom: PMMA cylinder centred at the origin, axis along z -----------
  auto* phantomSolid =
      new G4Tubs("Phantom", 0., fRadius, fHalfZ, 0., 360. * deg);
  fPhantomLV = new G4LogicalVolume(phantomSolid, pmma, "Phantom");
  fPhantomLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.3)));
  new G4PVPlacement(nullptr, {}, fPhantomLV, "Phantom", worldLV, false, 0, true);

  return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
  // Total energy deposited in the phantom (-> dose). MT-safe: G4 merges the
  // scorer hits map across threads automatically.
  auto* mfd = new G4MultiFunctionalDetector(stageA::kScorerMFD);
  G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
  mfd->RegisterPrimitive(new G4PSEnergyDeposit(stageA::kScorerEdep));
  SetSensitiveDetector(fPhantomLV, mfd);
}

G4double DetectorConstruction::PhantomMass() const {
  return fPhantomLV ? fPhantomLV->GetMass() : 0.;
}
