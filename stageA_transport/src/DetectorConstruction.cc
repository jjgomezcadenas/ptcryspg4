#include "DetectorConstruction.hh"

#include "StageAConfig.hh"
#include "DetectorMessenger.hh"

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
#include "G4PhysicalConstants.hh"

DetectorConstruction::DetectorConstruction()
    : fMaterialName(stageA::kPhantomMaterial),
      fTargetRadius(stageA::kTargetRadiusMM * mm),
      fTargetProxDepth(stageA::kTargetProxDepthMM * mm),
      fTargetDistDepth(stageA::kTargetDistDepthMM * mm) {
  fMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction() { delete fMessenger; }

G4VPhysicalVolume* DetectorConstruction::Construct() {
  // NIST manager builds standard materials (correct composition + density) by
  // name, so we never assemble elements/isotopes by hand.
  auto* nist = G4NistManager::Instance();
  G4Material* phantomMat = nist->FindOrBuildMaterial(fMaterialName);
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  fRadius = 0.5 * stageA::kPhantomDiameterMM * mm;
  fHalfZ = 0.5 * stageA::kPhantomLengthMM * mm;

  // --- world: an air box with a margin around the phantom -------------------
  // The world is the top volume that contains everything; tracks die when they
  // leave it.
  const G4double worldXY = fRadius + 20. * cm;
  const G4double worldHalfZ = fHalfZ + 20. * cm;
  auto* worldSolid = new G4Box("World", worldXY, worldXY, worldHalfZ);
  auto* worldLV = new G4LogicalVolume(worldSolid, air, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  // The single un-mothered placement (nullptr mother) defines the coordinate
  // origin; the "true" final arg enables overlap checking.
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0, true);

  // --- phantom: cylinder centred at the origin, axis along z ----------------
  auto* phantomSolid =
      new G4Tubs("Phantom", 0., fRadius, fHalfZ, 0., 360. * deg);
  fPhantomLV = new G4LogicalVolume(phantomSolid, phantomMat, "Phantom");
  fPhantomLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.3)));
  new G4PVPlacement(nullptr, {}, fPhantomLV, "Phantom", worldLV, false, 0, true);

  return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
  // A multifunctional detector is a container of "primitive scorers"; here a
  // single G4PSEnergyDeposit sums the energy deposited in the phantom each
  // event. G4 merges its hits map across threads automatically.
  auto* mfd = new G4MultiFunctionalDetector(stageA::kScorerMFD);
  G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
  mfd->RegisterPrimitive(new G4PSEnergyDeposit(stageA::kScorerEdep));
  SetSensitiveDetector(fPhantomLV, mfd);  // bind the scorer to the phantom volume
}

G4double DetectorConstruction::PhantomMass() const {
  // G4LogicalVolume::GetMass() integrates ρ·V over the volume (and daughters).
  return fPhantomLV ? fPhantomLV->GetMass() : 0.;
}

G4double DetectorConstruction::TargetMass() const {
  // The target box is conceptual (not a volume), so compute its mass directly
  // from its dimensions and the phantom material density.
  if (!fPhantomLV) return 0.;
  const G4double length = fTargetDistDepth - fTargetProxDepth;
  const G4double volume = pi * fTargetRadius * fTargetRadius * length;
  return volume * fPhantomLV->GetMaterial()->GetDensity();
}
