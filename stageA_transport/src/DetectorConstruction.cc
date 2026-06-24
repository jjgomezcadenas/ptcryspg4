#include "DetectorConstruction.hh"

#include "StageAConfig.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalConstants.hh"

DetectorConstruction::DetectorConstruction()
    : fMaterialName(stageA::kPhantomMaterial),
      fGeometry(stageA::kDefaultGeometry),
      fTargetRadius(stageA::kTargetRadiusMM * mm),
      fTargetProxDepth(stageA::kTargetProxDepthMM * mm),
      fTargetDistDepth(stageA::kTargetDistDepthMM * mm) {
  fMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction() { delete fMessenger; }

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  const bool head = (fGeometry == stageA::kGeometryMirdHead);

  // --- world: an air box sized to contain whichever geometry, with a margin.
  // The world is the top volume; tracks die when they leave it.
  G4double worldHalfXY, worldHalfZ;
  if (head) {
    worldHalfXY = worldHalfZ = 30. * cm;  // generous: head bounds ~±20 cm
  } else {
    fRadius = 0.5 * stageA::kPhantomDiameterMM * mm;
    fHalfZ = 0.5 * stageA::kPhantomLengthMM * mm;
    worldHalfXY = fRadius + 20. * cm;
    worldHalfZ = fHalfZ + 20. * cm;
  }
  auto* worldSolid = new G4Box("World", worldHalfXY, worldHalfXY, worldHalfZ);
  auto* worldLV = new G4LogicalVolume(worldSolid, air, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  // The single un-mothered placement (nullptr mother) defines the coordinate
  // origin; the "true" final arg enables overlap checking.
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0, true);

  if (head)
    BuildMirdHead(worldLV);
  else
    BuildCylinder(worldLV);

  return worldPV;
}

// Default Parodi phantom: a homogeneous cylinder centred at the origin, axis +z.
void DetectorConstruction::BuildCylinder(G4LogicalVolume* worldLV) {
  auto* nist = G4NistManager::Instance();
  G4Material* phantomMat = nist->FindOrBuildMaterial(fMaterialName);
  auto* phantomSolid =
      new G4Tubs("Phantom", 0., fRadius, fHalfZ, 0., 360. * deg);
  fPhantomLV = new G4LogicalVolume(phantomSolid, phantomMat, "Phantom");
  fPhantomLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.3)));
  new G4PVPlacement(nullptr, {}, fPhantomLV, "Phantom", worldLV, false, 0, true);
  fBeamHalfExtent = fHalfZ;
}

// Heterogeneous MIRD head (Phase 1): soft-tissue outer head ⊃ bone skull shell ⊃
// brain, replicated from advanced/human_phantom. Built in the head-local frame
// (x = L-R, y = A-P, z = S-I), then placed so the brain centre is at the origin
// and the L-R axis lies along the beam (+z) — a lateral field.
void DetectorConstruction::BuildMirdHead(G4LogicalVolume* worldLV) {
  using namespace stageA;
  auto* nist = G4NistManager::Instance();
  G4Material* scalp = nist->FindOrBuildMaterial(kScalpMaterial);
  G4Material* bone = nist->FindOrBuildMaterial(kSkullMaterial);
  G4Material* brainMat = nist->FindOrBuildMaterial(kBrainMaterial);

  // Outer head: elliptical tube (face/vault) ∪ ellipsoid cap (cranial dome).
  auto* tube = new G4EllipticalTube("HeadTube", kHeadTubeAxMM * mm,
                                    kHeadTubeByMM * mm, kHeadTubeDzMM * mm);
  auto* cap = new G4Ellipsoid("HeadCap", kHeadCapAxMM * mm, kHeadCapByMM * mm,
                              kHeadCapCzMM * mm, 0., kHeadCapCzMM * mm);
  auto* headSolid = new G4UnionSolid("Head", tube, cap, nullptr,
                                     G4ThreeVector(0, 0, kHeadCapDzMM * mm));
  auto* headLV = new G4LogicalVolume(headSolid, scalp, "Head");
  headLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.9, 0.8, 0.7, 0.2)));

  // Lateral placement: rotate local x → world z, then translate the brain
  // centre (head-local (0,0,kBrainPosZ)) to the world origin.
  const G4Transform3D tf =
      G4Translate3D(-kBrainPosZMM * mm, 0., 0.) * G4RotateY3D(90. * deg);
  new G4PVPlacement(tf, headLV, "Head", worldLV, false, 0, true);

  // Skull cranium shell (outer − inner ellipsoid), daughter of the head.
  auto* cOut = new G4Ellipsoid("CraniumOut", kSkullOutAxMM * mm,
                               kSkullOutByMM * mm, kSkullOutCzMM * mm);
  auto* cIn = new G4Ellipsoid("CraniumIn", kSkullInAxMM * mm, kSkullInByMM * mm,
                              kSkullInCzMM * mm);
  auto* skullSolid = new G4SubtractionSolid("Skull", cOut, cIn, nullptr,
                                            G4ThreeVector(0, 0, 10. * mm));
  auto* skullLV = new G4LogicalVolume(skullSolid, bone, "Skull");
  skullLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.95, 0.95, 0.9, 0.4)));
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kSkullPosZMM * mm), skullLV,
                    "Skull", headLV, false, 0, true);

  // Brain ellipsoid, filling the cranium cavity; daughter of the head.
  auto* brainSolid = new G4Ellipsoid("Brain", kBrainAxMM * mm, kBrainByMM * mm,
                                     kBrainCzMM * mm);
  auto* brainLV = new G4LogicalVolume(brainSolid, brainMat, "Brain");
  brainLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.4)));
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kBrainPosZMM * mm), brainLV,
                    "Brain", headLV, false, 0, true);

  fPhantomLV = headLV;  // scoring volume = whole head (GetMass includes daughters)
  fBeamHalfExtent = kHeadTubeAxMM * mm;  // L-R semi-axis, now along +z
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

G4String DetectorConstruction::PhantomLabel() const {
  // What run_meta.csv records as "phantom_material": the NIST name for the
  // cylinder, or a label for the heterogeneous head (per-region medium is Phase 2).
  return (fGeometry == stageA::kGeometryMirdHead) ? G4String("MIRD_head")
                                                  : fMaterialName;
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
