#include "DetectorConstruction.hh"

#include "StageAConfig.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
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

namespace {
// Map a head-local ellipsoid (centre, semi-axes) to a world-frame PhantomRegion.
// The lateral placement world = (z_loc - kBrainOffsetZ, y_loc, -x_loc) is a 90°
// rotation, so axis-aligned stays axis-aligned: world centre = (cz - off, cy,
// -cx), world semi-axes = (sz, sy, sx). Euler angles 0.
PhantomRegion HeadRegion(const char* name, const char* mat, G4double cx,
                         G4double cy, G4double cz, G4double sx, G4double sy,
                         G4double sz) {
  const G4double off = stageA::kBrainOffsetZMM;
  return PhantomRegion{name, mat, "ellipsoid", sz, sy, sx,
                       cz - off, cy, (cx == 0. ? 0. : -cx), 0., 0., 0.};
}
}  // namespace

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

  const bool head = (fGeometry == stageA::kGeometryMirdHead ||
                     fGeometry == stageA::kGeometryUniformHead);

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

  if (fGeometry == stageA::kGeometryMirdHead)
    BuildMirdHead(worldLV);
  else if (fGeometry == stageA::kGeometryUniformHead)
    BuildUniformHead(worldLV);
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

  // One homogeneous region (world frame): a cylinder centred at the origin.
  fRegions = {{"phantom", fMaterialName, "cylinder", fRadius / mm, fRadius / mm,
               fHalfZ / mm, 0., 0., 0., 0., 0., 0.}};
}

// Heterogeneous head (Phase 1): a soft-tissue scalp ellipsoid ⊃ bone skull shell
// ⊃ brain (skull/brain are the MIRD ellipsoids; the MIRD face/neck is dropped —
// the lateral field never crosses it). Built in the head-local frame (x = L-R,
// y = A-P, z = S-I) with the origin at the skull/scalp centre; then placed so the
// brain centre is at the world origin and the L-R axis lies along the beam (+z).
void DetectorConstruction::BuildMirdHead(G4LogicalVolume* worldLV) {
  using namespace stageA;
  auto* nist = G4NistManager::Instance();
  G4Material* scalp = nist->FindOrBuildMaterial(kScalpMaterial);
  G4Material* bone = nist->FindOrBuildMaterial(kSkullMaterial);
  G4Material* brainMat = nist->FindOrBuildMaterial(kBrainMaterial);

  // Outer head: a single soft-tissue scalp ellipsoid enclosing the skull.
  auto* headSolid = new G4Ellipsoid("Head", kScalpAxMM * mm, kScalpByMM * mm,
                                    kScalpCzMM * mm);
  auto* headLV = new G4LogicalVolume(headSolid, scalp, "Head");
  headLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.9, 0.8, 0.7, 0.2)));

  // Lateral placement: rotate local x → world z, then translate the brain centre
  // (head-local (0,0,kBrainOffsetZ)) to the world origin.
  const G4Transform3D tf =
      G4Translate3D(-kBrainOffsetZMM * mm, 0., 0.) * G4RotateY3D(90. * deg);
  new G4PVPlacement(tf, headLV, "Head", worldLV, false, 0, true);

  // Skull cranium shell (outer − inner ellipsoid), centred at the head origin;
  // the inner cavity is offset to the brain position.
  auto* cOut = new G4Ellipsoid("CraniumOut", kSkullOutAxMM * mm,
                               kSkullOutByMM * mm, kSkullOutCzMM * mm);
  auto* cIn = new G4Ellipsoid("CraniumIn", kSkullInAxMM * mm, kSkullInByMM * mm,
                              kSkullInCzMM * mm);
  auto* skullSolid = new G4SubtractionSolid(
      "Skull", cOut, cIn, nullptr, G4ThreeVector(0, 0, kBrainOffsetZMM * mm));
  auto* skullLV = new G4LogicalVolume(skullSolid, bone, "Skull");
  skullLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.95, 0.95, 0.9, 0.4)));
  new G4PVPlacement(nullptr, {}, skullLV, "Skull", headLV, false, 0, true);

  // Brain ellipsoid, filling the cranium cavity; daughter of the head.
  auto* brainSolid = new G4Ellipsoid("Brain", kBrainAxMM * mm, kBrainByMM * mm,
                                     kBrainCzMM * mm);
  auto* brainLV = new G4LogicalVolume(brainSolid, brainMat, "Brain");
  brainLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.4)));
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kBrainOffsetZMM * mm), brainLV,
                    "Brain", headLV, false, 0, true);

  fPhantomLV = headLV;  // scoring volume = whole head (GetMass includes daughters)
  fBeamHalfExtent = kScalpAxMM * mm;  // L-R semi-axis, now along +z
  // Depth reference = the entrance face. The beam enters the scalp at -fHalfZ, so
  // fHalfZ must be the beam-axis half-extent (not the cylinder length), else the
  // target box and the reported depths are mis-registered against the head.
  fHalfZ = fBeamHalfExtent;

  // Medium regions (world frame), priority-ordered: brain carves the skull
  // shell, skull carves the scalp.
  fRegions = {
      HeadRegion("brain", kBrainMaterial, 0., 0., kBrainOffsetZMM, kBrainAxMM,
                 kBrainByMM, kBrainCzMM),
      HeadRegion("skull", kSkullMaterial, 0., 0., 0., kSkullOutAxMM,
                 kSkullOutByMM, kSkullOutCzMM),
      HeadRegion("scalp", kScalpMaterial, 0., 0., 0., kScalpAxMM, kScalpByMM,
                 kScalpCzMM),
  };
}

// Uniform head (Phase 2): the SAME outer envelope as the MIRD head (the scalp
// ellipsoid) but a single homogeneous material (brain). Same shape as the
// 3-region head, so the two cases isolate the effect of the skull/scalp on the
// proton range and the isotope mix. One medium region.
void DetectorConstruction::BuildUniformHead(G4LogicalVolume* worldLV) {
  using namespace stageA;
  auto* nist = G4NistManager::Instance();
  G4Material* brainMat = nist->FindOrBuildMaterial(kBrainMaterial);

  auto* headSolid = new G4Ellipsoid("Head", kScalpAxMM * mm, kScalpByMM * mm,
                                    kScalpCzMM * mm);
  fPhantomLV = new G4LogicalVolume(headSolid, brainMat, "Head");
  fPhantomLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.3)));

  const G4Transform3D tf =
      G4Translate3D(-kBrainOffsetZMM * mm, 0., 0.) * G4RotateY3D(90. * deg);
  new G4PVPlacement(tf, fPhantomLV, "Head", worldLV, false, 0, true);
  fBeamHalfExtent = kScalpAxMM * mm;
  fHalfZ = fBeamHalfExtent;  // depth reference = entrance face (see BuildMirdHead)

  fRegions = {HeadRegion("head", kBrainMaterial, 0., 0., 0., kScalpAxMM,
                         kScalpByMM, kScalpCzMM)};
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
  // run_meta.csv "phantom_material": the single material when the phantom is one
  // region (cylinder, uniform head), else "multi" (the medium is per-region in
  // phantom_regions.csv). The case itself is the separate "geometry" column.
  return (fRegions.size() == 1) ? fRegions[0].material : G4String("multi");
}

G4double DetectorConstruction::PhantomMass() const {
  // G4LogicalVolume::GetMass() integrates ρ·V over the volume (and daughters).
  return fPhantomLV ? fPhantomLV->GetMass() : 0.;
}

// Material of the first priority-ordered region containing p, else nullptr (air).
// Mirrors the phantom_regions.csv point→material rule (ellipsoid / cylinder,
// axis-aligned). Lengths in PhantomRegion are mm.
const G4Material* DetectorConstruction::MaterialAt(const G4ThreeVector& p) const {
  auto* nist = G4NistManager::Instance();
  for (const auto& r : fRegions) {
    const G4double dx = p.x() - r.cx * mm;
    const G4double dy = p.y() - r.cy * mm;
    const G4double dz = p.z() - r.cz * mm;
    bool inside = false;
    if (r.solid == "cylinder") {  // (a,b,c) = (radius, radius, half-length)
      inside = (dx * dx + dy * dy) <= (r.a * mm) * (r.a * mm) &&
               std::abs(dz) <= r.c * mm;
    } else {  // ellipsoid: (a,b,c) = semi-axes
      const G4double fx = dx / (r.a * mm), fy = dy / (r.b * mm),
                     fz = dz / (r.c * mm);
      inside = (fx * fx + fy * fy + fz * fz) <= 1.0;
    }
    if (inside) return nist->FindOrBuildMaterial(r.material);
  }
  return nullptr;
}

G4double DetectorConstruction::TargetMass() const {
  // The target box is conceptual (not a volume), so compute its mass directly
  // from its dimensions and the density of the medium at the box centre (e.g.
  // brain for the head — not the scalp mother LV the box is nested in).
  if (!fPhantomLV) return 0.;
  const G4double zc = 0.5 * (TargetProxZ() + TargetDistZ());
  const G4Material* mat = MaterialAt(G4ThreeVector(0., 0., zc));
  if (!mat) mat = fPhantomLV->GetMaterial();  // fallback (box centre in air)
  const G4double length = fTargetDistDepth - fTargetProxDepth;
  const G4double volume = pi * fTargetRadius * fTargetRadius * length;
  return volume * mat->GetDensity();
}
