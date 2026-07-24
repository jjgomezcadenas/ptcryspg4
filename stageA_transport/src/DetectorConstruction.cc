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
// Head orientation relative to the fixed +z beam: which head-local axis lies
// along the beam. Lateral (L-R along +z) is the mird_head / uniform_head base
// case; Posterior (A-P along +z, beam through the occiput) is the headep case.
enum class HeadAxis { Lateral, Posterior };

// Map a head-local ellipsoid (centre c, semi-axes s) to a world-frame
// PhantomRegion, applying the same 90° rotation + translation T [mm] that places
// the head. A 90° rotation keeps an axis-aligned ellipsoid axis-aligned, so only
// the axes and centre are permuted/signed:
//   Lateral   world = ( cz, cy, -cx) + T ,  semi = (sz, sy, sx)
//   Posterior world = ( cx, -cz, cy) + T ,  semi = (sx, sz, sy)
PhantomRegion HeadRegionOriented(HeadAxis axis, G4double Tx, G4double Ty,
                                 G4double Tz, const char* name, const char* mat,
                                 G4double cx, G4double cy, G4double cz,
                                 G4double sx, G4double sy, G4double sz) {
  auto nz = [](G4double v) { return v == 0. ? 0. : v; };  // avoid -0.0 in the CSV
  if (axis == HeadAxis::Lateral)
    return PhantomRegion{name, mat, "ellipsoid", sz, sy, sx,
                         nz(cz + Tx), nz(cy + Ty), nz(-cx + Tz), 0., 0., 0.};
  return PhantomRegion{name, mat, "ellipsoid", sx, sz, sy,
                       nz(cx + Tx), nz(-cz + Ty), nz(cy + Tz), 0., 0., 0.};
}

// Lateral head region (the base case): translation T = (-kBrainOffsetZ, 0, 0),
// which centres the brain on the beam axis (world x = y = 0).
PhantomRegion HeadRegion(const char* name, const char* mat, G4double cx,
                         G4double cy, G4double cz, G4double sx, G4double sy,
                         G4double sz) {
  return HeadRegionOriented(HeadAxis::Lateral, -stageA::kBrainOffsetZMM, 0., 0.,
                            name, mat, cx, cy, cz, sx, sy, sz);
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
                     fGeometry == stageA::kGeometryUniformHead ||
                     fGeometry == stageA::kGeometryHeadEP ||
                     fGeometry == stageA::kGeometryUniformHeadEP);

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

  fScoringLVs.clear();
  if (fGeometry == stageA::kGeometryMirdHead)
    BuildMirdHead(worldLV);
  else if (fGeometry == stageA::kGeometryHeadEP)
    BuildHeadEP(worldLV);
  else if (fGeometry == stageA::kGeometryUniformHeadEP)
    BuildUniformHeadEP(worldLV);
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
  fScoringLVs.push_back(fPhantomLV);
  fBeamHalfExtent = fHalfZ;

  // One homogeneous region (world frame): a cylinder centred at the origin.
  fRegions = {{"phantom", fMaterialName, "cylinder", fRadius / mm, fRadius / mm,
               fHalfZ / mm, 0., 0., 0., 0., 0., 0.}};
}

// Build the scalp/skull/brain envelope into worldLV with the given placement
// transform (which sets the head's orientation relative to the +z beam). The
// head is a soft-tissue scalp ellipsoid ⊃ bone skull shell ⊃ brain (skull/brain
// are the MIRD ellipsoids; the MIRD face/neck is dropped). Built in the
// head-local frame (x = L-R, y = A-P, z = S-I), origin at the skull/scalp centre.
// Returns the head (scoring) and brain (for daughters) logical volumes.
void DetectorConstruction::BuildHeadEnvelope(G4LogicalVolume* worldLV,
                                             const G4Transform3D& tf,
                                             G4LogicalVolume** headLVout,
                                             G4LogicalVolume** brainLVout) {
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

  fScoringLVs.insert(fScoringLVs.end(), {headLV, skullLV, brainLV});
  *headLVout = headLV;
  *brainLVout = brainLV;
}

// Heterogeneous head: the scalp/skull/brain envelope placed laterally,
// so the L-R axis lies along the beam (+z) and the brain centre is at the world
// origin. The base case, unchanged.
void DetectorConstruction::BuildMirdHead(G4LogicalVolume* worldLV) {
  using namespace stageA;
  const G4Transform3D tf =
      G4Translate3D(-kBrainOffsetZMM * mm, 0., 0.) * G4RotateY3D(90. * deg);
  G4LogicalVolume *headLV, *brainLV;
  BuildHeadEnvelope(worldLV, tf, &headLV, &brainLV);

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

// HeadEP: the MIRD head plus a posterior-fossa tumour, oriented so the beam runs
// straight through the occiput — the head A-P axis lies along +z (RotateX 90°),
// beam entering the back of the head travelling anterior. The translation T
// centres the tumour on the beam axis: under world = (x, -z, y), the tumour's
// transverse coords are (x, -z), so they are cancelled; Tz = 0 keeps the head
// centred along the beam (entrance at -A-P semi-axis).
void DetectorConstruction::BuildHeadEP(G4LogicalVolume* worldLV) {
  using namespace stageA;
  auto* nist = G4NistManager::Instance();

  const G4double Tx = -kTumourPosXMM;
  const G4double Ty = kTumourPosZMM;  // = -(-kTumourPosZMM), cancels world y
  const G4double Tz = 0.;
  const G4Transform3D tf =
      G4Translate3D(Tx * mm, Ty * mm, Tz * mm) * G4RotateX3D(90. * deg);

  G4LogicalVolume *headLV, *brainLV;
  BuildHeadEnvelope(worldLV, tf, &headLV, &brainLV);

  // Tumour ellipsoid, fully inside the brain (daughter of the brain LV). The
  // brain LV frame is the head-local frame shifted by -kBrainOffsetZ in z.
  G4Material* tumourMat = nist->FindOrBuildMaterial(kTumourMaterial);
  auto* tumSolid = new G4Ellipsoid("Tumour", kTumourAxMM * mm, kTumourByMM * mm,
                                   kTumourCzMM * mm);
  auto* tumLV = new G4LogicalVolume(tumSolid, tumourMat, "Tumour");
  tumLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.9, 0.2, 0.2, 0.6)));
  new G4PVPlacement(
      nullptr,
      G4ThreeVector(kTumourPosXMM * mm, kTumourPosYMM * mm,
                    (kTumourPosZMM - kBrainOffsetZMM) * mm),
      tumLV, "Tumour", brainLV, false, 0, true);
  fScoringLVs.push_back(tumLV);

  fPhantomLV = headLV;
  fBeamHalfExtent = kScalpByMM * mm;  // A-P semi-axis, now along +z
  fHalfZ = fBeamHalfExtent;

  // Medium regions (world frame), priority-ordered: tumour carves the brain,
  // brain the skull shell, skull the scalp.
  fRegions = {
      HeadRegionOriented(HeadAxis::Posterior, Tx, Ty, Tz, "tumour",
                         kTumourMaterial, kTumourPosXMM, kTumourPosYMM,
                         kTumourPosZMM, kTumourAxMM, kTumourByMM, kTumourCzMM),
      HeadRegionOriented(HeadAxis::Posterior, Tx, Ty, Tz, "brain", kBrainMaterial,
                         0., 0., kBrainOffsetZMM, kBrainAxMM, kBrainByMM,
                         kBrainCzMM),
      HeadRegionOriented(HeadAxis::Posterior, Tx, Ty, Tz, "skull", kSkullMaterial,
                         0., 0., 0., kSkullOutAxMM, kSkullOutByMM, kSkullOutCzMM),
      HeadRegionOriented(HeadAxis::Posterior, Tx, Ty, Tz, "scalp", kScalpMaterial,
                         0., 0., 0., kScalpAxMM, kScalpByMM, kScalpCzMM),
  };
}

// Uniform HeadEP: the SAME envelope, posterior placement, and target as headep
// — the scalp ellipsoid with its A-P axis along the beam, translated so the
// headep tumour site sits on the beam axis — but brain throughout, no layers,
// no tumour region. The homogeneity control for the σ_R study: comparing
// headep against this isolates what the tissue layering does to the distal
// range, with the orientation, target window, and outer attenuation shell
// held identical. One medium region.
void DetectorConstruction::BuildUniformHeadEP(G4LogicalVolume* worldLV) {
  using namespace stageA;
  auto* nist = G4NistManager::Instance();
  G4Material* brainMat = nist->FindOrBuildMaterial(kBrainMaterial);

  auto* headSolid = new G4Ellipsoid("Head", kScalpAxMM * mm, kScalpByMM * mm,
                                    kScalpCzMM * mm);
  fPhantomLV = new G4LogicalVolume(headSolid, brainMat, "Head");
  fPhantomLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.9, 0.3)));

  // The headep placement, verbatim (see BuildHeadEP for the translation logic).
  const G4double Tx = -kTumourPosXMM;
  const G4double Ty = kTumourPosZMM;
  const G4Transform3D tf =
      G4Translate3D(Tx * mm, Ty * mm, 0.) * G4RotateX3D(90. * deg);
  new G4PVPlacement(tf, fPhantomLV, "Head", worldLV, false, 0, true);
  fScoringLVs.push_back(fPhantomLV);
  fBeamHalfExtent = kScalpByMM * mm;  // A-P semi-axis, now along +z
  fHalfZ = fBeamHalfExtent;

  fRegions = {HeadRegionOriented(HeadAxis::Posterior, Tx, Ty, 0., "head",
                                 kBrainMaterial, 0., 0., 0., kScalpAxMM,
                                 kScalpByMM, kScalpCzMM)};
}

// Uniform head: the SAME outer envelope as the MIRD head (the scalp
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
  fScoringLVs.push_back(fPhantomLV);
  fBeamHalfExtent = kScalpAxMM * mm;
  fHalfZ = fBeamHalfExtent;  // depth reference = entrance face (see BuildMirdHead)

  fRegions = {HeadRegion("head", kBrainMaterial, 0., 0., 0., kScalpAxMM,
                         kScalpByMM, kScalpCzMM)};
}

void DetectorConstruction::ConstructSDandField() {
  // A multifunctional detector is a container of "primitive scorers"; here a
  // single G4PSEnergyDeposit sums the energy deposited in the phantom each
  // event. G4 merges its hits map across threads automatically. Binding to a
  // volume covers that volume only, so every phantom volume (mother and
  // daughters) is bound.
  auto* mfd = new G4MultiFunctionalDetector(stageA::kScorerMFD);
  G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
  mfd->RegisterPrimitive(new G4PSEnergyDeposit(stageA::kScorerEdep));
  for (auto* lv : fScoringLVs) SetSensitiveDetector(lv, mfd);
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
