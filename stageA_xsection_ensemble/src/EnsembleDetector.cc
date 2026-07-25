#include "EnsembleDetector.hh"

#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"

EnsembleDetector::EnsembleDetector(const ensemble::EnsembleCli& cli)
    : fCli(cli),
      fTargetRadius(cli.target_radius_mm * mm),
      fTargetProx(cli.target_prox_mm * mm),
      fTargetDist(cli.target_dist_mm * mm) {}

G4VPhysicalVolume* EnsembleDetector::Construct() {
  using namespace ensemble;
  auto* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  auto* worldSolid = new G4Box("World", 30. * cm, 30. * cm, 30. * cm);
  auto* worldLV = new G4LogicalVolume(worldSolid, air, "World");
  auto* worldPV =
      new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

  if (fCli.geometry == kGeometryCylinder) {
    G4Material* mat = nist->FindOrBuildMaterial(fCli.material);
    const G4double radius = 0.5 * fCli.phantom_diameter_mm * mm;
    fHalfZ = 0.5 * fCli.phantom_length_mm * mm;
    auto* solid = new G4Tubs("Phantom", 0., radius, fHalfZ, 0., 360. * deg);
    fPhantomLV = new G4LogicalVolume(solid, mat, "Phantom");
    new G4PVPlacement(nullptr, {}, fPhantomLV, "Phantom", worldLV, false, 0,
                      true);
  } else {
    G4Material* brain = nist->FindOrBuildMaterial(kBrainMaterial);
    auto* solid = new G4Ellipsoid("Head", kScalpAxMM * mm, kScalpByMM * mm,
                                  kScalpCzMM * mm);
    fPhantomLV = new G4LogicalVolume(solid, brain, "Head");
    G4Transform3D tf;
    if (fCli.geometry == kGeometryUniformHeadEP) {
      // Posterior placement: A-P axis along the beam, tumour site on the axis
      // (stageA_transport BuildUniformHeadEP, verbatim constants).
      tf = G4Translate3D(-kTumourPosXMM * mm, kTumourPosZMM * mm, 0.) *
           G4RotateX3D(90. * deg);
      fHalfZ = kScalpByMM * mm;
    } else {
      // Lateral placement: L-R axis along the beam (uniform_head).
      tf = G4Translate3D(-kBrainOffsetZMM * mm, 0., 0.) *
           G4RotateY3D(90. * deg);
      fHalfZ = kScalpAxMM * mm;
    }
    new G4PVPlacement(tf, fPhantomLV, "Head", worldLV, false, 0, true);
  }
  return worldPV;
}

G4double EnsembleDetector::PhantomMass() const {
  return fPhantomLV ? fPhantomLV->GetMass() : 0.;
}

G4double EnsembleDetector::TargetMass() const {
  if (!fPhantomLV) return 0.;
  const G4double length = fTargetDist - fTargetProx;
  const G4double volume = pi * fTargetRadius * fTargetRadius * length;
  return volume * fPhantomLV->GetMaterial()->GetDensity();
}
