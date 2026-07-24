#include "ThinTargetDetector.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

ThinTargetDetector::ThinTargetDetector(const ThinTargetConfig& config)
    : config_(config) {}

G4VPhysicalVolume* ThinTargetDetector::Construct() {
  G4String symbol;
  G4double atomic_mass = 0.0;
  if (config_.target == "C12") {
    symbol = "C"; target_z_ = 6; target_a_ = 12; atomic_mass = 12.0; density_g_cm3_ = 2.00;
  } else if (config_.target == "N14") {
    symbol = "N"; target_z_ = 7; target_a_ = 14; atomic_mass = 14.003074; density_g_cm3_ = 0.81;
  } else if (config_.target == "O16") {
    symbol = "O"; target_z_ = 8; target_a_ = 16; atomic_mass = 15.994915; density_g_cm3_ = 1.14;
  } else {
    G4Exception("ThinTargetDetector::Construct", "InvalidTarget", FatalException,
                "Target must be C12, N14 or O16");
  }

  target_isotope_ = new G4Isotope(
      config_.target, target_z_, target_a_, atomic_mass * g / mole);
  target_element_ = new G4Element(
      "isotopically pure " + config_.target, symbol, 1);
  target_element_->AddIsotope(target_isotope_, 100.0 * perCent);
  target_material_ = new G4Material(
      config_.target + "_target", density_g_cm3_ * g / cm3, 1);
  target_material_->AddElement(target_element_, 1.0);

  thickness_mm_ = (config_.areal_mg_cm2 * 1.0e-3 / density_g_cm3_) * 10.0;
  auto* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto* world_solid = new G4Box("world", 10.0 * cm, 10.0 * cm, 10.0 * cm);
  auto* world_logical = new G4LogicalVolume(world_solid, vacuum, "world");
  auto* world = new G4PVPlacement(nullptr, {}, world_logical, "world", nullptr, false, 0);
  auto* target_solid = new G4Box("target", 5.0 * cm, 5.0 * cm,
                                 0.5 * thickness_mm_ * mm);
  target_logical_ = new G4LogicalVolume(target_solid, target_material_, "target");
  new G4PVPlacement(nullptr, {}, target_logical_, "target", world_logical, false, 0);
  return world;
}
