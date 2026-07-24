#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "ThinTargetConfig.hh"

class G4LogicalVolume;
class G4Element;
class G4Isotope;
class G4Material;
class G4VPhysicalVolume;

class ThinTargetDetector : public G4VUserDetectorConstruction {
 public:
  explicit ThinTargetDetector(const ThinTargetConfig& config);
  G4VPhysicalVolume* Construct() override;
  G4LogicalVolume* TargetLogical() const { return target_logical_; }
  G4Material* TargetMaterial() const { return target_material_; }
  G4Element* TargetElement() const { return target_element_; }
  G4Isotope* TargetIsotope() const { return target_isotope_; }
  G4double DensityGcm3() const { return density_g_cm3_; }
  G4double ThicknessMm() const { return thickness_mm_; }
  G4int TargetZ() const { return target_z_; }
  G4int TargetA() const { return target_a_; }

 private:
  ThinTargetConfig config_;
  G4LogicalVolume* target_logical_ = nullptr;
  G4Material* target_material_ = nullptr;
  G4Element* target_element_ = nullptr;
  G4Isotope* target_isotope_ = nullptr;
  G4double density_g_cm3_ = 0.0;
  G4double thickness_mm_ = 0.0;
  G4int target_z_ = 0;
  G4int target_a_ = 0;
};
