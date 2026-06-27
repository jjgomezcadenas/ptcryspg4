#ifndef STAGEA_DETECTORCONSTRUCTION_HH
#define STAGEA_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

// One medium region for phantom_regions.csv (Phase 2): a world-frame solid with
// its NIST material. Ellipsoid uses (a,b,c) = semi-axes; cylinder uses
// (a,b,c) = (radius, radius, half-length). Euler angles are intrinsic X-Y-Z in
// degrees and are 0 here (the regions are axis-aligned in the world frame).
struct PhantomRegion {
  G4String name, material, solid;
  G4double a, b, c;        // mm
  G4double cx, cy, cz;     // mm, world centre
  G4double ex, ey, ez;     // deg, intrinsic XYZ Euler (0 = axis-aligned)
};

// Stage A geometry (spec sec 2.2): a homogeneous right-circular cylinder, axis
// along +z, that serves both as the β⁺ production target and as the passive
// 511 keV attenuator. Geant4 calls the two virtuals below automatically: it is
// the mandatory "what does the apparatus look like" user class.
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction();             // builds the messenger; sets default material
  ~DetectorConstruction() override;   // deletes the messenger

  // Build the volume tree (world + phantom) and return the world physical
  // volume. Called once by the run manager at /run/initialize.
  G4VPhysicalVolume* Construct() override;

  // Attach the energy-deposit scorer to the phantom. Called after Construct(),
  // once per worker thread in MT, so each thread owns its sensitive detector.
  void ConstructSDandField() override;

  // Phantom material (NIST name, e.g. G4_TISSUE_SOFT_ICRP). Settable from a
  // macro before /run/initialize via DetectorMessenger; default from StageAConfig.
  void SetMaterialName(const G4String& name) { fMaterialName = name; }
  const G4String& MaterialName() const { return fMaterialName; }

  // Geometry selector: "cylinder" (default) or "mird_head". Set from a macro
  // before /run/initialize. The label written to run_meta.csv is the NIST
  // material for the cylinder, or "MIRD_head" for the head.
  void SetGeometry(const G4String& name) { fGeometry = name; }
  const G4String& Geometry() const { return fGeometry; }
  G4String PhantomLabel() const;

  // The medium regions (world frame), for phantom_regions.csv. Priority-ordered:
  // the first region containing a point owns it. Valid after Construct().
  const std::vector<PhantomRegion>& Regions() const { return fRegions; }

  // Beam placement helpers, read by the primary generator.
  G4double PhantomHalfLength() const { return fHalfZ; }
  G4double PhantomRadius() const { return fRadius; }
  // Half-extent of the phantom along the beam axis (+z): the gun starts just
  // upstream of this. Cylinder → half-length; MIRD head → its L-R semi-axis.
  G4double BeamAxisHalfExtent() const { return fBeamHalfExtent; }
  // Phantom mass (ρ·V), for the dose normalization. Valid only after Construct().
  G4double PhantomMass() const;

  // Target (tumour) box for dose normalization. It is NOT a geometry volume — it
  // is a scoring region the SteppingAction tests each step against (transparent,
  // no carving). Depths are from the entrance face; settable by card
  // (/stageA/target/...), defaults from StageAConfig.
  void SetTargetRadius(G4double r) { fTargetRadius = r; }
  void SetTargetProxDepth(G4double d) { fTargetProxDepth = d; }
  void SetTargetDistDepth(G4double d) { fTargetDistDepth = d; }
  G4double TargetRadius() const { return fTargetRadius; }
  G4double TargetProxZ() const { return -fHalfZ + fTargetProxDepth; }  // beam at -fHalfZ
  G4double TargetDistZ() const { return -fHalfZ + fTargetDistDepth; }
  G4double TargetMass() const;  // π·r²·L·ρ (ρ from the medium at the box centre)
  // Medium at the target-box centre, or nullptr if it is in air — a mis-placed
  // box (depths beyond the phantom) makes the dose normalization meaningless.
  const G4Material* TargetMaterial() const {
    return MaterialAt(G4ThreeVector(0., 0., 0.5 * (TargetProxZ() + TargetDistZ())));
  }

  // Material of the first (priority-ordered) region containing a world point,
  // else nullptr (air). Same rule as phantom_regions.csv; used for the target-box
  // density and the per-bin density of the central-axis dose profile.
  const G4Material* MaterialAt(const G4ThreeVector& p) const;

 private:
  // Build the homogeneous cylinder (default) or the MIRD head into worldLV; each
  // sets fPhantomLV (the scoring volume) and fBeamHalfExtent.
  void BuildCylinder(G4LogicalVolume* worldLV);
  void BuildUniformHead(G4LogicalVolume* worldLV);  // 1-region brain ellipsoid
  void BuildMirdHead(G4LogicalVolume* worldLV);     // 3-region scalp/skull/brain

  G4double fRadius = 0.;   // cylinder radius, set in Construct()
  G4double fHalfZ = 0.;    // cylinder half-length (G4Tubs uses half-z)
  G4double fBeamHalfExtent = 0.;          // phantom half-extent along +z
  G4LogicalVolume* fPhantomLV = nullptr;  // scoring volume (cylinder or head)
  G4String fMaterialName;                 // NIST name; default set in the ctor
  G4String fGeometry;                     // "cylinder" | "mird_head"
  G4double fTargetRadius = 0.;            // target box, defaults set in the ctor
  G4double fTargetProxDepth = 0.;
  G4double fTargetDistDepth = 0.;
  std::vector<PhantomRegion> fRegions;   // medium regions, filled in the builders
  DetectorMessenger* fMessenger = nullptr;
};

#endif  // STAGEA_DETECTORCONSTRUCTION_HH
