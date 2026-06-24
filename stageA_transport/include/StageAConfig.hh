#ifndef STAGEA_CONFIG_HH
#define STAGEA_CONFIG_HH

// Single source of truth for Stage-A configuration constants. Read by the
// geometry, the gun, and the run-metadata writer (so run_meta.csv always
// matches what was simulated). Bare numbers carry the unit named in the suffix;
// multiply by the CLHEP unit (mm, MeV, ...) at the point of use.

namespace stageA {

// Beam (spec sec 2.1)
inline constexpr double kBeamEnergyMeV = 100.0;
inline constexpr double kBeamSigmaMM = 3.0;

// Phantom — Parodi standard: the head approximated as a Ø16 cm × 16 cm cylinder.
inline constexpr double kPhantomDiameterMM = 160.0;
inline constexpr double kPhantomLengthMM = 160.0;
// Brain (ICRP) is the Parodi skull-base reference material. Macro-selectable
// variants: "G4_TISSUE_SOFT_ICRP", "G4_PLEXIGLASS" (PMMA), "G4_WATER".
inline constexpr const char* kPhantomMaterial = "G4_BRAIN_ICRP";

// Target (tumour) box — Parodi standard. Used as the SOBP lateral-fluence disk
// radius and (next step) the nested dose-scoring volume. Depths are measured
// from the beam-entrance face (z = -kPhantomLengthMM/2).
inline constexpr double kTargetRadiusMM = 30.0;      // Ø6 cm
inline constexpr double kTargetProxDepthMM = 55.0;   // proximal depth
inline constexpr double kTargetDistDepthMM = 105.0;  // distal depth (5 cm thick)

// --- MIRD-head variant (Phase 1) ------------------------------------------
// A heterogeneous stylized head: brain inside a skull shell inside a soft-tissue
// outer head, replicated from the Geant4 advanced/human_phantom MIRD model.
// Semi-axes [mm], in the head-local frame (x = left-right, y = anterior-
// posterior, z = superior-inferior). Selected by /stageA/phantom/geometry
// mird_head; the cylinder above remains the default.
inline constexpr const char* kGeometryCylinder = "cylinder";
inline constexpr const char* kGeometryMirdHead = "mird_head";
inline constexpr const char* kDefaultGeometry = kGeometryCylinder;

// Outer head: elliptical tube (face/vault) unioned with an ellipsoid cap.
inline constexpr double kHeadTubeAxMM = 70.0;   // x semi (L-R)
inline constexpr double kHeadTubeByMM = 100.0;  // y semi (A-P)
inline constexpr double kHeadTubeDzMM = 77.5;   // z half-height (tube)
inline constexpr double kHeadCapAxMM = 70.0;
inline constexpr double kHeadCapByMM = 100.0;
inline constexpr double kHeadCapCzMM = 85.0;
inline constexpr double kHeadCapDzMM = 77.5;    // cap union offset (+z)
// Skull cranium shell = outer ellipsoid minus inner ellipsoid (offset +10 mm z).
inline constexpr double kSkullOutAxMM = 68.0;
inline constexpr double kSkullOutByMM = 98.0;
inline constexpr double kSkullOutCzMM = 83.0;
inline constexpr double kSkullInAxMM = 60.0;
inline constexpr double kSkullInByMM = 90.0;
inline constexpr double kSkullInCzMM = 65.0;
inline constexpr double kSkullPosZMM = 77.5;    // skull placement in head frame
// Brain ellipsoid (fills the cranium cavity).
inline constexpr double kBrainAxMM = 60.0;
inline constexpr double kBrainByMM = 90.0;
inline constexpr double kBrainCzMM = 65.0;
inline constexpr double kBrainPosZMM = 87.5;    // brain placement in head frame
// Materials (NIST). Brain matches the cylinder reference.
inline constexpr const char* kScalpMaterial = "G4_TISSUE_SOFT_ICRP";
inline constexpr const char* kSkullMaterial = "G4_BONE_CORTICAL_ICRP";
inline constexpr const char* kBrainMaterial = "G4_BRAIN_ICRP";

// Physics
inline constexpr const char* kPhysicsList = "QGSP_BIC_HP";

// Dose scorer names: collection is "<MFD>/<primitive>" = "phantomMFD/edep".
inline constexpr const char* kScorerMFD = "phantomMFD";
inline constexpr const char* kScorerEdep = "edep";

}  // namespace stageA

#endif  // STAGEA_CONFIG_HH
