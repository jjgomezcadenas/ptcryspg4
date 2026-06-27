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

// Central-axis depth-dose core (Phase 3 Step 1): a thin on-axis cylinder whose
// fixed cross-section makes edep(z) ∝ dose within constant-density material —
// the clean profile graded for SOBP plateau flatness / R80 through the head
// (the full-plane edep tally conflates the SOBP shape with the ellipsoid's
// varying cross-section and the bone shells).
inline constexpr double kCoreRadiusMM = 5.0;

// --- MIRD-head variant (Phase 1) ------------------------------------------
// A heterogeneous stylized head: brain inside a skull shell inside a soft-tissue
// outer head, replicated from the Geant4 advanced/human_phantom MIRD model.
// Semi-axes [mm], in the head-local frame (x = left-right, y = anterior-
// posterior, z = superior-inferior). Selected by /stageA/phantom/geometry
// mird_head; the cylinder above remains the default.
inline constexpr const char* kGeometryCylinder = "cylinder";
inline constexpr const char* kGeometryUniformHead = "uniform_head";  // 1 region
inline constexpr const char* kGeometryMirdHead = "mird_head";        // 3 regions
inline constexpr const char* kDefaultGeometry = kGeometryCylinder;

// Head-local frame origin = skull/scalp centre; the brain sits +kBrainOffsetZMM
// toward the crown (the MIRD cranium cavity offset). The face/neck of the full
// MIRD head is dropped: the lateral cranial field never traverses it, so the
// outer head is a single soft-tissue scalp ellipsoid enclosing the skull.
// Scalp = skull-outer semi-axes + ~4 mm.
inline constexpr double kScalpAxMM = 72.0;
inline constexpr double kScalpByMM = 102.0;
inline constexpr double kScalpCzMM = 87.0;
// Skull cranium shell = outer ellipsoid minus inner cavity (offset by the brain).
inline constexpr double kSkullOutAxMM = 68.0;
inline constexpr double kSkullOutByMM = 98.0;
inline constexpr double kSkullOutCzMM = 83.0;
inline constexpr double kSkullInAxMM = 60.0;
inline constexpr double kSkullInByMM = 90.0;
inline constexpr double kSkullInCzMM = 65.0;
// Brain ellipsoid (fills the cranium cavity), offset from the skull centre.
inline constexpr double kBrainAxMM = 60.0;
inline constexpr double kBrainByMM = 90.0;
inline constexpr double kBrainCzMM = 65.0;
inline constexpr double kBrainOffsetZMM = 10.0;  // brain centre vs skull centre
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
