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

// Physics
inline constexpr const char* kPhysicsList = "QGSP_BIC_HP";

// Dose scorer names: collection is "<MFD>/<primitive>" = "phantomMFD/edep".
inline constexpr const char* kScorerMFD = "phantomMFD";
inline constexpr const char* kScorerEdep = "edep";

}  // namespace stageA

#endif  // STAGEA_CONFIG_HH
