#ifndef PTCRYSPG4_COMMON_ISOTOPES_HH
#define PTCRYSPG4_COMMON_ISOTOPES_HH

// Isotope encoding for the ptcryspg4 simulation chain.
//
// MIRRORS common/SCHEMA.md (the authoritative contract) — do not edit
// independently. The Python mirror is common/isotopes.py.
//
// Stage A never needs half-lives -- Geant4's radioactive-decay process handles
// the decay internally; what C++ needs is the mapping: residual nucleus (Z, A)
// -> listed beta+ emitter? -> isotope_id. Half-lives and endpoints live on the
// Python side (handoff bookkeeping).

#include <cstdint>

namespace ptcrysp {

struct EmitterDef {
  int Z;
  int A;
  std::int8_t id;
  const char* name;
};

// isotope_id encoding (SCHEMA.md): 0=15O  1=11C  2=13N  3=10C  4=14O
inline constexpr EmitterDef kEmitters[] = {
    {8, 15, 0, "O15"},
    {6, 11, 1, "C11"},
    {7, 13, 2, "N13"},
    {6, 10, 3, "C10"},
    {8, 14, 4, "O14"},
};

inline constexpr int kNEmitters = sizeof(kEmitters) / sizeof(kEmitters[0]);

// isotope_id for residual nucleus (Z, A); -1 if not a listed emitter.
inline constexpr std::int8_t EmitterId(int Z, int A) {
  for (const auto& e : kEmitters) {
    if (e.Z == Z && e.A == A) return e.id;
  }
  return -1;
}

// Name for an isotope_id; nullptr if invalid.
inline constexpr const char* EmitterName(std::int8_t id) {
  for (const auto& e : kEmitters) {
    if (e.id == id) return e.name;
  }
  return nullptr;
}

}  // namespace ptcrysp

#endif  // PTCRYSPG4_COMMON_ISOTOPES_HH
