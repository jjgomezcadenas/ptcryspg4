#!/usr/bin/env python3
"""Phantom material data for the downstream PET simulation and reconstruction.

The source this repo produces is detector-independent but **not medium-
independent** (SCHEMA.md, "Medium / attenuation map"): the downstream sim must
propagate the 511 keV annihilation photons through the phantom and needs the
linear attenuation coefficient mu(511 keV) for reconstruction attenuation
correction. For the homogeneous standard scenario the medium is fully specified
by the phantom material name + cylinder geometry recorded in run_meta.csv; this
script turns the material name into the numbers a consumer actually needs:
density, elemental composition, and mu/mu-over-rho at the annihilation energy.

The compositions are the authoritative Geant4 NIST definitions (verbatim from a
Geant4 11.4.1 material dump). G4_BRAIN_ICRP was read from the B3 example output
(`examples/basic/B3/.../exampleB3.out`); re-dump with `G4NistManager::Instance()
->PrintG4Material("G4_BRAIN_ICRP")` to verify against another build.

mu/rho is computed from first principles: at 511 keV in low-Z biological matter
the photon interaction is Compton-dominated, so mu/rho = sum_i w_i (N_A Z_i/A_i)
sigma_KN(E), with sigma_KN the Klein-Nishina total cross section per electron.
Coherent (Rayleigh) + photoelectric add ~1-2% at 511 keV and are not included
here; the Compton value lands within ~0.5% of the NIST XCOM total for water
(0.0958 cm^2/g), which is adequate for an attenuation-correction map. Use a full
XCOM/NIST table if sub-percent accuracy is required.

Usage:
    python common/phantom_material.py [--material NAME] [--energy-keV E]
        [--meta run_meta.csv] [--csv [DIR]]

    # default: print the G4_BRAIN_ICRP card at 511 keV
    python common/phantom_material.py
    # pick the material from a frozen run and write CSVs into data/
    python common/phantom_material.py --meta data/run_meta.csv --csv data
"""

import argparse
import math
import os
import re
import sys
from dataclasses import dataclass

# Avogadro's number [1/mol] and classical electron radius [cm].
_N_A = 6.02214076e23
_R_E = 2.8179403262e-13
# Electron rest energy [keV] (sets the Klein-Nishina energy scale).
_MEC2_keV = 510.99895


@dataclass(frozen=True)
class Element:
    Z: int
    A: float          # standard atomic weight [g/mol]
    mass_fraction: float


@dataclass(frozen=True)
class Material:
    name: str                 # Geant4 NIST name
    density_g_cm3: float
    mean_excitation_eV: float  # Geant4 Imean (diagnostic; not used for mu)
    composition: dict          # element symbol -> Element
    note: str = ""


# Authoritative Geant4 11.4.1 NIST definitions. Mass fractions sum to 1.
# Add a material by dumping it from Geant4 (see module docstring) -- do not guess.
MATERIALS: dict[str, Material] = {
    # Standard scenario phantom -- read verbatim from the 11.4.1 B3 dump.
    "G4_BRAIN_ICRP": Material(
        name="G4_BRAIN_ICRP", density_g_cm3=1.040, mean_excitation_eV=73.3,
        composition={
            "H":  Element(1,  1.008,  0.107),
            "C":  Element(6,  12.011, 0.145),
            "N":  Element(7,  14.007, 0.022),
            "O":  Element(8,  15.999, 0.712),
            "Na": Element(11, 22.990, 0.002),
            "P":  Element(15, 30.974, 0.004),
            "S":  Element(16, 32.060, 0.002),
            "Cl": Element(17, 35.450, 0.003),
            "K":  Element(19, 39.098, 0.003),
        },
        note="ICRP brain; standard ptcryspg4 phantom.",
    ),
    # Pure-15O extreme reference (Geant4 default Imean = 78 eV).
    "G4_WATER": Material(
        name="G4_WATER", density_g_cm3=1.000, mean_excitation_eV=78.0,
        composition={
            "H": Element(1, 1.008,  0.111894),
            "O": Element(8, 15.999, 0.888106),
        },
        note="Reference medium; mu/rho here matches the canonical PET 0.096 cm^2/g.",
    ),
    # Carbon-rich PMMA benchmark, stoichiometry C5H8O2.
    "G4_PLEXIGLASS": Material(
        name="G4_PLEXIGLASS", density_g_cm3=1.190, mean_excitation_eV=74.0,
        composition={
            "H": Element(1, 1.008,  0.080538),
            "C": Element(6, 12.011, 0.599848),
            "O": Element(8, 15.999, 0.319614),
        },
        note="PMMA; carbon-rich (11C-dominated) cross-check material.",
    ),
}

DEFAULT_MATERIAL = "G4_BRAIN_ICRP"


def klein_nishina_cm2(energy_keV: float) -> float:
    """Klein-Nishina total cross section per electron [cm^2] at photon energy E."""
    a = energy_keV / _MEC2_keV
    term1 = (1 + a) / a**2 * (2 * (1 + a) / (1 + 2 * a) - math.log(1 + 2 * a) / a)
    term2 = math.log(1 + 2 * a) / (2 * a)
    term3 = (1 + 3 * a) / (1 + 2 * a) ** 2
    return 2 * math.pi * _R_E**2 * (term1 + term2 - term3)


def mu_over_rho_cm2_g(material: Material, energy_keV: float) -> float:
    """Compton-dominated mass attenuation coefficient mu/rho [cm^2/g]."""
    sigma = klein_nishina_cm2(energy_keV)
    return sum(e.mass_fraction * sigma * _N_A * (e.Z / e.A)
              for e in material.composition.values())


def attenuation(material: Material, energy_keV: float) -> dict:
    """mu/rho, mu, and mean free path for a material at a photon energy."""
    mu_rho = mu_over_rho_cm2_g(material, energy_keV)
    mu_cm = mu_rho * material.density_g_cm3
    return {
        "material": material.name,
        "energy_keV": energy_keV,
        "density_g_cm3": material.density_g_cm3,
        "mean_excitation_eV": material.mean_excitation_eV,
        "mu_rho_cm2_g": mu_rho,
        "mu_cm_inv": mu_cm,
        "mu_mm_inv": mu_cm / 10.0,
        "mean_free_path_cm": float("inf") if mu_cm == 0 else 1.0 / mu_cm,
    }


def _material_from_meta(meta_path: str) -> str:
    """Read the phantom_material column from a run_meta.csv."""
    import pandas as pd
    return str(pd.read_csv(meta_path).iloc[0]["phantom_material"])


def _print_card(material: Material, energy_keV: float) -> None:
    frac_sum = sum(e.mass_fraction for e in material.composition.values())
    att = attenuation(material, energy_keV)
    print(f"\n{material.name}   (rho = {material.density_g_cm3:g} g/cm^3, "
          f"Imean = {material.mean_excitation_eV:g} eV)")
    if material.note:
        print(f"  {material.note}")
    print(f"\n  {'elem':>4} {'Z':>3} {'A [g/mol]':>10} {'mass frac':>10}")
    print("  " + "-" * 31)
    for sym, e in material.composition.items():
        print(f"  {sym:>4} {e.Z:>3} {e.A:>10.3f} {e.mass_fraction:>10.6f}")
    print(f"  {'sum':>4} {'':>3} {'':>10} {frac_sum:>10.6f}")
    print(f"\n  photon attenuation @ {energy_keV:g} keV (Compton/Klein-Nishina):")
    print(f"    mu/rho           = {att['mu_rho_cm2_g']:.5f} cm^2/g")
    print(f"    mu               = {att['mu_cm_inv']:.5f} cm^-1  "
          f"(= {att['mu_mm_inv']:.6f} mm^-1)")
    print(f"    mean free path   = {att['mean_free_path_cm']:.2f} cm")
    print("  (+ ~1-2% from coherent + photoelectric, not included)")


def material_tag(name: str) -> str:
    """Filename-safe tag for a material name, e.g. G4_BRAIN_ICRP -> g4_brain_icrp."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def write_material_csv(material: Material, energy_keV: float, out_dir: str) -> tuple:
    """Write composition + meta CSV for a material; return the two paths."""
    import pandas as pd
    tag = material_tag(material.name)
    comp = pd.DataFrame(
        [{"element": s, "Z": e.Z, "A_g_mol": e.A, "mass_fraction": e.mass_fraction}
         for s, e in material.composition.items()])
    cpath = os.path.join(out_dir, f"phantom_material_{tag}.csv")
    comp.to_csv(cpath, index=False)

    meta = attenuation(material, energy_keV)
    meta["note"] = material.note
    mpath = os.path.join(out_dir, f"phantom_material_{tag}_meta.csv")
    pd.DataFrame([meta]).to_csv(mpath, index=False, float_format="%.6e")
    return cpath, mpath


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--material", default=None,
                    help=f"Geant4 NIST name (default {DEFAULT_MATERIAL}); "
                         f"known: {', '.join(MATERIALS)}")
    ap.add_argument("--energy-keV", type=float, default=511.0,
                    help="photon energy [keV] (default 511, annihilation)")
    ap.add_argument("--meta", default=None,
                    help="run_meta.csv to read phantom_material from")
    ap.add_argument("--csv", nargs="?", const=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "data"),
        default=None, metavar="DIR",
        help="also write composition + meta CSV (default dir: data/)")
    args = ap.parse_args()

    name = args.material or (
        _material_from_meta(args.meta) if args.meta else DEFAULT_MATERIAL)
    if name not in MATERIALS:
        sys.exit(f"error: '{name}' not in the registry ({', '.join(MATERIALS)}); "
                 f"add it by dumping the Geant4 NIST definition (see docstring).")
    material = MATERIALS[name]

    _print_card(material, args.energy_keV)
    if args.csv is not None:
        cpath, mpath = write_material_csv(material, args.energy_keV, args.csv)
        print(f"\nwrote composition -> {cpath}")
        print(f"      meta        -> {mpath}")


if __name__ == "__main__":
    main()
