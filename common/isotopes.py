"""Isotope table for the ptcryspg4 simulation chain.

MIRRORS common/SCHEMA.md (the authoritative contract) -- do not edit
independently. The C++ mirror is common/Isotopes.hh.

The handoff (Stage B0) needs isotope_id -> (name, half-life, endpoint) to
compute the decay constants lambda_j and the measured-window counts N_j
(spec Eq. 7). The analysis uses it for labelling.
"""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Isotope:
    name: str
    Z: int
    A: int
    half_life_s: float    # physical half-life [s]
    endpoint_MeV: float   # beta+ endpoint energy [MeV] (governs positron range)
    prompt_gamma: bool    # emits prompt de-excitation gamma in coincidence

    @property
    def lam(self) -> float:
        """Decay constant lambda = ln2 / T_half  [1/s]."""
        return math.log(2.0) / self.half_life_s


# isotope_id encoding (SCHEMA.md): 0=15O  1=11C  2=13N  3=10C  4=14O
ISOTOPES: dict[int, Isotope] = {
    0: Isotope("O15", Z=8, A=15, half_life_s=122.24, endpoint_MeV=1.74, prompt_gamma=False),
    1: Isotope("C11", Z=6, A=11, half_life_s=1223.4, endpoint_MeV=0.96, prompt_gamma=False),
    2: Isotope("N13", Z=7, A=13, half_life_s=597.9, endpoint_MeV=1.19, prompt_gamma=False),
    3: Isotope("C10", Z=6, A=10, half_life_s=19.29, endpoint_MeV=1.91, prompt_gamma=True),
    4: Isotope("O14", Z=8, A=14, half_life_s=70.62, endpoint_MeV=1.81, prompt_gamma=True),
}

N_ISOTOPES = len(ISOTOPES)

NAME_TO_ID = {iso.name: iid for iid, iso in ISOTOPES.items()}
