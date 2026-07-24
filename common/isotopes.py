"""The five positron emitters a proton beam creates in tissue.

One table, keyed by ``isotope_id`` — the integer code that labels every row of
``emitters.csv``: name, nuclear charge and mass, half-life, beta+ endpoint
energy, and whether the decay comes with a prompt de-excitation gamma. The
decay budget (``decay_sampling/budget.py``) takes the half-lives to compute
the decay constants lambda_j and from them the measured decays N_j; the
analysis scripts take the names and endpoints for labelling and ordering.

The authoritative contract is ``common/SCHEMA.md``; this table and its C++
mirror ``common/Isotopes.hh`` both follow it.
"""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Isotope:
    """One beta+ emitter: identity, decay timing, and the positron-range driver."""

    name: str
    """Short label, e.g. ``O15``."""

    Z: int
    """Nuclear charge."""

    A: int
    """Mass number."""

    half_life_s: float
    """Physical half-life [s]."""

    endpoint_MeV: float
    """Beta+ endpoint energy [MeV]; governs the positron range."""

    prompt_gamma: bool
    """The decay emits a prompt de-excitation gamma in coincidence."""

    @property
    def lam(self) -> float:
        """Decay constant lambda = ln2 / T_half  [1/s]."""
        return math.log(2.0) / self.half_life_s


ISOTOPES: dict[int, Isotope] = {
    0: Isotope("O15", Z=8, A=15, half_life_s=122.24, endpoint_MeV=1.74, prompt_gamma=False),
    1: Isotope("C11", Z=6, A=11, half_life_s=1223.4, endpoint_MeV=0.96, prompt_gamma=False),
    2: Isotope("N13", Z=7, A=13, half_life_s=597.9, endpoint_MeV=1.19, prompt_gamma=False),
    3: Isotope("C10", Z=6, A=10, half_life_s=19.29, endpoint_MeV=1.91, prompt_gamma=True),
    4: Isotope("O14", Z=8, A=14, half_life_s=70.62, endpoint_MeV=1.81, prompt_gamma=True),
}
"""The table, keyed by ``isotope_id``: 0=15O, 1=11C, 2=13N, 3=10C, 4=14O."""

N_ISOTOPES = len(ISOTOPES)
"""Number of listed emitters."""

NAME_TO_ID = {iso.name: iid for iid, iso in ISOTOPES.items()}
"""Reverse lookup, name -> ``isotope_id``."""
