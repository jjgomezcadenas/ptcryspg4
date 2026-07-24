"""Calculate selected laboratory reaction thresholds from AME2020 masses."""

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .config import load


U_C2_MEV = 931.49410242
MASS_SOURCE = "AME2020 atomic-mass evaluation"
MASS_SOURCE_URL = "https://amdc.impcas.ac.cn/web/masseval.html"

# Neutral-atom masses in u.  Electron numbers balance in both reactions.
ATOMIC_MASS_U = {
    "H1": 1.007825031898,
    "H2": 2.014101777844,
    "He4": 4.00260325413,
    "C11": 11.011432597,
    "O15": 15.003065636,
    "O16": 15.99491461957,
}


@dataclass(frozen=True)
class ThresholdReaction:
    channel_id: str
    reaction: str
    projectile: str
    target: str
    products: tuple[str, ...]


REACTIONS = (
    ThresholdReaction(
        "p_O16_x_O15", "O-16(p,d)O-15", "H1", "O16", ("H2", "O15")),
    ThresholdReaction(
        "p_O16_x_C11", "O-16(p,d+alpha)C-11", "H1", "O16",
        ("H2", "He4", "C11")),
)


def q_value_mev(reaction: ThresholdReaction) -> float:
    initial = ATOMIC_MASS_U[reaction.projectile] + ATOMIC_MASS_U[reaction.target]
    final = sum(ATOMIC_MASS_U[product] for product in reaction.products)
    return (initial - final) * U_C2_MEV


def laboratory_threshold_mev(reaction: ThresholdReaction) -> float:
    """Exact invariant-mass threshold for a projectile on a target at rest."""
    projectile = ATOMIC_MASS_U[reaction.projectile]
    target = ATOMIC_MASS_U[reaction.target]
    final = sum(ATOMIC_MASS_U[product] for product in reaction.products)
    return ((final * final - (projectile + target) ** 2)
            / (2.0 * target) * U_C2_MEV)


def threshold_table() -> pd.DataFrame:
    return pd.DataFrame([{
        "channel_id": reaction.channel_id,
        "lowest_channel": reaction.reaction,
        "q_value_MeV": q_value_mev(reaction),
        "laboratory_threshold_MeV": laboratory_threshold_mev(reaction),
        "mass_source": MASS_SOURCE,
        "mass_source_url": MASS_SOURCE_URL,
    } for reaction in REACTIONS])


def generate(repo: Path, config_path: Path) -> pd.DataFrame:
    table = threshold_table()
    config = load(config_path)
    configured = config["threshold_MeV"]
    for row in table.itertuples(index=False):
        if abs(float(configured[row.channel_id])
               - row.laboratory_threshold_MeV) > 0.01:
            raise ValueError(
                f"Configured threshold for {row.channel_id} differs from "
                f"the AME2020 calculation: {configured[row.channel_id]} versus "
                f"{row.laboratory_threshold_MeV:.5f} MeV")
    output = repo / "data/xsections/reaction_thresholds.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False, lineterminator="\n")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--config", type=Path,
                        default=Path("config/xsection_fit.toml"))
    args = parser.parse_args()
    table = generate(args.repo.resolve(), args.config.resolve())
    for row in table.itertuples(index=False):
        print(f"{row.lowest_channel}: Q = {row.q_value_MeV:.5f} MeV; "
              f"laboratory threshold = {row.laboratory_threshold_MeV:.5f} MeV")


if __name__ == "__main__":
    main()
