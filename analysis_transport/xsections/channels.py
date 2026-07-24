"""Frozen channel definitions shared by readers, simulation and plots."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Channel:
    channel_id: str
    target: str
    target_z: int
    target_a: int
    residual: str
    residual_z: int
    residual_a: int
    title: str


CHANNELS = (
    Channel("p_C12_x_C11", "C12", 6, 12, "C11", 6, 11, "C-12(p,x)C-11"),
    Channel("p_O16_x_O15", "O16", 8, 16, "O15", 8, 15, "O-16(p,x)O-15"),
    Channel("p_O16_x_C11", "O16", 8, 16, "C11", 6, 11, "O-16(p,x)C-11"),
    Channel("p_N14_x_N13", "N14", 7, 14, "N13", 7, 13, "N-14(p,x)N-13"),
    Channel("p_O16_x_N13", "O16", 8, 16, "N13", 7, 13, "O-16(p,x)N-13"),
)

BY_ID = {channel.channel_id: channel for channel in CHANNELS}
BY_PAIR = {(channel.target, channel.residual): channel for channel in CHANNELS}


def channel_for(target: str, residual: str) -> Channel:
    return BY_PAIR[(target, residual)]
