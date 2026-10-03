"""Palmgren-Miner linear cumulative damage: D = sum n_i / N_i, failure predicted at D = 1."""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np

from structures.fatigue.rainflow import Cycle


def miner_damage(cycles: Iterable[Cycle], life: Callable) -> float:
    """Damage of rainflow ``cycles``; ``life(amplitude, mean)`` returns cycles to failure
    (``inf`` for non-damaging cycles)."""
    cyc = list(cycles)
    if not cyc:
        return 0.0
    amp = np.array([c.amplitude for c in cyc])
    mean = np.array([c.mean for c in cyc])
    n = np.array([c.count for c in cyc])
    N = np.asarray(life(amp, mean), float)
    return float(np.sum(n / N))


def blocks_to_failure(damage_per_block: float, critical_damage: float = 1.0) -> float:
    return np.inf if damage_per_block == 0 else critical_damage / damage_per_block
