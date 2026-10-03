"""Synthetic variable-amplitude load spectra for transport-wing-like stress histories.

This is **not** a standardised spectrum (such as TWIST or FALSTAFF). It is a simple, seeded
generator that has the main features of a lower-wing-skin history:

* a ground-air-ground (GAG) cycle each flight: from the on-ground stress (wing hanging down,
  negative) up to the 1 g in-flight mean;
* gust/manoeuvre increments about the 1 g mean. Their magnitudes follow an exponential
  exceedance law, so the number of increments larger than s falls off as exp(-s / beta),
  which looks like a straight line on the semi-log exceedance plots used for gust spectra.

All stresses are in Pa. Used to exercise rainflow counting and Miner's rule, and as a stand-in
for the load histories of a structural-health-monitoring study.
"""

from __future__ import annotations

import numpy as np


def flight_history(
    n_flights: int,
    s_1g: float,
    s_ground: float,
    beta: float,
    cycles_per_flight: int = 40,
    seed: int = 0,
) -> np.ndarray:
    """Reversal sequence of ``n_flights`` flights.

    Args:
        s_1g: in-flight 1 g mean stress.
        s_ground: on-ground stress (usually negative for a lower wing skin).
        beta: scale of the exponential distribution of gust increments (mean increment).
        cycles_per_flight: gust cycles per flight (each an up and a down excursion).
    """
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_flights):
        out.append(s_ground)
        up = rng.exponential(beta, cycles_per_flight)
        down = rng.exponential(beta, cycles_per_flight)
        for a, b in zip(up, down, strict=True):
            out.extend([s_1g + a, s_1g - b])
        out.append(s_ground)
    return np.array(out)
