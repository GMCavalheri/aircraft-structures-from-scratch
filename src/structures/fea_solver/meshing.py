"""Mesh helpers."""

from __future__ import annotations

import numpy as np

from structures.fea_solver.model import Model


def beam_model(
    length: float, n_elements: int, E: float, A: float, I: float, angle: float = 0.0
) -> Model:
    """Straight line of ``n_elements`` equal frame elements from the origin, inclined by
    ``angle`` (radians). Nodes are numbered 0..n from the start."""
    m = Model()
    for x in np.linspace(0, length, n_elements + 1):
        m.add_node(x * np.cos(angle), x * np.sin(angle))
    for i in range(n_elements):
        m.add_frame(i, i + 1, E, A, I)
    return m
