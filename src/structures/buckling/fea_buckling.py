"""Linear (eigenvalue) buckling with the from-scratch FEA solver.

1. Solve the linear static problem under a reference load and read each element's axial
   force N_e (tension positive).
2. Assemble the geometric stiffness K_G(N) and solve (K + lambda K_G) phi = 0 on the free DOFs.

We solve -K_G phi = mu K phi, with K symmetric positive definite once the supports are
applied, so ``scipy.linalg.eigh`` applies directly. Then lambda = 1 / mu, and the smallest
positive load factor comes from the largest positive mu.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh

from structures.fea_solver import Model, beam_model


@dataclass
class BucklingResult:
    load_factors: np.ndarray  # ascending positive critical load multipliers
    modes: np.ndarray  # (n_dof, n_modes) full displacement vectors, max |v| = 1
    model: Model


def linear_buckling(model: Model, n_modes: int = 4) -> BucklingResult:
    static = model.solve()
    N = np.array([static.axial_force(i) for i in range(len(model.elements))])
    K = model.stiffness().toarray()
    KG = model.geometric_stiffness(N).toarray()
    restrained = model.restrained_dofs()
    free = np.setdiff1d(np.arange(model.n_dof), np.array(sorted(restrained), dtype=int))
    mu, phi = eigh(-KG[np.ix_(free, free)], K[np.ix_(free, free)])
    pos = mu > 1e-12 * np.abs(mu).max()
    order = np.argsort(mu[pos])[::-1][:n_modes]
    lam = 1.0 / mu[pos][order]
    modes = np.zeros((model.n_dof, len(lam)))
    modes[free] = phi[:, pos][:, order]
    for j in range(modes.shape[1]):
        v = modes[1::3, j]
        modes[:, j] /= v[np.argmax(np.abs(v))]
    return BucklingResult(lam, modes, model)


def column_model(
    length: float,
    n_elements: int,
    E: float,
    A: float,
    I: float,
    end: str = "pinned-pinned",
    load: float = 1.0,
) -> Model:
    """Column along x loaded by compressive force ``load`` at x = L. End conditions:
    ``pinned-pinned``, ``fixed-free``, ``fixed-fixed``, ``fixed-pinned`` (fixed end at x=0)."""
    m = beam_model(length, n_elements, E, A, I)
    n = n_elements
    supports = {
        "pinned-pinned": (("xy", 0), ("y", n)),
        "fixed-free": (("xyr", 0),),
        "fixed-fixed": (("xyr", 0), ("yr", n)),
        "fixed-pinned": (("xyr", 0), ("y", n)),
    }
    if end not in supports:
        raise ValueError(f"unknown end condition {end!r}")
    for dofs, node in supports[end]:
        m.fix(node, dofs)
    m.load(n, fx=-load)
    return m
