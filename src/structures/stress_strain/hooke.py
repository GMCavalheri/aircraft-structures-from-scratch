"""Generalised Hooke's law for a linear-elastic isotropic solid (Voigt notation,
engineering shear strains)."""

from __future__ import annotations

import numpy as np


def compliance_matrix(E: float, nu: float) -> np.ndarray:
    """6x6 compliance [S] with strain = S stress, order [x, y, z, yz, xz, xy]."""
    G = E / (2 * (1 + nu))
    S = np.zeros((6, 6))
    S[:3, :3] = -nu / E
    np.fill_diagonal(S[:3, :3], 1 / E)
    S[3:, 3:] = np.eye(3) / G
    return S


def stiffness_matrix(E: float, nu: float) -> np.ndarray:
    """6x6 stiffness [C] = [S]^-1, written out with the Lame constants."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[:3, :3] += 2 * G * np.eye(3)
    C[3:, 3:] = G * np.eye(3)
    return C


def plane_stress_stiffness(E: float, nu: float) -> np.ndarray:
    """3x3 [D] for plane stress (sz = 0): [sx, sy, txy] = D [ex, ey, gxy]."""
    return E / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])


def plane_strain_stiffness(E: float, nu: float) -> np.ndarray:
    """3x3 [D] for plane strain (ez = 0)."""
    f = E / ((1 + nu) * (1 - 2 * nu))
    return f * np.array([[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2]])
