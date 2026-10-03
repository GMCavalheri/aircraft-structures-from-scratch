"""Plane-stress stiffness of a unidirectional lamina and its rotation to laminate axes.

Material axes: 1 along the fibres, 2 transverse. Ply angle theta is measured from the laminate
x-axis to the fibre direction, counter-clockwise positive. Stress/strain vectors are
[s1, s2, t12] and [e1, e2, g12] (engineering shear). Notation follows Kaw, *Mechanics of
Composite Materials*, 2nd ed., Ch. 2:

    Q11 = E1 / (1 - nu12 nu21),  Q12 = nu12 E2 / (1 - nu12 nu21),
    Q22 = E2 / (1 - nu12 nu21),  Q66 = G12,

    [sigma]_12 = [T] [sigma]_xy,   Qbar = T^-1 Q R T R^-1,  R = diag(1, 1, 2).
"""

from __future__ import annotations

import numpy as np

from structures.materials import OrthotropicPly


def reduced_stiffness(ply: OrthotropicPly) -> np.ndarray:
    d = 1 - ply.nu12 * ply.nu21
    return np.array(
        [
            [ply.E1 / d, ply.nu12 * ply.E2 / d, 0.0],
            [ply.nu12 * ply.E2 / d, ply.E2 / d, 0.0],
            [0.0, 0.0, ply.G12],
        ]
    )


def reduced_compliance(ply: OrthotropicPly) -> np.ndarray:
    return np.array(
        [
            [1 / ply.E1, -ply.nu12 / ply.E1, 0.0],
            [-ply.nu12 / ply.E1, 1 / ply.E2, 0.0],
            [0.0, 0.0, 1 / ply.G12],
        ]
    )


def stress_transformation(theta: float) -> np.ndarray:
    """[T] mapping global stresses [sx, sy, txy] to material stresses [s1, s2, t12]."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array(
        [[c * c, s * s, 2 * s * c], [s * s, c * c, -2 * s * c], [-s * c, s * c, c * c - s * s]]
    )


_R = np.diag([1.0, 1.0, 2.0])
_R_INV = np.diag([1.0, 1.0, 0.5])


def strain_transformation(theta: float) -> np.ndarray:
    """Maps global engineering strains [ex, ey, gxy] to [e1, e2, g12]: R T R^-1."""
    return _R @ stress_transformation(theta) @ _R_INV


def qbar(ply: OrthotropicPly, theta: float) -> np.ndarray:
    """Transformed reduced stiffness of a ply at angle ``theta`` (radians)."""
    T = stress_transformation(theta)
    return np.linalg.solve(T, reduced_stiffness(ply)) @ strain_transformation(theta)


def invariants(ply: OrthotropicPly) -> np.ndarray:
    """Tsai-Pagano stiffness invariants U1..U5 (Kaw Eq. 2.131)."""
    Q = reduced_stiffness(ply)
    q11, q12, q22, q66 = Q[0, 0], Q[0, 1], Q[1, 1], Q[2, 2]
    return np.array(
        [
            (3 * q11 + 3 * q22 + 2 * q12 + 4 * q66) / 8,
            (q11 - q22) / 2,
            (q11 + q22 - 2 * q12 - 4 * q66) / 8,
            (q11 + q22 + 6 * q12 - 4 * q66) / 8,
            (q11 + q22 - 2 * q12 + 4 * q66) / 8,
        ]
    )


def qbar_from_invariants(U: np.ndarray, theta: float) -> np.ndarray:
    """Qbar from the invariants: Qbar11 = U1 + U2 cos 2t + U3 cos 4t, etc."""
    U1, U2, U3, U4, U5 = U
    c2, c4, s2, s4 = np.cos(2 * theta), np.cos(4 * theta), np.sin(2 * theta), np.sin(4 * theta)
    q11 = U1 + U2 * c2 + U3 * c4
    q12 = U4 - U3 * c4
    q22 = U1 - U2 * c2 + U3 * c4
    q16 = U2 * s2 / 2 + U3 * s4
    q26 = U2 * s2 / 2 - U3 * s4
    q66 = U5 - U3 * c4
    return np.array([[q11, q12, q16], [q12, q22, q26], [q16, q26, q66]])
