"""Element matrices for 2D bar (truss) and Euler-Bernoulli frame elements.

Every node carries three DOFs (u, v, theta). Local element DOF order is
[u1, v1, theta1, u2, v2, theta2], with local x running from node 1 to node 2. A bar element
uses the same 6x6 layout with zero bending terms, so bars and frames assemble identically
(Cook et al., *Concepts and Applications of Finite Element Analysis*, Chs. 2 and 4;
Przemieniecki, *Theory of Matrix Structural Analysis*, Ch. 5).
"""

from __future__ import annotations

import numpy as np


def rotation(c: float, s: float) -> np.ndarray:
    """6x6 transformation T with u_local = T u_global for direction cosines (c, s)."""
    r = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    T = np.zeros((6, 6))
    T[:3, :3] = r
    T[3:, 3:] = r
    return T


def bar_stiffness(E: float, A: float, L: float) -> np.ndarray:
    """Local stiffness of an axial bar, EA/L [[1, -1], [-1, 1]] on the u DOFs."""
    k = np.zeros((6, 6))
    k[np.ix_([0, 3], [0, 3])] = E * A / L * np.array([[1, -1], [-1, 1]])
    return k


def frame_stiffness(E: float, A: float, I: float, L: float) -> np.ndarray:
    """Local stiffness of a 2D Euler-Bernoulli frame element (axial + cubic Hermite bending)."""
    k = bar_stiffness(E, A, L)
    b = (
        E
        * I
        / L**3
        * np.array(
            [
                [12, 6 * L, -12, 6 * L],
                [6 * L, 4 * L**2, -6 * L, 2 * L**2],
                [-12, -6 * L, 12, -6 * L],
                [6 * L, 2 * L**2, -6 * L, 4 * L**2],
            ]
        )
    )
    k[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = b
    return k


def frame_geometric_stiffness(N: float, L: float) -> np.ndarray:
    """Consistent geometric stiffness of a frame element under axial force N (tension +).

    K_G = N / (30 L) [[36, 3L, -36, 3L], [3L, 4L^2, -3L, -L^2], ...] on the bending DOFs, plus
    N / L [[1, -1], [-1, 1]] on the axial DOFs (Cook et al., Sec. 18.2).
    """
    k = np.zeros((6, 6))
    k[np.ix_([0, 3], [0, 3])] = N / L * np.array([[1, -1], [-1, 1]])
    g = (
        N
        / (30 * L)
        * np.array(
            [
                [36, 3 * L, -36, 3 * L],
                [3 * L, 4 * L**2, -3 * L, -(L**2)],
                [-36, -3 * L, 36, -3 * L],
                [3 * L, -(L**2), -3 * L, 4 * L**2],
            ]
        )
    )
    k[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = g
    return k


def consistent_load(L: float, qy1: float, qy2: float, qx1: float = 0.0, qx2: float = 0.0):
    """Work-equivalent local nodal loads for linearly varying distributed loads (N/m):
    transverse qy from qy1 (node 1) to qy2 (node 2), axial qx likewise.

    Uniform qy gives [0, qL/2, qL^2/12, 0, qL/2, -qL^2/12].
    """
    return np.array(
        [
            L * (2 * qx1 + qx2) / 6,
            L * (7 * qy1 + 3 * qy2) / 20,
            L**2 * (3 * qy1 + 2 * qy2) / 60,
            L * (qx1 + 2 * qx2) / 6,
            L * (3 * qy1 + 7 * qy2) / 20,
            -(L**2) * (2 * qy1 + 3 * qy2) / 60,
        ]
    )


def hermite(xi):
    """Cubic Hermite shape functions at xi in [0, 1] for [v1, L theta1, v2, L theta2]."""
    xi = np.asarray(xi, float)
    return np.array(
        [1 - 3 * xi**2 + 2 * xi**3, xi - 2 * xi**2 + xi**3, 3 * xi**2 - 2 * xi**3, -(xi**2) + xi**3]
    )
