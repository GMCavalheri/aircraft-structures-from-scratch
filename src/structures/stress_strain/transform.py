"""Stress/strain transformation and principal values.

Voigt order: 2D ``[sx, sy, txy]``, 3D ``[sx, sy, sz, tyz, txz, txy]``. Strain vectors use
engineering shear strain gamma = 2 eps. Rotations are counter-clockwise by ``theta`` (radians)
from the x-axis to the new x'-axis. For the 2D (plane) case:

    sx' = (sx + sy)/2 + (sx - sy)/2 cos 2t + txy sin 2t
    sy' = (sx + sy)/2 - (sx - sy)/2 cos 2t - txy sin 2t
    txy' =            - (sx - sy)/2 sin 2t + txy cos 2t

(Gere & Goodno, *Mechanics of Materials*, Sec. 7.2; Hibbeler, *Mechanics of Materials*, Ch. 9).
"""

from __future__ import annotations

import numpy as np


def stress_tensor(sx=0.0, sy=0.0, sz=0.0, tyz=0.0, txz=0.0, txy=0.0) -> np.ndarray:
    """Symmetric 3x3 Cauchy stress tensor from its six components."""
    return np.array([[sx, txy, txz], [txy, sy, tyz], [txz, tyz, sz]], dtype=float)


def strain_tensor(ex=0.0, ey=0.0, ez=0.0, gyz=0.0, gxz=0.0, gxy=0.0) -> np.ndarray:
    """Symmetric 3x3 (tensorial) strain tensor from engineering strains (gamma = 2 eps)."""
    return np.array(
        [[ex, gxy / 2, gxz / 2], [gxy / 2, ey, gyz / 2], [gxz / 2, gyz / 2, ez]], dtype=float
    )


def as_tensor(s) -> np.ndarray:
    """Accept a 3x3 tensor, a 2D Voigt vector (plane stress) or a 3D Voigt vector."""
    s = np.asarray(s, dtype=float)
    if s.shape == (3, 3):
        return s
    if s.shape == (3,):
        return stress_tensor(sx=s[0], sy=s[1], txy=s[2])
    if s.shape == (6,):
        return stress_tensor(*s)
    raise ValueError(f"expected a 3x3 tensor or a Voigt vector of length 3 or 6, got {s.shape}")


def rotation_2d(theta: float) -> np.ndarray:
    """Direction cosines [[l1, m1], [l2, m2]] of axes rotated by ``theta`` about z."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s], [-s, c]])


def transform_tensor(sigma, R) -> np.ndarray:
    """sigma' = R sigma R^T for a 3x3 rotation matrix whose rows are the new axes."""
    sigma, R = np.asarray(sigma, float), np.asarray(R, float)
    return R @ sigma @ R.T


def transform_stress_2d(s, theta):
    """Plane stress ``[sx, sy, txy]`` in axes rotated by ``theta``. Broadcasts over ``theta``."""
    sx, sy, txy = (np.asarray(v, float) for v in s)
    theta = np.asarray(theta, float)
    avg, half = (sx + sy) / 2, (sx - sy) / 2
    c, n = np.cos(2 * theta), np.sin(2 * theta)
    return np.array([avg + half * c + txy * n, avg - half * c - txy * n, -half * n + txy * c])


def transform_strain_2d(e, theta):
    """Plane strain ``[ex, ey, gxy]`` (engineering shear) in axes rotated by ``theta``."""
    ex, ey, gxy = (np.asarray(v, float) for v in e)
    out = transform_stress_2d([ex, ey, gxy / 2], theta)
    out[2] *= 2
    return out


def principal_2d(s) -> tuple[float, float, float]:
    """In-plane principal stresses (s1 >= s2) and the angle ``theta_p`` of the s1 direction.

    s1,2 = (sx + sy)/2 +- sqrt(((sx - sy)/2)^2 + txy^2),  tan 2 theta_p = 2 txy / (sx - sy).
    """
    sx, sy, txy = (float(v) for v in s)
    avg, radius = (sx + sy) / 2, float(np.hypot((sx - sy) / 2, txy))
    theta_p = 0.5 * float(np.arctan2(2 * txy, sx - sy))
    return avg + radius, avg - radius, theta_p


def max_shear_2d(s) -> tuple[float, float]:
    """Maximum in-plane shear stress (Mohr's circle radius) and the angle at which it acts."""
    sx, sy, txy = (float(v) for v in s)
    _, _, theta_p = principal_2d(s)
    return float(np.hypot((sx - sy) / 2, txy)), theta_p - np.pi / 4


def principal_3d(s) -> tuple[np.ndarray, np.ndarray]:
    """Principal stresses sorted s1 >= s2 >= s3 and unit principal directions (columns)."""
    w, v = np.linalg.eigh(as_tensor(s))
    order = np.argsort(w)[::-1]
    return w[order], v[:, order]


def stress_invariants(s) -> tuple[float, float, float]:
    """I1 = tr(sigma), I2 = (tr^2 - tr(sigma^2)) / 2, I3 = det(sigma)."""
    t = as_tensor(s)
    tr = np.trace(t)
    return float(tr), float((tr**2 - np.trace(t @ t)) / 2), float(np.linalg.det(t))
