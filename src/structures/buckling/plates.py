"""Elastic buckling of flat rectangular plates (aircraft skin panels).

    sigma_cr = k pi^2 E / (12 (1 - nu^2)) (t / b)^2

b is the loaded edge width (compression) or the short side (shear). For a plate simply
supported on all four edges under uniaxial compression, with a the length in the load
direction, buckling into m half-waves gives

    k = (m b / a + a / (m b))^2,   minimised over m = 1, 2, ...,

so k = 4 whenever a/b is an integer and k -> 4 for long plates. Shear buckling, all edges
simply supported, a >= b: k_s = 5.35 + 4 (b / a)^2 (Timoshenko & Gere, *Theory of Elastic
Stability*, Secs. 9.2 and 9.7; Bruhn, Sec. C5).
"""

from __future__ import annotations

import numpy as np


def compression_k(aspect_ratio, max_half_waves: int = 50) -> np.ndarray:
    """Buckling coefficient k for a simply supported plate in uniaxial compression; returns
    (k, m) with the governing number of half-waves m."""
    r = np.atleast_1d(np.asarray(aspect_ratio, float))[:, None]
    m = np.arange(1, max_half_waves + 1)[None, :]
    k = (m / r + r / m) ** 2
    i = np.argmin(k, axis=1)
    kmin, mmin = k[np.arange(len(r)), i], m[0, i]
    if np.ndim(aspect_ratio) == 0:
        return float(kmin[0]), int(mmin[0])
    return kmin, mmin


def shear_k(aspect_ratio: float) -> float:
    """Shear buckling coefficient, simply supported edges (a/b >= 1 approximation)."""
    r = max(aspect_ratio, 1 / aspect_ratio)
    return 5.35 + 4 / r**2


def plate_critical_stress(E: float, nu: float, t: float, b: float, k: float) -> float:
    return k * np.pi**2 * E / (12 * (1 - nu**2)) * (t / b) ** 2
