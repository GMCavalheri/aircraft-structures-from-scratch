"""Column buckling: Euler, Johnson parabola and tangent-modulus (Engesser) columns.

Euler (MIL-HDBK-5J Eq. 1.3.8(b)):  P_cr = pi^2 E I / (K L)^2,  F_c = pi^2 E / (L'/rho)^2.

Effective-length factors K are the theoretical values for ideal end conditions; fixed-pinned
is pi / 4.4934 = 0.6992, where 4.4934 is the first root of tan(x) = x (Timoshenko & Gere,
*Theory of Elastic Stability*, Ch. 2).

Johnson parabola for intermediate columns (Bruhn, *Analysis and Design of Flight Vehicle
Structures*, Sec. C2; Megson, Ch. 8):

    F_c = F_cy - F_cy^2 / (4 pi^2 E) (L'/rho)^2,   valid for L'/rho < pi sqrt(2 E / F_cy),

where it meets the Euler curve tangentially at F_c = F_cy / 2.

Tangent-modulus column (MIL-HDBK-5J Eq. 1.3.8(a)): F_c = pi^2 E_t / (L'/rho)^2 with E_t from the
Ramberg-Osgood curve of MIL-HDBK-5J Eq. 1.3.9: e = f / E + 0.002 (f / F_0.2)^n, so
1 / E_t = 1 / E + 0.002 n f^(n-1) / F_0.2^n.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

EFFECTIVE_LENGTH = {
    "pinned-pinned": 1.0,
    "fixed-free": 2.0,
    "fixed-fixed": 0.5,
    "fixed-pinned": np.pi / 4.493409457909064,
}


def effective_length_factor(end: str) -> float:
    try:
        return EFFECTIVE_LENGTH[end]
    except KeyError:
        raise ValueError(f"end condition must be one of {list(EFFECTIVE_LENGTH)}") from None


def euler_load(E: float, I: float, L: float, end: str = "pinned-pinned") -> float:
    """Critical load pi^2 E I / (K L)^2."""
    return np.pi**2 * E * I / (effective_length_factor(end) * L) ** 2


def euler_stress(E: float, slenderness):
    """Euler critical stress pi^2 E / (L'/rho)^2 for effective slenderness L'/rho."""
    return np.pi**2 * E / np.asarray(slenderness, float) ** 2


def johnson_transition(E: float, Fcy: float) -> float:
    """Slenderness at which the Johnson parabola meets the Euler hyperbola."""
    return float(np.pi * np.sqrt(2 * E / Fcy))


def johnson_stress(E: float, Fcy: float, slenderness):
    """Johnson parabola below the transition slenderness, Euler above it."""
    lam = np.asarray(slenderness, float)
    johnson = Fcy - Fcy**2 / (4 * np.pi**2 * E) * lam**2
    with np.errstate(divide="ignore"):
        euler = euler_stress(E, lam)
    return np.where(lam < johnson_transition(E, Fcy), johnson, euler)


def tangent_modulus(stress, E: float, F02: float, n: float):
    """Ramberg-Osgood tangent modulus E_t at stress ``stress``."""
    f = np.asarray(stress, float)
    return 1.0 / (1.0 / E + 0.002 * n * f ** (n - 1) / F02**n)


def tangent_modulus_stress(E: float, F02: float, n: float, slenderness):
    """Column stress F_c solving F_c = pi^2 E_t(F_c) / (L'/rho)^2 (Engesser)."""
    lam = np.atleast_1d(np.asarray(slenderness, float))
    out = np.empty_like(lam)
    for i, s in enumerate(lam):
        g = lambda f, s=s: f - np.pi**2 * tangent_modulus(f, E, F02, n) / s**2  # noqa: E731
        hi = min(np.pi**2 * E / s**2, 3 * F02)  # the elastic Euler stress bounds the root
        out[i] = brentq(g, 1e-9 * F02, hi) if g(hi) > 0 else hi
    return out if np.ndim(slenderness) else float(out[0])
