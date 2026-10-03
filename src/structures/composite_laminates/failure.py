"""Lamina failure criteria expressed as strength ratios.

The strength ratio SR is the factor by which the applied stress state can be multiplied
before the criterion reaches failure (SR < 1 means failed). Strengths are positive numbers.
Formulas follow Kaw, *Mechanics of Composite Materials*, Sec. 2.8:

* Maximum stress / maximum strain: the smallest ratio of allowable to applied component,
  choosing tension or compression allowables by sign. Ultimate strains are strengths divided
  by the corresponding modulus.
* Tsai-Hill: (s1/X)^2 - s1 s2 / X^2 + (s2/Y)^2 + (t12/S)^2 = 1. The original form uses the
  tensile strengths X = Xt, Y = Yt whatever the signs; the *modified* Tsai-Hill picks tensile
  or compressive strengths according to the signs of s1 and s2.
* Tsai-Wu: F1 s1 + F2 s2 + F11 s1^2 + F22 s2^2 + F66 t12^2 + 2 F12 s1 s2 = 1, with
  F1 = 1/Xt - 1/Xc, F11 = 1/(Xt Xc), F2 = 1/Yt - 1/Yc, F22 = 1/(Yt Yc), F66 = 1/S^2 and
  F12 = F12* sqrt(F11 F22); F12* = -1/2 is Kaw's (Mises-Hencky) default. Solving
  a SR^2 + b SR - 1 = 0 gives SR = (-b + sqrt(b^2 + 4 a)) / (2 a).
"""

from __future__ import annotations

import numpy as np

from structures.materials import OrthotropicPly


def _ratio(value, tension, compression):
    if value > 0:
        return tension / value
    if value < 0:
        return compression / -value
    return np.inf


def max_stress(sigma12, ply: OrthotropicPly) -> tuple[float, str]:
    """(SR, mode) with mode one of 1T, 1C, 2T, 2C, 12S."""
    s1, s2, t = sigma12
    cand = {
        "1T" if s1 >= 0 else "1C": _ratio(s1, ply.Xt, ply.Xc),
        "2T" if s2 >= 0 else "2C": _ratio(s2, ply.Yt, ply.Yc),
        "12S": ply.S / abs(t) if t else np.inf,
    }
    mode = min(cand, key=cand.get)
    return float(cand[mode]), mode


def max_strain(eps12, ply: OrthotropicPly) -> tuple[float, str]:
    e1, e2, g = eps12
    cand = {
        "1T" if e1 >= 0 else "1C": _ratio(e1, ply.Xt / ply.E1, ply.Xc / ply.E1),
        "2T" if e2 >= 0 else "2C": _ratio(e2, ply.Yt / ply.E2, ply.Yc / ply.E2),
        "12S": ply.S / ply.G12 / abs(g) if g else np.inf,
    }
    mode = min(cand, key=cand.get)
    return float(cand[mode]), mode


def tsai_hill(sigma12, ply: OrthotropicPly, modified: bool = False) -> float:
    s1, s2, t = sigma12
    X = ply.Xt if (s1 >= 0 or not modified) else ply.Xc
    Y = ply.Yt if (s2 >= 0 or not modified) else ply.Yc
    index = (s1 / X) ** 2 - s1 * s2 / X**2 + (s2 / Y) ** 2 + (t / ply.S) ** 2
    return float(np.inf if index == 0 else 1 / np.sqrt(index))


def tsai_wu_coefficients(ply: OrthotropicPly, f12_star: float = -0.5):
    F1 = 1 / ply.Xt - 1 / ply.Xc
    F2 = 1 / ply.Yt - 1 / ply.Yc
    F11 = 1 / (ply.Xt * ply.Xc)
    F22 = 1 / (ply.Yt * ply.Yc)
    F66 = 1 / ply.S**2
    F12 = f12_star * np.sqrt(F11 * F22)
    return F1, F2, F11, F22, F66, F12


def tsai_wu(sigma12, ply: OrthotropicPly, f12_star: float = -0.5) -> float:
    s1, s2, t = sigma12
    F1, F2, F11, F22, F66, F12 = tsai_wu_coefficients(ply, f12_star)
    a = F11 * s1**2 + F22 * s2**2 + F66 * t**2 + 2 * F12 * s1 * s2
    b = F1 * s1 + F2 * s2
    if a == 0:
        return float(np.inf if b <= 0 else 1 / b)
    return float((-b + np.sqrt(b * b + 4 * a)) / (2 * a))


def tsai_wu_index(sigma12, ply: OrthotropicPly, f12_star: float = -0.5) -> float:
    """Left-hand side of the Tsai-Wu criterion (failure when >= 1)."""
    s1, s2, t = sigma12
    F1, F2, F11, F22, F66, F12 = tsai_wu_coefficients(ply, f12_star)
    return float(F1 * s1 + F2 * s2 + F11 * s1**2 + F22 * s2**2 + F66 * t**2 + 2 * F12 * s1 * s2)


CRITERIA = ("tsai_wu", "tsai_hill", "modified_tsai_hill", "max_stress", "max_strain")
