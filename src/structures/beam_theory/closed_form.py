"""Closed-form Euler-Bernoulli results for standard load cases, used as validation
references (Gere & Goodno, *Mechanics of Materials*, Appendix G; Roark's *Formulas for Stress
and Strain*, Table 8.1). Magnitudes; P and w are positive *downward* here, and deflections are
returned as positive downward values."""

from __future__ import annotations

import numpy as np


def cantilever_tip_load(P, L, EI):
    """Tip deflection PL^3/3EI, tip slope PL^2/2EI, root moment PL."""
    return {"deflection": P * L**3 / (3 * EI), "slope": P * L**2 / (2 * EI), "moment": P * L}


def cantilever_udl(w, L, EI):
    """Tip deflection wL^4/8EI, tip slope wL^3/6EI, root moment wL^2/2."""
    return {"deflection": w * L**4 / (8 * EI), "slope": w * L**3 / (6 * EI), "moment": w * L**2 / 2}


def simply_supported_center_load(P, L, EI):
    """Midspan deflection PL^3/48EI, end slope PL^2/16EI, max moment PL/4."""
    return {"deflection": P * L**3 / (48 * EI), "slope": P * L**2 / (16 * EI), "moment": P * L / 4}


def simply_supported_udl(w, L, EI):
    """Midspan deflection 5wL^4/384EI, end slope wL^3/24EI, max moment wL^2/8."""
    return {
        "deflection": 5 * w * L**4 / (384 * EI),
        "slope": w * L**3 / (24 * EI),
        "moment": w * L**2 / 8,
    }


def simply_supported_point_load(P, a, L, EI):
    """Load P at distance a from the left support. With b_s = min(a, L - a) the distance to
    the nearer support, the maximum deflection P b_s (L^2 - b_s^2)^1.5 / (9 sqrt(3) L EI) occurs
    sqrt((L^2 - b_s^2) / 3) from the *farther* support. ``x_max`` is measured from the left.
    Reactions P (L - a) / L (left) and P a / L (right); maximum moment P a (L - a) / L."""
    b_s = min(a, L - a)
    x_far = np.sqrt((L**2 - b_s**2) / 3)
    return {
        "deflection": P * b_s * (L**2 - b_s**2) ** 1.5 / (9 * np.sqrt(3) * L * EI),
        "x_max": x_far if a >= L - a else L - x_far,
        "R_left": P * (L - a) / L,
        "R_right": P * a / L,
        "moment": P * a * (L - a) / L,
    }


def simply_supported_triangular(w0, L, EI):
    """Load rising linearly from 0 at the left to w0 at the right: reactions w0L/6 (left) and
    w0L/3 (right); v = w0 x (7L^4 - 10L^2 x^2 + 3x^4) / (360 L EI), maximum at
    x = L sqrt(1 - sqrt(8/15)) = 0.5193 L, where it equals 0.006522 w0 L^4 / EI."""
    x = L * np.sqrt(1 - np.sqrt(8 / 15))
    v = w0 * x * (7 * L**4 - 10 * L**2 * x**2 + 3 * x**4) / (360 * L * EI)
    return {
        "R_left": w0 * L / 6,
        "R_right": w0 * L / 3,
        "deflection": v,
        "x_max": x,
        "moment": w0 * L**2 / (9 * np.sqrt(3)),
    }


def propped_cantilever_udl(w, L, EI):
    """Fixed at x = 0, simple support at x = L: reactions 5wL/8 (fixed), 3wL/8 (prop), fixed-end
    moment wL^2/8; v = w x^2 (3L^2 - 5Lx + 2x^2) / (48 EI), maximum at x = (15 - sqrt(33)) L / 16
    = 0.5785 L, where it equals 0.005416 w L^4 / EI (often quoted as wL^4 / 185 EI)."""
    x = (15 - np.sqrt(33)) * L / 16
    v = w * x**2 * (3 * L**2 - 5 * L * x + 2 * x**2) / (48 * EI)
    return {
        "R_fixed": 5 * w * L / 8,
        "R_prop": 3 * w * L / 8,
        "moment": w * L**2 / 8,
        "deflection": v,
        "x_max": x,
    }


def fixed_fixed_udl(w, L, EI):
    """End moments wL^2/12, midspan moment wL^2/24, midspan deflection wL^4/384EI."""
    return {
        "moment": w * L**2 / 12,
        "moment_mid": w * L**2 / 24,
        "deflection": w * L**4 / (384 * EI),
    }


def fixed_fixed_center_load(P, L, EI):
    """End moments PL/8, midspan deflection PL^3/192EI."""
    return {"moment": P * L / 8, "deflection": P * L**3 / (192 * EI)}


def two_span_continuous_udl(w, L, EI):
    """Two equal spans L on three supports: reactions 3wL/8, 10wL/8, 3wL/8; moment over the
    middle support wL^2/8."""
    return {"R_end": 3 * w * L / 8, "R_mid": 10 * w * L / 8, "moment": w * L**2 / 8}


def simply_supported_end_couple(M0, L, EI):
    """Couple M0 at the right end: end slopes M0 L / 3EI (loaded end), M0 L / 6EI (far end),
    reactions M0 / L."""
    return {"slope_loaded": M0 * L / (3 * EI), "slope_far": M0 * L / (6 * EI), "R": M0 / L}
