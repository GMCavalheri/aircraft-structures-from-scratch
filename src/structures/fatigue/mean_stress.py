"""Mean-stress corrections: the fully reversed (R = -1) amplitude S_ar that is as damaging
as amplitude S_a at mean S_m.

    Goodman    S_ar = S_a / (1 - S_m / S_u)
    Gerber     S_ar = S_a / (1 - (S_m / S_u)^2)
    Soderberg  S_ar = S_a / (1 - S_m / S_y)
    SWT        S_ar = sqrt(S_max S_a)              (Smith, Watson & Topper, 1970)
    Walker     S_ar = S_max^(1 - g) S_a^g          (g = 0.5 recovers SWT)

(Dowling, *Mechanical Behavior of Materials*, Secs. 9.6-9.7.) Gerber treats compressive and
tensile means alike, so it should not be used for compressive means.
"""

from __future__ import annotations

import numpy as np


def goodman(s_a, s_m, s_u):
    return np.asarray(s_a, float) / (1 - np.asarray(s_m, float) / s_u)


def gerber(s_a, s_m, s_u):
    return np.asarray(s_a, float) / (1 - (np.asarray(s_m, float) / s_u) ** 2)


def soderberg(s_a, s_m, s_y):
    return np.asarray(s_a, float) / (1 - np.asarray(s_m, float) / s_y)


def smith_watson_topper(s_a, s_m):
    s_a, s_m = np.asarray(s_a, float), np.asarray(s_m, float)
    s_max = s_a + s_m
    return np.sqrt(np.maximum(s_max, 0) * s_a)


def walker(s_a, s_m, gamma):
    s_a, s_m = np.asarray(s_a, float), np.asarray(s_m, float)
    s_max = np.maximum(s_a + s_m, 0)
    return s_max ** (1 - gamma) * s_a**gamma
