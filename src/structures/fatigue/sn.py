"""Stress-life (S-N) curves.

* Basquin power law N = C S^-m, fitted by least squares on log N versus log S.
* MIL-HDBK-5J equivalent-stress model, which folds every stress ratio R onto one curve:

      log10 Nf = A1 + A2 log10(Seq - A4),   Seq = Smax (1 - R)^A3

  A4 acts as a fatigue limit: lives are infinite for Seq <= A4. Since Sa = Smax (1 - R) / 2,
  Seq = 2^A3 Smax^(1 - A3) Sa^A3 is a Walker-type mean-stress correction with exponent A3
  (MIL-HDBK-5J Secs. 1.4.9.2 and 9.3.4; Walker, ASTM STP 462, 1970).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from structures.materials import KSI

DATA = Path(__file__).resolve().parent / "data"


@dataclass(frozen=True)
class BasquinCurve:
    """N = C S^-m (S in Pa, or any consistent unit; S is an amplitude or range as fitted)."""

    C: float
    m: float

    def life(self, S):
        return self.C * np.asarray(S, float) ** (-self.m)

    def strength(self, N):
        return (self.C / np.asarray(N, float)) ** (1 / self.m)


def fit_basquin(S, N) -> BasquinCurve:
    """Least-squares fit of log10 N = log10 C - m log10 S (life is the dependent variable,
    as in ASTM E739)."""
    x, y = np.log10(np.asarray(S, float)), np.log10(np.asarray(N, float))
    slope, intercept = np.polyfit(x, y, 1)
    return BasquinCurve(C=10**intercept, m=-slope)


@dataclass(frozen=True)
class EquivalentStressCurve:
    """MIL-HDBK-5J best-fit S/N equation (coefficients in ksi as printed)."""

    A1: float
    A2: float
    A3: float
    A4: float  # ksi
    std_error: float = 0.0  # standard error of log10 life
    name: str = ""
    source: str = ""

    def equivalent_stress(self, s_max, R):
        """Seq [Pa] for maximum stress ``s_max`` [Pa] and stress ratio ``R`` (R = Smin/Smax)."""
        return np.asarray(s_max, float) * (1 - np.asarray(R, float)) ** self.A3

    def life(self, s_max, R):
        """Median cycles to failure; ``inf`` at or below the fatigue limit A4."""
        seq = self.equivalent_stress(s_max, R) / KSI
        excess = seq - self.A4
        with np.errstate(divide="ignore", invalid="ignore"):
            logn = self.A1 + self.A2 * np.log10(np.where(excess > 0, excess, np.nan))
        return np.where(excess > 0, 10.0**logn, np.inf)

    def life_from_amplitude(self, s_a, s_m):
        """Life for amplitude ``s_a`` and mean ``s_m`` [Pa]. Cycles with Smax <= 0 (wholly
        compressive) are treated as non-damaging."""
        s_a, s_m = np.asarray(s_a, float), np.asarray(s_m, float)
        s_max = s_m + s_a
        with np.errstate(divide="ignore", invalid="ignore"):
            R = np.where(s_max > 0, (s_m - s_a) / s_max, 0.0)
        return np.where(s_max > 0, self.life(np.maximum(s_max, 0), R), np.inf)

    def max_stress_for_life(self, N, R):
        """Inverse: Smax [Pa] giving median life ``N`` at stress ratio ``R``."""
        seq = 10 ** ((np.log10(np.asarray(N, float)) - self.A1) / self.A2) + self.A4
        return seq * KSI / (1 - np.asarray(R, float)) ** self.A3


@cache
def handbook_curve(key: str) -> EquivalentStressCurve:
    """``"2024-T3_Kt1"``, ``"2024-T3_Kt2"``, ``"7075-T6_Kt1"`` or ``"7075-T6_Kt2"``."""
    with open(DATA / "mil_hdbk_5j_sn.csv") as f:
        for row in csv.DictReader(line for line in f if not line.startswith("#")):
            if row["key"] == key:
                return EquivalentStressCurve(
                    A1=float(row["A1"]),
                    A2=float(row["A2"]),
                    A3=float(row["A3"]),
                    A4=float(row["A4_ksi"]),
                    std_error=float(row["std_error_log_life"]),
                    name=f"{row['material']} Kt = {row['Kt']}",
                    source=row["source"],
                )
    raise KeyError(key)


def fit_equivalent_stress(s_max, R, N, guess=(10.0, -4.0, 0.5, 10.0)) -> EquivalentStressCurve:
    """Fit A1..A4 by nonlinear least squares on log10 life (s_max in Pa, A4 in ksi). Runouts
    must be removed beforehand."""
    s_max_ksi = np.asarray(s_max, float) / KSI
    R, logN = np.asarray(R, float), np.log10(np.asarray(N, float))

    def residual(p):
        A1, A2, A3, A4 = p
        excess = np.maximum(s_max_ksi * (1 - R) ** A3 - A4, 1e-9)
        return A1 + A2 * np.log10(excess) - logN

    seq_min = float(np.min(s_max_ksi * (1 - R) ** guess[2]))
    sol = least_squares(
        residual, guess, bounds=([-np.inf, -np.inf, 0.0, 0.0], [np.inf, 0.0, 1.0, 0.999 * seq_min])
    )
    A1, A2, A3, A4 = sol.x
    dof = max(len(logN) - 4, 1)
    return EquivalentStressCurve(
        A1, A2, A3, A4, std_error=float(np.sqrt(np.sum(sol.fun**2) / dof)), name="fitted"
    )
