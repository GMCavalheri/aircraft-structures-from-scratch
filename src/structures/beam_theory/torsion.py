"""Torsion of shafts and thin-walled sections.

* Solid/hollow circular shaft: tau = T r / J, twist phi = T L / (G J).
* Closed thin-walled cell (Bredt-Batho): shear flow q = T / (2 A_enc), tau = q / t,
  twist rate dphi/dx = T / (4 A_enc^2 G) * sum(l_i / t_i).
* Open thin-walled section: J = sum(b t^3) / 3, tau_max = T t_max / J.

(Megson, *Aircraft Structures for Engineering Students*, Chs. 3 and 18.)
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TorsionResult:
    tau_max: float  # maximum shear stress, Pa
    twist_rate: float  # rad per metre
    J: float  # torsion constant, m^4

    def twist(self, length: float) -> float:
        return self.twist_rate * length


def circular_shaft(T: float, G: float, d_outer: float, d_inner: float = 0.0) -> TorsionResult:
    J = np.pi * (d_outer**4 - d_inner**4) / 32
    return TorsionResult(tau_max=T * d_outer / 2 / J, twist_rate=T / (G * J), J=J)


def bredt_batho(
    T: float, G: float, enclosed_area: float, walls: Sequence[tuple[float, float]]
) -> TorsionResult:
    """Single closed cell; ``walls`` is a list of (length, thickness) around the cell."""
    q = T / (2 * enclosed_area)
    line_integral = sum(length / t for length, t in walls)
    t_min = min(t for _, t in walls)
    J = 4 * enclosed_area**2 / line_integral
    return TorsionResult(tau_max=q / t_min, twist_rate=T / (G * J), J=J)


def shear_flow_closed(T: float, enclosed_area: float) -> float:
    """Constant shear flow q = T / (2 A) in a single-cell closed section (N/m)."""
    return T / (2 * enclosed_area)


def open_thin_walled(T: float, G: float, walls: Sequence[tuple[float, float]]) -> TorsionResult:
    """Open section made of thin rectangles ``(length, thickness)``."""
    J = sum(b * t**3 for b, t in walls) / 3
    t_max = max(t for _, t in walls)
    return TorsionResult(tau_max=T * t_max / J, twist_rate=T / (G * J), J=J)
