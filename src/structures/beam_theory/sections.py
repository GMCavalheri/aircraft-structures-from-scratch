"""Cross-section properties.

Section coordinates: ``z`` horizontal, ``y`` vertical (up), so that ``I`` (= I_z, the second
moment about the horizontal centroidal axis) governs bending in the vertical x-y plane where the
deflection v is measured. ``y_top``/``y_bottom`` are distances from the centroid to the extreme
fibres, both positive.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Section:
    """Geometric properties of a beam cross-section (SI: m, m^2, m^4)."""

    A: float
    I: float
    y_top: float
    y_bottom: float
    J: float = 0.0
    I_lateral: float = 0.0
    I_product: float = 0.0
    name: str = ""

    @property
    def S_top(self) -> float:
        """Elastic section modulus to the top fibre, I / y_top."""
        return self.I / self.y_top

    @property
    def S_bottom(self) -> float:
        return self.I / self.y_bottom

    @property
    def radius_of_gyration(self) -> float:
        return float(np.sqrt(self.I / self.A))


def rectangle(b: float, h: float) -> Section:
    """Solid rectangle of width ``b`` and depth ``h``. Torsion constant from the series
    solution of Saint-Venant's problem (Roark's Formulas, Table 10.7 case 4):
    J = a c^3 [1/3 - 0.21 c/a (1 - c^4 / (12 a^4))], a >= c the long and short sides."""
    a, c = max(b, h), min(b, h)
    J = a * c**3 * (1 / 3 - 0.21 * c / a * (1 - c**4 / (12 * a**4)))
    return Section(
        A=b * h,
        I=b * h**3 / 12,
        y_top=h / 2,
        y_bottom=h / 2,
        J=J,
        I_lateral=h * b**3 / 12,
        name=f"rectangle {b:g} x {h:g}",
    )


def circle(d: float) -> Section:
    """Solid circular section of diameter ``d``: I = pi d^4 / 64, J = pi d^4 / 32."""
    I = np.pi * d**4 / 64
    return Section(
        A=np.pi * d**2 / 4,
        I=I,
        y_top=d / 2,
        y_bottom=d / 2,
        J=2 * I,
        I_lateral=I,
        name=f"circle d={d:g}",
    )


def tube(d_outer: float, d_inner: float) -> Section:
    """Circular tube: I = pi (do^4 - di^4) / 64, J = 2 I."""
    I = np.pi * (d_outer**4 - d_inner**4) / 64
    return Section(
        A=np.pi * (d_outer**2 - d_inner**2) / 4,
        I=I,
        y_top=d_outer / 2,
        y_bottom=d_outer / 2,
        J=2 * I,
        I_lateral=I,
        name=f"tube {d_outer:g}/{d_inner:g}",
    )


def i_section(b_flange: float, t_flange: float, h: float, t_web: float) -> Section:
    """Doubly symmetric I-section of overall depth ``h`` (exact rectangles; J from thin-walled
    open-section theory, sum b t^3 / 3)."""
    hw = h - 2 * t_flange
    A = 2 * b_flange * t_flange + hw * t_web
    I = (b_flange * h**3 - (b_flange - t_web) * hw**3) / 12
    I_lat = 2 * t_flange * b_flange**3 / 12 + hw * t_web**3 / 12
    J = (2 * b_flange * t_flange**3 + hw * t_web**3) / 3
    return Section(
        A=A, I=I, y_top=h / 2, y_bottom=h / 2, J=J, I_lateral=I_lat, name=f"I {h:g} x {b_flange:g}"
    )


def thin_walled(
    segments: Sequence[tuple[tuple[float, float], tuple[float, float], float]], closed: bool = False
) -> Section:
    """Section built from straight thin-walled segments ``((z1, y1), (z2, y2), t)``.

    Each segment is treated as a thin rectangle of length l and thickness t along its
    mid-line. For an open section J = sum l t^3 / 3. For ``closed=True`` the segments must form
    one closed cell, listed in order, and J follows Bredt-Batho: J = 4 A_enc^2 / sum(l / t).
    """
    segs = [(np.asarray(p, float), np.asarray(q, float), float(t)) for p, q, t in segments]
    lengths = np.array([np.linalg.norm(q - p) for p, q, _ in segs])
    areas = np.array([ln * t for ln, (_, _, t) in zip(lengths, segs, strict=True)])
    mids = np.array([(p + q) / 2 for p, q, _ in segs])
    A = areas.sum()
    zc, yc = (areas[:, None] * mids).sum(axis=0) / A
    Iz = Iy = Izy = 0.0
    for ln, a, (p, q, t), m in zip(lengths, areas, segs, mids, strict=True):
        c, s = (q - p) / ln  # direction cosines of the segment (z, y)
        # rotated thin rectangle about its own centroid, plus parallel-axis terms
        Iz += ln * t / 12 * (ln**2 * s**2 + t**2 * c**2) + a * (m[1] - yc) ** 2
        Iy += ln * t / 12 * (ln**2 * c**2 + t**2 * s**2) + a * (m[0] - zc) ** 2
        Izy += ln * t / 12 * (ln**2 - t**2) * c * s + a * (m[0] - zc) * (m[1] - yc)
    if closed:
        pts = np.array([p for p, _, _ in segs])
        enclosed = 0.5 * abs(
            np.dot(pts[:, 0], np.roll(pts[:, 1], -1)) - np.dot(pts[:, 1], np.roll(pts[:, 0], -1))
        )
        J = 4 * enclosed**2 / sum(ln / t for ln, (_, _, t) in zip(lengths, segs, strict=True))
    else:
        J = sum(ln * t**3 / 3 for ln, (_, _, t) in zip(lengths, segs, strict=True))
    ys = np.concatenate([[p[1], q[1]] for p, q, _ in segs])
    return Section(
        A=A,
        I=Iz,
        y_top=ys.max() - yc,
        y_bottom=yc - ys.min(),
        J=J,
        I_lateral=Iy,
        I_product=Izy,
        name="thin-walled",
    )
