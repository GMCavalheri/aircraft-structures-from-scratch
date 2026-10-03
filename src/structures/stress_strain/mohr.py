"""Mohr's circle for plane stress, and the three circles of a 3D stress state.

Sign convention for plotting: normal stress on the horizontal axis, shear positive *down*
(the clockwise-shear-up convention of Gere & Goodno), so that a rotation of the element by
theta counter-clockwise moves the point by 2 theta counter-clockwise on the circle.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from structures.stress_strain.transform import principal_3d, transform_stress_2d


@dataclass(frozen=True)
class MohrCircle:
    """Mohr's circle of a plane stress state ``[sx, sy, txy]``."""

    sx: float
    sy: float
    txy: float

    @property
    def center(self) -> float:
        return (self.sx + self.sy) / 2

    @property
    def radius(self) -> float:
        return float(np.hypot((self.sx - self.sy) / 2, self.txy))

    @property
    def principal(self) -> tuple[float, float]:
        return self.center + self.radius, self.center - self.radius

    def point(self, theta: float) -> tuple[float, float]:
        """(sigma_x', tau_x'y') on a face rotated by ``theta``; it lies on the circle."""
        sx, _, t = transform_stress_2d([self.sx, self.sy, self.txy], theta)
        return float(sx), float(t)

    def outline(self, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
        phi = np.linspace(0, 2 * np.pi, n)
        return self.center + self.radius * np.cos(phi), self.radius * np.sin(phi)


def mohr_circles_3d(s) -> list[tuple[float, float]]:
    """(center, radius) of the three Mohr circles built on s1 >= s2 >= s3."""
    p, _ = principal_3d(s)
    pairs = [(p[0], p[2]), (p[0], p[1]), (p[1], p[2])]
    return [((a + b) / 2, (a - b) / 2) for a, b in pairs]


def plot_mohr(ax, circle: MohrCircle, color="#2a78d6", scale: float = 1.0, unit: str = "Pa"):
    """Draw a plane-stress Mohr circle with the x and y faces and principal points marked."""
    x, y = circle.outline()
    ax.plot(x / scale, y / scale, color=color, lw=1.5)
    X = (circle.sx / scale, circle.txy / scale)  # shear positive down: plot +txy below axis
    Y = (circle.sy / scale, -circle.txy / scale)
    ax.plot([X[0], Y[0]], [X[1], Y[1]], color=color, lw=1.0, ls="--")
    ax.plot(*X, "o", color=color)
    ax.plot(*Y, "o", color=color, mfc="white")
    ax.annotate("x face", X, textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.annotate("y face", Y, textcoords="offset points", xytext=(6, 6), fontsize=8)
    for p in circle.principal:
        ax.plot(p / scale, 0, "s", color="#0b0b0b", ms=4)
    ax.axhline(0, color="#52514e", lw=0.8)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel(f"Normal stress σ [{unit}]")
    ax.set_ylabel(f"Shear stress τ [{unit}] (positive down)")
