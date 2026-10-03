"""Classical Lamination Theory (CLT).

Plies are listed in order of increasing z, starting at the face z = -h/2 (Kaw's convention:
ply 1 spans h0 = -h/2 to h1). Then

    A = sum Qbar_k (h_k - h_{k-1}),  B = 1/2 sum Qbar_k (h_k^2 - h_{k-1}^2),
    D = 1/3 sum Qbar_k (h_k^3 - h_{k-1}^3),

    [N; M] = [[A, B], [B, D]] [eps0; kappa],   eps(z) = eps0 + z kappa.

N in N/m, M in N m/m (Kaw, Ch. 4). Ply-by-ply failure discounts failed plies completely
(Qbar = 0), as in Kaw Sec. 5.3.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from structures.composite_laminates import failure
from structures.composite_laminates.lamina import qbar, strain_transformation, stress_transformation
from structures.materials import OrthotropicPly


@dataclass(frozen=True)
class Ply:
    material: OrthotropicPly
    angle: float  # radians
    thickness: float  # m


@dataclass
class PlyPoint:
    """Results at one through-thickness point of one ply."""

    ply: int
    position: str  # "top", "middle" or "bottom" ("top" = the face at larger z)
    z: float
    strain_xy: np.ndarray
    stress_xy: np.ndarray
    strain_12: np.ndarray
    stress_12: np.ndarray


class Laminate:
    def __init__(self, plies: Sequence[Ply]):
        self.plies = list(plies)
        self.active = np.ones(len(self.plies), dtype=bool)  # False once a ply has failed

    @classmethod
    def from_angles(
        cls,
        material: OrthotropicPly,
        angles_deg: Sequence[float],
        ply_thickness: float,
        symmetric: bool = False,
    ) -> Laminate:
        """E.g. ``from_angles(m, [0, 45, -45, 90], t, symmetric=True)`` for [0/45/-45/90]s."""
        angles = list(angles_deg) + (list(angles_deg)[::-1] if symmetric else [])
        return cls([Ply(material, np.radians(a), ply_thickness) for a in angles])

    @property
    def thickness(self) -> float:
        return sum(p.thickness for p in self.plies)

    @property
    def z(self) -> np.ndarray:
        """Ply interface coordinates h_0 = -h/2, ..., h_n = +h/2."""
        return np.concatenate([[0.0], np.cumsum([p.thickness for p in self.plies])]) - (
            self.thickness / 2
        )

    def qbar(self, k: int) -> np.ndarray:
        if not self.active[k]:
            return np.zeros((3, 3))
        p = self.plies[k]
        return qbar(p.material, p.angle)

    def stiffness(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        h = self.z
        A, B, D = np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))
        for k in range(len(self.plies)):
            Q = self.qbar(k)
            A += Q * (h[k + 1] - h[k])
            B += Q * (h[k + 1] ** 2 - h[k] ** 2) / 2
            D += Q * (h[k + 1] ** 3 - h[k] ** 3) / 3
        return A, B, D

    @property
    def A(self):
        return self.stiffness()[0]

    @property
    def B(self):
        return self.stiffness()[1]

    @property
    def D(self):
        return self.stiffness()[2]

    def abd(self) -> np.ndarray:
        A, B, D = self.stiffness()
        return np.block([[A, B], [B, D]])

    def compliance(self) -> np.ndarray:
        return np.linalg.inv(self.abd())

    def engineering_constants(self) -> dict[str, float]:
        """In-plane constants from A* = A^-1 and flexural constants from D* = D^-1 (valid as
        written for symmetric laminates, where B = 0; Kaw Sec. 4.6)."""
        h = self.thickness
        As = np.linalg.inv(self.A)
        Ds = np.linalg.inv(self.D)
        return {
            "Ex": 1 / (h * As[0, 0]),
            "Ey": 1 / (h * As[1, 1]),
            "Gxy": 1 / (h * As[2, 2]),
            "nuxy": -As[0, 1] / As[0, 0],
            "nuyx": -As[0, 1] / As[1, 1],
            "Ex_flex": 12 / (h**3 * Ds[0, 0]),
            "Ey_flex": 12 / (h**3 * Ds[1, 1]),
            "Gxy_flex": 12 / (h**3 * Ds[2, 2]),
            "nuxy_flex": -Ds[0, 1] / Ds[0, 0],
            "nuyx_flex": -Ds[0, 1] / Ds[1, 1],
        }

    def response(self, N=(0, 0, 0), M=(0, 0, 0)) -> tuple[np.ndarray, np.ndarray]:
        """Mid-plane strains eps0 and curvatures kappa under resultants N (N/m), M (N m/m)."""
        x = np.linalg.solve(self.abd(), np.concatenate([N, M]).astype(float))
        return x[:3], x[3:]

    def ply_points(self, N=(0, 0, 0), M=(0, 0, 0)) -> list[PlyPoint]:
        eps0, kap = self.response(N, M)
        h = self.z
        out = []
        for k, p in enumerate(self.plies):
            for pos, z in (("bottom", h[k]), ("middle", (h[k] + h[k + 1]) / 2), ("top", h[k + 1])):
                exy = eps0 + z * kap
                sxy = self.qbar(k) @ exy
                out.append(
                    PlyPoint(
                        k,
                        pos,
                        float(z),
                        exy,
                        sxy,
                        strain_transformation(p.angle) @ exy,
                        stress_transformation(p.angle) @ sxy,
                    )
                )
        return out

    def strength_ratios(
        self, N=(0, 0, 0), M=(0, 0, 0), criterion: str = "tsai_wu", f12_star: float = -0.5
    ) -> list[tuple[PlyPoint, float]]:
        """Strength ratio at every ply point of every active ply."""
        out = []
        for pt in self.ply_points(N, M):
            if not self.active[pt.ply]:
                continue
            mat = self.plies[pt.ply].material
            if criterion == "tsai_wu":
                sr = failure.tsai_wu(pt.stress_12, mat, f12_star)
            elif criterion == "tsai_hill":
                sr = failure.tsai_hill(pt.stress_12, mat)
            elif criterion == "modified_tsai_hill":
                sr = failure.tsai_hill(pt.stress_12, mat, modified=True)
            elif criterion == "max_stress":
                sr = failure.max_stress(pt.stress_12, mat)[0]
            elif criterion == "max_strain":
                sr = failure.max_strain(pt.strain_12, mat)[0]
            else:
                raise ValueError(f"criterion must be one of {failure.CRITERIA}")
            out.append((pt, sr))
        return out

    def first_ply_failure(
        self, N=(0, 0, 0), M=(0, 0, 0), criterion: str = "tsai_wu", f12_star: float = -0.5
    ) -> tuple[float, list[int]]:
        """(load factor, plies failing at it) for the load (N, M) scaled proportionally."""
        srs = self.strength_ratios(N, M, criterion, f12_star)
        sr_min = min(sr for _, sr in srs)
        failed = sorted({pt.ply for pt, sr in srs if sr <= sr_min * (1 + 1e-9)})
        return sr_min, failed

    def ply_by_ply_failure(
        self, N=(0, 0, 0), M=(0, 0, 0), criterion: str = "tsai_wu", f12_star: float = -0.5
    ) -> list[tuple[float, list[int]]]:
        """Progressive failure with complete ply discount. Returns successive (load factor,
        newly failed plies); the load factor applies to the proportional load (N, M). Plies
        already over-stressed when another ply fails are failed at that same load."""
        saved = self.active.copy()
        history = []
        try:
            while self.active.any():
                try:
                    sr, failed = self.first_ply_failure(N, M, criterion, f12_star)
                except np.linalg.LinAlgError:
                    break  # remaining laminate cannot carry load
                if history and sr < history[-1][0]:
                    sr = history[-1][0]  # damage cascades at the previous load level
                self.active[failed] = False
                history.append((float(sr), failed))
        finally:
            self.active = saved
        return history
