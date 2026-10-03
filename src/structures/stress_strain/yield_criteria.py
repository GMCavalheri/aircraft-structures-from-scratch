"""Yield criteria for ductile isotropic metals.

von Mises (distortion energy):
    s_vm = sqrt(((s1 - s2)^2 + (s2 - s3)^2 + (s3 - s1)^2) / 2)
         = sqrt(sx^2 + sy^2 + sz^2 - sx sy - sy sz - sz sx + 3 (txy^2 + tyz^2 + txz^2))
Tresca (maximum shear): s_tresca = s1 - s3, using all three principal stresses (the
out-of-plane zero counts for plane stress).
"""

from __future__ import annotations

import numpy as np

from structures.stress_strain.transform import as_tensor, principal_3d


def von_mises(s) -> float:
    """Von Mises equivalent stress of a 3x3 tensor or a 2D/3D Voigt vector."""
    t = as_tensor(s)
    sx, sy, sz = t[0, 0], t[1, 1], t[2, 2]
    shear = t[0, 1] ** 2 + t[1, 2] ** 2 + t[0, 2] ** 2
    return float(np.sqrt(sx**2 + sy**2 + sz**2 - sx * sy - sy * sz - sz * sx + 3 * shear))


def tresca(s) -> float:
    """Tresca equivalent stress s1 - s3 (twice the absolute maximum shear stress)."""
    p, _ = principal_3d(s)
    return float(p[0] - p[2])


def safety_factor(s, yield_strength: float, criterion: str = "von_mises") -> float:
    """Factor of safety against yield, ``Fy / s_eq``, for ``"von_mises"`` or ``"tresca"``."""
    eq = {"von_mises": von_mises, "tresca": tresca}[criterion](s)
    return float("inf") if eq == 0 else yield_strength / eq
