"""Isotropic (metallic) and orthotropic (ply) material records, in SI units.

The numbers live in ``data/*.csv`` with their citations; this module only converts units.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"

KSI = 6.894757e6  # Pa per ksi (1000 lbf/in^2)
MSI = 1e3 * KSI  # Pa per 10^6 psi
LB_PER_IN3 = 27679.90  # kg/m^3 per lb/in^3

_UNITS = {"GPa": 1e9, "MPa": 1e6, "Msi": MSI, "ksi": KSI}


def _rows(name: str) -> list[dict[str, str]]:
    with open(DATA / name) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


@dataclass(frozen=True)
class IsotropicMaterial:
    """Linear-elastic isotropic material with design allowables (all SI: Pa, kg/m^3).

    ``Ftu``/``Fty``/``Fcy``/``Fsu`` are ultimate tensile, tensile yield, compressive yield and
    ultimate shear allowables; ``basis`` is the statistical basis (A, B or S).
    """

    name: str
    E: float
    nu: float
    G: float | None = None
    Ec: float | None = None
    density: float | None = None
    Ftu: float | None = None
    Fty: float | None = None
    Fcy: float | None = None
    Fsu: float | None = None
    basis: str = ""
    source: str = ""

    @property
    def shear_modulus(self) -> float:
        """G as tabulated, else the isotropic relation E / (2 (1 + nu))."""
        return self.G if self.G is not None else self.E / (2 * (1 + self.nu))


@dataclass(frozen=True)
class OrthotropicPly:
    """Unidirectional ply: in-plane engineering constants (Pa) and strengths (Pa, positive).

    ``Xt``/``Xc`` fibre-direction tension/compression strength, ``Yt``/``Yc`` transverse, ``S``
    in-plane shear strength.
    """

    name: str
    E1: float
    E2: float
    G12: float
    nu12: float
    Xt: float = float("inf")
    Xc: float = float("inf")
    Yt: float = float("inf")
    Yc: float = float("inf")
    S: float = float("inf")
    source: str = ""

    @property
    def nu21(self) -> float:
        return self.nu12 * self.E2 / self.E1


@cache
def isotropic(key: str, basis: str = "B") -> IsotropicMaterial:
    """MIL-HDBK-5J allowables for ``key`` (``"2024-T3"``, ``"7075-T6"``, ``"Ti-6Al-4V"``,
    ``"4340"``) on the requested statistical basis (falls back to S-basis rows)."""
    rows = [r for r in _rows("mil_hdbk_5j.csv") if r["key"] == key]
    if not rows:
        raise KeyError(key)
    row = next((r for r in rows if r["basis"] == basis), None)
    if row is None:
        row = next(r for r in rows if r["basis"] == "S")
    return IsotropicMaterial(
        name=row["description"],
        E=float(row["E_msi"]) * MSI,
        Ec=float(row["Ec_msi"]) * MSI,
        G=float(row["G_msi"]) * MSI,
        nu=float(row["nu"]),
        density=float(row["density_lb_in3"]) * LB_PER_IN3,
        Ftu=float(row["Ftu_ksi"]) * KSI,
        Fty=float(row["Fty_ksi"]) * KSI,
        Fcy=float(row["Fcy_ksi"]) * KSI,
        Fsu=float(row["Fsu_ksi"]) * KSI,
        basis=row["basis"],
        source=row["source"],
    )


@cache
def ply(key: str) -> OrthotropicPly:
    """Unidirectional ply properties for ``key`` (``"T300/5208"``, ``"T300/976"``)."""
    for row in _rows("plies.csv"):
        if row["key"] == key:
            m, s = _UNITS[row["modulus_unit"]], _UNITS[row["strength_unit"]]
            return OrthotropicPly(
                name=key,
                E1=float(row["E1"]) * m,
                E2=float(row["E2"]) * m,
                G12=float(row["G12"]) * m,
                nu12=float(row["nu12"]),
                Xt=float(row["Xt"]) * s,
                Xc=float(row["Xc"]) * s,
                Yt=float(row["Yt"]) * s,
                Yc=float(row["Yc"]) * s,
                S=float(row["S"]) * s,
                source=row["source"],
            )
    raise KeyError(key)
