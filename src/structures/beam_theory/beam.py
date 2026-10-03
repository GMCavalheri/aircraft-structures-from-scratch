"""Euler-Bernoulli beam solver using Macaulay (singularity) functions.

Conventions: x from the left end, loads and deflection v positive **up**, couples positive
**counter-clockwise**, sagging bending moment positive, so that

    dV/dx = q,   dM/dx = V,   EI v'' = M.

The internal forces at a section are built from everything to its left as sums of
<x - a>^n terms. Unknown support reactions are point forces (and couples, at fixed supports)
with unknown magnitudes. Together with the two integration constants of EI v they are found
from one linear system:

    v = 0 at every support,  v' = 0 at every fixed support,  V(L+) = 0,  M(L+) = 0.

The last two are overall force and moment equilibrium, so one routine handles statically
determinate and indeterminate beams alike (Gere & Goodno, *Mechanics of Materials*, Secs.
9.5 and 10.4).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import factorial

import numpy as np


@dataclass(frozen=True)
class PointLoad:
    x: float
    P: float  # N, positive up


@dataclass(frozen=True)
class Couple:
    x: float
    M: float  # N m, positive counter-clockwise


@dataclass(frozen=True)
class DistributedLoad:
    """Linearly varying load from ``x1`` (intensity ``w1``) to ``x2`` (``w2``), N/m, up +."""

    x1: float
    x2: float
    w1: float
    w2: float | None = None

    @property
    def w_end(self) -> float:
        return self.w1 if self.w2 is None else self.w2


SUPPORT_KINDS = ("pin", "roller", "fixed")


@dataclass(frozen=True)
class Support:
    x: float
    kind: str = "pin"  # "pin"/"roller" restrain v; "fixed" restrains v and v'

    def __post_init__(self):
        if self.kind not in SUPPORT_KINDS:
            raise ValueError(f"support kind must be one of {SUPPORT_KINDS}")


def _mac(x, a, n):
    """<x - a>^n for n >= 0, with <0>^0 = 1 (a load at x = a is included at x = a)."""
    x = np.asarray(x, float)
    return np.where(x >= a, (x - a) ** n if n > 0 else 1.0, 0.0)


@dataclass
class _Terms:
    """Moment expressed as sum c <x - a>^k (k >= 0); shear as sum c <x - a>^k."""

    M: list[tuple[float, float, int]] = field(default_factory=list)
    V: list[tuple[float, float, int]] = field(default_factory=list)

    def eval(self, which: str, x, integrations: int = 0):
        out = np.zeros_like(np.asarray(x, float))
        for c, a, k in getattr(self, which):
            out = out + c * factorial(k) / factorial(k + integrations) * _mac(
                x, a, k + integrations
            )
        return out


def _load_terms(load) -> _Terms:
    t = _Terms()
    if isinstance(load, PointLoad):
        t.V.append((load.P, load.x, 0))
        t.M.append((load.P, load.x, 1))
    elif isinstance(load, Couple):
        t.M.append((-load.M, load.x, 0))  # a CCW couple reduces the sagging moment to its right
    elif isinstance(load, DistributedLoad):
        a, b, w1, w2 = load.x1, load.x2, load.w1, load.w_end
        slope = (w2 - w1) / (b - a)
        # q = w1<x-a>^0 + slope<x-a>^1 - w2<x-b>^0 - slope<x-b>^1
        for c, x0, n in [(w1, a, 0), (slope, a, 1), (-w2, b, 0), (-slope, b, 1)]:
            if c != 0:
                t.V.append((c / (n + 1), x0, n + 1))
                t.M.append((c / ((n + 1) * (n + 2)), x0, n + 2))
    else:
        raise TypeError(f"unknown load type {type(load).__name__}")
    return t


@dataclass
class BeamSolution:
    length: float
    EI: float
    reactions: list[tuple[Support, float, float]]  # (support, force, couple)
    _terms: _Terms
    _C1: float
    _C2: float

    def _inside(self, x):
        # Values at x = a are right-hand limits (a load at a is included), except at the right
        # end, where the left-hand limit is the meaningful one (a reaction there zeroes V, M).
        return np.minimum(np.asarray(x, float), self.length * (1 - 1e-12))

    def shear(self, x):
        return self._terms.eval("V", self._inside(x))

    def moment(self, x):
        return self._terms.eval("M", self._inside(x))

    def slope(self, x):
        return (self._terms.eval("M", x, 1) + self._C1) / self.EI

    def deflection(self, x):
        return (self._terms.eval("M", x, 2) + self._C1 * np.asarray(x, float) + self._C2) / self.EI

    def bending_stress(self, x, y, I: float):
        """sigma_x = -M y / I at height ``y`` above the neutral axis (tension positive)."""
        return -self.moment(x) * y / I

    def extreme(self, quantity: str, n: int = 2001) -> tuple[float, float]:
        """(x, value) of the largest |quantity| for ``"shear"``, ``"moment"``, ``"slope"`` or
        ``"deflection"`` on a dense grid."""
        x = np.linspace(0, self.length, n)
        f = getattr(self, quantity)(x)
        i = int(np.argmax(np.abs(f)))
        return float(x[i]), float(f[i])


class Beam:
    """Prismatic Euler-Bernoulli beam of length ``L`` and flexural rigidity ``EI``."""

    def __init__(self, length: float, EI: float):
        self.length = float(length)
        self.EI = float(EI)
        self.loads: list = []
        self.supports: list[Support] = []

    def add(self, *items) -> Beam:
        for item in items:
            (self.supports if isinstance(item, Support) else self.loads).append(item)
        return self

    def solve(self) -> BeamSolution:
        L = self.length
        known = _Terms()
        for load in self.loads:
            t = _load_terms(load)
            known.M += t.M
            known.V += t.V
        # unknown columns: one force per support, one couple per fixed support, C1, C2
        unknowns: list[tuple[str, Support | None]] = [("F", s) for s in self.supports]
        unknowns += [("M", s) for s in self.supports if s.kind == "fixed"]
        unknowns += [("C1", None), ("C2", None)]
        n = len(unknowns)

        def unit_terms(kind, s):
            if kind == "F":
                return _load_terms(PointLoad(s.x, 1.0))
            return _load_terms(Couple(s.x, 1.0))

        rows, rhs = [], []

        def condition(fn):
            # fn(terms, C1, C2) -> value; linear in the unknowns
            row = np.zeros(n)
            for j, (kind, s) in enumerate(unknowns):
                if kind in ("F", "M"):
                    row[j] = fn(unit_terms(kind, s), 0.0, 0.0)
                else:
                    row[j] = fn(
                        _Terms(), 1.0 if kind == "C1" else 0.0, 1.0 if kind == "C2" else 0.0
                    )
            rows.append(row)
            rhs.append(-fn(known, 0.0, 0.0))

        for s in self.supports:
            condition(lambda t, c1, c2, s=s: float(t.eval("M", s.x, 2) + c1 * s.x + c2))
            if s.kind == "fixed":
                condition(lambda t, c1, c2, s=s: float(t.eval("M", s.x, 1) + c1))
        condition(lambda t, c1, c2: float(t.eval("V", L)))
        condition(lambda t, c1, c2: float(t.eval("M", L)))

        A, b = np.array(rows), np.array(rhs)
        if A.shape[0] != n or np.linalg.matrix_rank(A) < n:
            raise ValueError("beam is a mechanism (insufficient or ill-placed supports)")
        sol = np.linalg.solve(A, b)

        terms = _Terms(M=list(known.M), V=list(known.V))
        forces = {id(s): 0.0 for s in self.supports}
        couples = {id(s): 0.0 for s in self.supports}
        for (kind, s), value in zip(unknowns, sol, strict=True):
            if kind in ("F", "M"):
                ut = unit_terms(kind, s)
                terms.M += [(c * value, a, k) for c, a, k in ut.M]
                terms.V += [(c * value, a, k) for c, a, k in ut.V]
                (forces if kind == "F" else couples)[id(s)] = float(value)
        C1, C2 = sol[-2], sol[-1]
        reactions = [(s, forces[id(s)], couples[id(s)]) for s in self.supports]
        return BeamSolution(L, self.EI, reactions, terms, float(C1), float(C2))
