"""A small 2D finite element model: nodes, bar/frame elements, supports and loads.

Workflow: build a :class:`Model`, call :meth:`Model.solve`, read displacements, reactions and
element end forces from the returned :class:`Solution`. The global stiffness matrix is
assembled in sparse COO form, converted to CSR, partitioned into free and restrained DOFs
and solved with ``scipy.sparse.linalg.spsolve``:

    K_ff u_f = F_f - K_fr u_r,     R = K_rf u_f + K_rr u_r - F_r.

Rotational DOFs of nodes that only touch bar elements carry no stiffness; they are
restrained automatically so pure trusses solve without user bookkeeping.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from structures.fea_solver.elements import (
    bar_stiffness,
    consistent_load,
    frame_geometric_stiffness,
    frame_stiffness,
    hermite,
    rotation,
)

DOF = {"x": 0, "y": 1, "r": 2}


@dataclass
class Element:
    n1: int
    n2: int
    E: float
    A: float
    I: float = 0.0  # zero for a bar
    kind: str = "frame"  # "bar" or "frame"
    q: np.ndarray = field(default_factory=lambda: np.zeros(4))  # qy1, qy2, qx1, qx2 (local)
    lumped: bool = False  # apply q as simple nodal forces (no moments) instead of consistently


class Model:
    def __init__(self):
        self.nodes: list[tuple[float, float]] = []
        self.elements: list[Element] = []
        self.restraints: dict[int, float] = {}  # global DOF -> prescribed displacement
        self.nodal_loads: dict[int, float] = {}  # global DOF -> force

    # --- building -----------------------------------------------------------------------
    def add_node(self, x: float, y: float) -> int:
        self.nodes.append((float(x), float(y)))
        return len(self.nodes) - 1

    def add_bar(self, n1: int, n2: int, E: float, A: float) -> int:
        self.elements.append(Element(n1, n2, E, A, 0.0, "bar"))
        return len(self.elements) - 1

    def add_frame(self, n1: int, n2: int, E: float, A: float, I: float) -> int:
        self.elements.append(Element(n1, n2, E, A, I, "frame"))
        return len(self.elements) - 1

    def fix(self, node: int, dofs: str = "xyr", value: float = 0.0) -> Model:
        """Restrain ``dofs`` (any of "x", "y", "r") of ``node`` to ``value`` (default 0)."""
        for d in dofs:
            self.restraints[3 * node + DOF[d]] = float(value)
        return self

    def load(self, node: int, fx: float = 0.0, fy: float = 0.0, mz: float = 0.0) -> Model:
        for d, v in zip("xyr", (fx, fy, mz), strict=True):
            if v:
                i = 3 * node + DOF[d]
                self.nodal_loads[i] = self.nodal_loads.get(i, 0.0) + float(v)
        return self

    def distributed(
        self,
        element: int,
        qy1: float,
        qy2: float | None = None,
        qx1: float = 0.0,
        qx2: float | None = None,
        lumped: bool = False,
    ) -> Model:
        """Linearly varying load on a frame element in *local* axes (qy transverse, qx axial).

        By default it is converted to work-equivalent (consistent) nodal forces and moments.
        ``lumped=True`` splits the resultant between the two nodes as forces only, the cruder
        scheme whose error is shown in the convergence study.
        """
        e = self.elements[element]
        e.q += [qy1, qy1 if qy2 is None else qy2, qx1, qx1 if qx2 is None else qx2]
        e.lumped = lumped
        return self

    # --- geometry helpers -----------------------------------------------------------------
    def geometry(self, e: Element) -> tuple[float, float, float]:
        (x1, y1), (x2, y2) = self.nodes[e.n1], self.nodes[e.n2]
        L = float(np.hypot(x2 - x1, y2 - y1))
        return L, (x2 - x1) / L, (y2 - y1) / L

    @staticmethod
    def dofs(e: Element) -> np.ndarray:
        return np.r_[3 * e.n1 : 3 * e.n1 + 3, 3 * e.n2 : 3 * e.n2 + 3]

    def local_stiffness(self, e: Element) -> np.ndarray:
        L, _, _ = self.geometry(e)
        if e.kind == "bar":
            return bar_stiffness(e.E, e.A, L)
        return frame_stiffness(e.E, e.A, e.I, L)

    def local_load(self, e: Element) -> np.ndarray:
        L, _, _ = self.geometry(e)
        if e.lumped:
            qy, qx = (e.q[0] + e.q[1]) / 2 * L / 2, (e.q[2] + e.q[3]) / 2 * L / 2
            return np.array([qx, qy, 0.0, qx, qy, 0.0])
        return consistent_load(L, *e.q)

    @property
    def n_dof(self) -> int:
        return 3 * len(self.nodes)

    # --- assembly -------------------------------------------------------------------------
    def _assemble(self, element_matrix) -> sp.csr_matrix:
        rows, cols, vals = [], [], []
        for i, e in enumerate(self.elements):
            _, c, s = self.geometry(e)
            T = rotation(c, s)
            k = T.T @ element_matrix(i, e) @ T
            d = self.dofs(e)
            rows.append(np.repeat(d, 6))
            cols.append(np.tile(d, 6))
            vals.append(k.ravel())
        n = self.n_dof
        return sp.coo_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
        ).tocsr()

    def stiffness(self) -> sp.csr_matrix:
        return self._assemble(lambda i, e: self.local_stiffness(e))

    def geometric_stiffness(self, axial_forces) -> sp.csr_matrix:
        """Global K_G for element axial forces ``axial_forces[i]`` (tension positive)."""
        forces = np.asarray(axial_forces, float)

        def kg(i, e):
            L, _, _ = self.geometry(e)
            return frame_geometric_stiffness(forces[i], L)

        return self._assemble(kg)

    def load_vector(self) -> np.ndarray:
        F = np.zeros(self.n_dof)
        for i, v in self.nodal_loads.items():
            F[i] += v
        for e in self.elements:
            if np.any(e.q):
                _, c, s = self.geometry(e)
                F[self.dofs(e)] += rotation(c, s).T @ self.local_load(e)
        return F

    def restrained_dofs(self) -> dict[int, float]:
        """User restraints plus rotations of nodes that are not attached to a frame element."""
        r = dict(self.restraints)
        framed = {n for e in self.elements if e.kind == "frame" for n in (e.n1, e.n2)}
        for node in range(len(self.nodes)):
            if node not in framed:
                r.setdefault(3 * node + 2, 0.0)
        return r

    # --- solve ----------------------------------------------------------------------------
    def solve(self) -> Solution:
        K = self.stiffness()
        F = self.load_vector()
        restrained = self.restrained_dofs()
        r = np.array(sorted(restrained), dtype=int)
        f = np.setdiff1d(np.arange(self.n_dof), r)
        u = np.zeros(self.n_dof)
        u[r] = [restrained[i] for i in r]
        Kff = K[f][:, f]
        rhs = F[f] - K[f][:, r] @ u[r]
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            try:
                u[f] = spsolve(Kff.tocsc(), rhs)
            except MatrixRankWarning:
                u[f] = np.nan
        if not np.all(np.isfinite(u)):
            raise np.linalg.LinAlgError("singular stiffness matrix: the model is a mechanism")
        R = np.zeros(self.n_dof)
        R[r] = K[r] @ u - F[r]
        return Solution(self, u, R)


@dataclass
class Solution:
    model: Model
    u: np.ndarray  # global displacements, 3 per node
    R: np.ndarray  # reactions at restrained DOFs (zero elsewhere)

    def displacement(self, node: int) -> np.ndarray:
        """(u, v, theta) of a node."""
        return self.u[3 * node : 3 * node + 3]

    def reaction(self, node: int) -> np.ndarray:
        return self.R[3 * node : 3 * node + 3]

    def end_forces(self, element: int) -> np.ndarray:
        """Local end forces on the element [N1, V1, M1, N2, V2, M2] (forces the nodes exert on
        the element): f = k T u_e - f_eq."""
        m = self.model
        e = m.elements[element]
        _, c, s = m.geometry(e)
        return m.local_stiffness(e) @ rotation(c, s) @ self.u[m.dofs(e)] - m.local_load(e)

    def axial_force(self, element: int) -> float:
        """Axial force, tension positive."""
        return float(self.end_forces(element)[3])

    def axial_stress(self, element: int) -> float:
        return self.axial_force(element) / self.model.elements[element].A

    def moment(self, element: int, xi) -> np.ndarray:
        """Internal (sagging positive) bending moment at xi in [0, 1] along a frame element,
        exact for linearly varying loads: M(s) = -M1 + V1 s + integral of q (s - t) dt."""
        m = self.model
        e = m.elements[element]
        L, _, _ = m.geometry(e)
        f = self.end_forces(element)
        s = np.asarray(xi, float) * L
        q1, q2 = e.q[0], e.q[1]
        load = q1 * s**2 / 2 + (q2 - q1) / L * s**3 / 6
        return -f[2] + f[1] * s + load

    def deflected_shape(self, element: int, n: int = 21, scale: float = 1.0) -> np.ndarray:
        """Global coordinates (n x 2) of the deformed element: linear axial and cubic
        Hermite transverse interpolation of the nodal displacements."""
        m = self.model
        e = m.elements[element]
        L, c, s = m.geometry(e)
        ul = rotation(c, s) @ self.u[m.dofs(e)]
        xi = np.linspace(0, 1, n)
        ax = (1 - xi) * ul[0] + xi * ul[3]
        if e.kind == "frame":
            tr = hermite(xi).T @ np.array([ul[1], L * ul[2], ul[4], L * ul[5]])
        else:
            tr = (1 - xi) * ul[1] + xi * ul[4]
        x0 = np.array(m.nodes[e.n1])
        local = np.column_stack([xi * L + scale * ax, scale * tr])
        R = np.array([[c, -s], [s, c]])
        return x0 + local @ R.T
