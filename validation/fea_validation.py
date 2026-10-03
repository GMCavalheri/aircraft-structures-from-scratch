"""Phase 3 validation: the finite element solver against closed forms and the Phase 2
Macaulay solver, plus a mesh-convergence study.

Writes ``docs/figures/fea_convergence.png``, ``docs/figures/fea_truss.png``,
``docs/figures/fea_frame.png`` and ``validation/results/fea_solver.md``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from report import rel_err, write_table

from structures import plotting as P
from structures.beam_theory import Beam, DistributedLoad, PointLoad, Support, closed_form
from structures.fea_solver import Model, beam_model, hermite

E, A, I = 70e9, 4e-4, 2e-6
EI = E * I
L, w, Pl = 3.0, 2e3, 5e3


def ss_udl(n, lumped=False):
    m = beam_model(L, n, E, A, I)
    m.fix(0, "xy").fix(n, "y")
    for e in range(n):
        m.distributed(e, -w, lumped=lumped)
    return m, m.solve()


def midspan_interpolated(m, s, n):
    """Midspan deflection by Hermite interpolation inside the element containing L/2."""
    e = int(np.floor(n / 2)) if n % 2 else n // 2 - 1
    xi = (L / 2 - e * L / n) / (L / n)
    el = m.elements[e]
    ul = s.u[m.dofs(el)]
    return float(hermite(xi) @ np.array([ul[1], L / n * ul[2], ul[4], L / n * ul[5]]))


def convergence():
    ref = closed_form.simply_supported_udl(w, L, EI)["deflection"]
    ns = np.array([2, 4, 8, 16, 32, 64])
    nodal_consistent = []
    nodal_lumped = []
    odd = np.array([1, 3, 7, 15, 31, 63])
    interp = []
    for n in ns:
        _, s = ss_udl(n)
        nodal_consistent.append(abs(-s.displacement(n // 2)[1] - ref) / ref)
        _, s = ss_udl(n, lumped=True)
        nodal_lumped.append(abs(-s.displacement(n // 2)[1] - ref) / ref)
    for n in odd:
        m, s = ss_udl(n)
        interp.append(abs(-midspan_interpolated(m, s, n) - ref) / ref)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    floor = 1e-16
    ax.loglog(
        ns,
        np.maximum(nodal_consistent, floor),
        "o-",
        color=P.SERIES[0],
        label="nodal value, consistent loads (exact)",
    )
    ax.loglog(ns, nodal_lumped, "s-", color=P.SERIES[1], label="nodal value, lumped loads")
    ax.loglog(odd, interp, "^-", color=P.SERIES[2], label="between nodes (Hermite interpolation)")
    h = 1 / ns
    ax.loglog(
        ns, nodal_lumped[0] * (h / h[0]) ** 2, ls="--", color=P.REFERENCE, lw=1, label="slope −2"
    )
    ax.loglog(
        odd, interp[1] * (odd / odd[1]) ** -4.0, ls=":", color=P.REFERENCE, lw=1, label="slope −4"
    )
    ax.set_xlabel("Number of elements")
    ax.set_ylabel("Relative error, midspan deflection")
    ax.set_title("Mesh convergence: simply supported beam, uniform load")
    ax.legend(fontsize=8)
    P.save(fig, "fea_convergence.png")
    rate = lambda e: np.log2(e[-2] / e[-1])  # noqa: E731
    rate_odd = np.log(interp[-2] / interp[-1]) / np.log(odd[-1] / odd[-2])
    return [
        [
            "Consistent loads, nodal midspan deflection (n = 2…64)",
            f"max error {max(nodal_consistent):.1e}",
            "exact nodal values",
            "round-off",
        ],
        ["Lumped loads, nodal midspan deflection", f"order {rate(nodal_lumped):.2f}", "2", "✔"],
        ["Hermite interpolation between nodes", f"order {rate_odd:.2f}", "4", "✔"],
    ]


def comparisons():
    rows = []

    def row(name, fe, ref, label):
        rows.append([name, f"{fe:.6g}", f"{ref:.6g}", label, f"{rel_err(fe, ref):+.1e}"])

    m = beam_model(L, 4, E, A, I)
    m.fix(0).load(4, fy=-Pl)
    row(
        "Cantilever tip deflection (4 el.)",
        -m.solve().displacement(4)[1],
        closed_form.cantilever_tip_load(Pl, L, EI)["deflection"],
        "PL³/3EI",
    )
    _, s = ss_udl(6)
    row(
        "Simply supported UDL, midspan (6 el.)",
        -s.displacement(3)[1],
        closed_form.simply_supported_udl(w, L, EI)["deflection"],
        "5wL⁴/384EI",
    )
    m = beam_model(L, 8, E, A, I)
    m.fix(0).fix(8, "y")
    for e in range(8):
        m.distributed(e, -w)
    row(
        "Propped cantilever, prop reaction (8 el.)",
        m.solve().reaction(8)[1],
        closed_form.propped_cantilever_udl(w, L, EI)["R_prop"],
        "3wL/8",
    )
    # overhanging beam with mixed loads vs the Phase 2 Macaulay solver
    mac = (
        Beam(3.0, EI)
        .add(
            Support(0, "fixed"),
            Support(2.2),
            DistributedLoad(0, 3.0, -600.0),
            PointLoad(3.0, -900.0),
        )
        .solve()
    )
    xs = np.array([0, 0.55, 1.1, 1.65, 2.2, 2.6, 3.0])
    m = Model()
    for x in xs:
        m.add_node(x, 0)
    for i in range(len(xs) - 1):
        m.add_frame(i, i + 1, E, A, I)
        m.distributed(i, -600.0)
    m.fix(0).fix(4, "y").load(6, fy=-900.0)
    s = m.solve()
    row(
        "Overhang beam, tip deflection vs Macaulay solver",
        s.displacement(6)[1],
        float(mac.deflection(3.0)),
        "Phase 2",
    )
    row(
        "Overhang beam, moment over roller vs Macaulay",
        s.moment(3, 1.0),
        float(mac.moment(2.2)),
        "Phase 2",
    )
    # L-frame
    h, b = 2.0, 1.2
    m = Model()
    n0, n1, n2 = m.add_node(0, 0), m.add_node(0, h), m.add_node(b, h)
    m.add_frame(n0, n1, E, A, I)
    m.add_frame(n1, n2, E, A, I)
    m.fix(n0).load(n2, fy=-Pl)
    ref = Pl * b**3 / (3 * EI) + Pl * b**2 * h / EI + Pl * h / (E * A)
    row(
        "L-frame tip deflection (2 el.)",
        -m.solve().displacement(n2)[1],
        ref,
        "Pb³/3EI + Pb²h/EI + Ph/EA",
    )
    return rows


def pratt_truss(panels=6, a=2.0, h=2.0, load=10e3):
    m = Model()
    bottom = [m.add_node(i * a, 0) for i in range(panels + 1)]
    top = [m.add_node(i * a, h) for i in range(1, panels)]
    Ab = 1e-3
    for i in range(panels):
        m.add_bar(bottom[i], bottom[i + 1], 200e9, Ab)
    for i in range(panels - 2):
        m.add_bar(top[i], top[i + 1], 200e9, Ab)
    m.add_bar(bottom[0], top[0], 200e9, Ab)
    m.add_bar(top[-1], bottom[-1], 200e9, Ab)
    for i in range(panels - 1):
        m.add_bar(bottom[i + 1], top[i], 200e9, Ab)  # verticals
    for i in range(1, panels - 1):  # Pratt diagonals slope down towards midspan
        if i < panels // 2:
            m.add_bar(top[i - 1], bottom[i + 1], 200e9, Ab)
        else:
            m.add_bar(top[i], bottom[i], 200e9, Ab)
    m.fix(bottom[0], "xy").fix(bottom[-1], "y")
    for n in bottom[1:-1]:
        m.load(n, fy=-load)
    return m, m.solve(), bottom, load


def figure_truss():
    m, s, bottom, load = pratt_truss()
    forces = np.array([s.axial_force(i) for i in range(len(m.elements))]) / 1e3
    scale = 200.0
    segs0 = [np.array([m.nodes[e.n1], m.nodes[e.n2]]) for e in m.elements]
    segs = [s.deflected_shape(i, n=2, scale=scale) for i in range(len(m.elements))]
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.add_collection(LineCollection(segs0, colors=P.TEXT_MUTED, linewidths=0.8, linestyles="--"))
    lim = np.abs(forces).max()
    lc = LineCollection(segs, cmap="coolwarm_r", norm=TwoSlopeNorm(0, -lim, lim), linewidths=2.5)
    lc.set_array(forces)
    ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax, pad=0.01)
    cb.set_label("Axial force [kN] (tension +)")
    ax.autoscale()
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_title(
        f"Pratt truss, {load / 1e3:.0f} kN at each bottom joint (displacements ×{scale:.0f})"
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    P.save(fig, "fea_truss.png")
    mid = len(bottom) // 2
    R = s.reaction(bottom[0])[1]
    return [
        [
            "Pratt truss: support reaction",
            f"{R / 1e3:.3f} kN",
            f"{(len(bottom) - 2) * load / 2 / 1e3:.3f} kN",
            "statics",
            f"{rel_err(R, (len(bottom) - 2) * load / 2):+.1e}",
        ],
        [
            "Pratt truss: midspan bottom-chord force",
            f"{s.axial_force(mid - 1) / 1e3:.3f} kN",
            f"{_pratt_chord(mid, load) / 1e3:.3f} kN",
            "method of sections",
            f"{rel_err(s.axial_force(mid - 1), _pratt_chord(mid, load)):+.1e}",
        ],
    ]


def _pratt_chord(mid, load, a=2.0, h=2.0, panels=6):
    # Section through the bottom chord of the panel just left of midspan (nodes mid-1 to mid).
    # The cut top chord and diagonal meet at the top joint above node mid-1, so moments about
    # it give N_bottom h = R (mid - 1) a - sum of the loads to its left times their lever arms.
    R = (panels - 1) * load / 2
    j = mid - 1
    M = R * j * a - sum(load * (j - i) * a for i in range(1, j))
    return M / h


def figure_frame():
    # Portal frame, fixed bases, lateral and gravity load on the beam.
    h, span = 3.0, 4.0
    m = Model()
    nodes = [m.add_node(0, 0)]
    for y in np.linspace(0, h, 5)[1:]:
        nodes.append(m.add_node(0, y))
    for x in np.linspace(0, span, 7)[1:]:
        nodes.append(m.add_node(x, h))
    for y in np.linspace(h, 0, 5)[1:]:
        nodes.append(m.add_node(span, y))
    for i in range(len(nodes) - 1):
        m.add_frame(nodes[i], nodes[i + 1], 200e9, 5e-3, 8e-5)
    m.fix(nodes[0]).fix(nodes[-1])
    for i in range(4, 10):
        m.distributed(i, -20e3)
    m.load(nodes[4], fx=40e3)
    s = m.solve()
    scale = 50.0
    fig, ax = plt.subplots(figsize=(6, 4.6))
    for i in range(len(m.elements)):
        e = m.elements[i]
        p0 = np.array([m.nodes[e.n1], m.nodes[e.n2]])
        ax.plot(p0[:, 0], p0[:, 1], color=P.TEXT_MUTED, lw=0.8, ls="--")
        d = s.deflected_shape(i, n=21, scale=scale)
        ax.plot(d[:, 0], d[:, 1], color=P.SERIES[0], lw=2)
    ax.set_aspect("equal")
    ax.set_title(f"Portal frame: 40 kN sway + 20 kN/m (displacements ×{scale:.0f})")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    P.save(fig, "fea_frame.png")
    Rx = s.reaction(nodes[0])[0] + s.reaction(nodes[-1])[0]
    Ry = s.reaction(nodes[0])[1] + s.reaction(nodes[-1])[1]
    return [
        [
            "Portal frame: ΣRx, ΣRy",
            f"{Rx / 1e3:.3f}, {Ry / 1e3:.3f} kN",
            "−40.000, 80.000 kN",
            "equilibrium",
            f"{max(abs(rel_err(Rx, -40e3)), abs(rel_err(Ry, 80e3))):.1e}",
        ]
    ]


if __name__ == "__main__":
    P.use_style()
    rows = comparisons() + figure_truss() + figure_frame()
    write_table(
        "fea_solver",
        "FEA solver validation",
        ["Case", "FEA", "Reference", "Source", "Rel. error"],
        rows,
    )
    write_table(
        "fea_convergence",
        "FEA mesh convergence",
        ["Study", "Result", "Expected", ""],
        convergence(),
    )
    figure_frame()
