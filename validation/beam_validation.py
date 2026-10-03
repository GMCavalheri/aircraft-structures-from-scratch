"""Phase 2 validation: Macaulay beam solver and torsion against closed-form results.

Writes ``docs/figures/beam_diagrams.png`` and ``validation/results/beam_theory.md``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from report import rel_err, write_table

from structures import plotting as P
from structures.beam_theory import (
    Beam,
    Couple,
    DistributedLoad,
    PointLoad,
    Support,
    bredt_batho,
    circular_shaft,
    closed_form,
    open_thin_walled,
)

L, EI, Pm, w = 2.0, 3.0e5, 1.2e3, 800.0


def cases():
    """(name, solved beam, quantity extractor, reference value, reference label)."""
    out = []

    def add(name, beam, get, ref, label):
        s = beam.solve()
        out.append([name, get(s), ref, label])

    add(
        "Cantilever, tip load — tip deflection",
        Beam(L, EI).add(Support(0, "fixed"), PointLoad(L, -Pm)),
        lambda s: -s.deflection(L),
        closed_form.cantilever_tip_load(Pm, L, EI)["deflection"],
        "PL³/3EI",
    )
    add(
        "Cantilever, UDL — tip deflection",
        Beam(L, EI).add(Support(0, "fixed"), DistributedLoad(0, L, -w)),
        lambda s: -s.deflection(L),
        closed_form.cantilever_udl(w, L, EI)["deflection"],
        "wL⁴/8EI",
    )
    add(
        "Simply supported, centre load — midspan deflection",
        Beam(L, EI).add(Support(0), Support(L), PointLoad(L / 2, -Pm)),
        lambda s: -s.deflection(L / 2),
        closed_form.simply_supported_center_load(Pm, L, EI)["deflection"],
        "PL³/48EI",
    )
    add(
        "Simply supported, UDL — midspan deflection",
        Beam(L, EI).add(Support(0), Support(L), DistributedLoad(0, L, -w)),
        lambda s: -s.deflection(L / 2),
        closed_form.simply_supported_udl(w, L, EI)["deflection"],
        "5wL⁴/384EI",
    )
    add(
        "Simply supported, triangular load — max deflection",
        Beam(L, EI).add(Support(0), Support(L), DistributedLoad(0, L, 0, -w)),
        lambda s: -s.extreme("deflection", 20001)[1],
        closed_form.simply_supported_triangular(w, L, EI)["deflection"],
        "0.006522 w₀L⁴/EI",
    )
    add(
        "Simply supported, end couple — far-end slope",
        Beam(L, EI).add(Support(0), Support(L), Couple(L, 500.0)),
        lambda s: abs(s.slope(0)),
        closed_form.simply_supported_end_couple(500.0, L, EI)["slope_far"],
        "M₀L/6EI",
    )
    add(
        "Propped cantilever, UDL — prop reaction (indeterminate)",
        Beam(L, EI).add(Support(0, "fixed"), Support(L), DistributedLoad(0, L, -w)),
        lambda s: s.reactions[1][1],
        closed_form.propped_cantilever_udl(w, L, EI)["R_prop"],
        "3wL/8",
    )
    add(
        "Fixed-fixed, UDL — end moment (indeterminate)",
        Beam(L, EI).add(Support(0, "fixed"), Support(L, "fixed"), DistributedLoad(0, L, -w)),
        lambda s: -s.moment(0),
        closed_form.fixed_fixed_udl(w, L, EI)["moment"],
        "wL²/12",
    )
    add(
        "Two-span continuous, UDL — middle reaction (indeterminate)",
        Beam(2 * L, EI).add(Support(0), Support(L), Support(2 * L), DistributedLoad(0, 2 * L, -w)),
        lambda s: s.reactions[1][1],
        closed_form.two_span_continuous_udl(w, L, EI)["R_mid"],
        "10wL/8",
    )
    rows = [[n, f"{c:.6g}", f"{r:.6g}", lab, f"{rel_err(c, r):+.1e}"] for n, c, r, lab in out]
    return rows


def torsion_rows():
    r, t, G, T = 0.05, 0.001, 27e9, 400.0
    exact = circular_shaft(T, G, 2 * r + t, 2 * r - t)
    bb = bredt_batho(T, G, np.pi * r**2, [(2 * np.pi * r, t)])
    op = open_thin_walled(T, G, [(2 * np.pi * r, t)])
    return [
        [
            "Thin tube r = 50 mm, t = 1 mm: J, Bredt–Batho vs exact annulus",
            f"{bb.J:.4e} m⁴",
            f"{exact.J:.4e} m⁴",
            "π(dₒ⁴−dᵢ⁴)/32",
            f"{rel_err(bb.J, exact.J):+.1e}",
        ],
        [
            "Same tube slit open: stiffness ratio closed/open",
            f"{bb.J / op.J:.0f}",
            f"{3 * r**2 / t**2:.0f}",
            "3r²/t²",
            f"{rel_err(bb.J / op.J, 3 * r**2 / t**2):+.1e}",
        ],
    ]


def figure():
    beam = Beam(3.0, EI).add(
        Support(0, "fixed"), Support(2.2), DistributedLoad(0, 3.0, -600.0), PointLoad(3.0, -900.0)
    )
    s = beam.solve()
    x = np.linspace(0, 3.0, 1201)
    fig, axes = plt.subplots(3, 1, figsize=(7, 6.6), sharex=True)
    for ax, f, lab, k in [
        (axes[0], s.shear(x) / 1e3, "Shear V [kN]", 0),
        (axes[1], s.moment(x) / 1e3, "Moment M [kN·m]", 1),
        (axes[2], s.deflection(x) * 1e3, "Deflection v [mm]", 2),
    ]:
        ax.plot(x, f, color=P.SERIES[k])
        ax.fill_between(x, f, color=P.SERIES[k], alpha=0.12, lw=0)
        ax.axhline(0, color=P.TEXT_MUTED, lw=0.8)
        ax.set_ylabel(lab)
    for xs in (0, 2.2):
        for ax in axes:
            ax.axvline(xs, color=P.GRID, lw=1.0, zorder=0)
    axes[0].set_title("Fixed at x = 0, roller at 2.2 m, 600 N/m UDL + 900 N tip load")
    axes[2].set_xlabel("x [m]")
    fig.tight_layout()
    P.save(fig, "beam_diagrams.png")


if __name__ == "__main__":
    P.use_style()
    header = ["Case", "Computed", "Closed form", "Formula", "Rel. error"]
    write_table(
        "beam_theory",
        "Beam theory validation",
        header,
        cases() + torsion_rows(),
        notes="EI = 3×10⁵ N·m², L = 2 m, P = 1.2 kN, w = 800 N/m; SI units.",
    )
    figure()
