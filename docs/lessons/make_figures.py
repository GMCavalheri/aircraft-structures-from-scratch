"""Teaching figures for the lesson series (docs/figures/lessons/).

Run with ``uv run python docs/lessons/make_figures.py``. Validation figures live in
``validation/``; these illustrate ideas.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from structures import plotting as P
from structures.buckling import column_model, euler_load, linear_buckling
from structures.composite_laminates import failure, qbar, stress_transformation
from structures.fatigue import count_cycles, handbook_curve
from structures.fea_solver import Model, beam_model, hermite
from structures.materials import KSI, isotropic, ply
from structures.stress_strain import principal_2d, transform_stress_2d

OUT = P.FIGURES / "lessons"
DEG = np.pi / 180
MPA = 1e6


def save(fig, name):
    return P.save(fig, name, OUT)


def l1_transformation():
    s = np.array([-20, 90, 60]) * MPA
    th = np.linspace(0, 180, 361) * DEG
    sx, sy, txy = transform_stress_2d(s, th)
    s1, s2, tp = principal_2d(s)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(th / DEG, sx / MPA, color=P.SERIES[0], label="σx′")
    ax.plot(th / DEG, sy / MPA, color=P.SERIES[1], label="σy′")
    ax.plot(th / DEG, txy / MPA, color=P.SERIES[2], label="τx′y′")
    for t, lab in [(tp, "principal"), (tp - 45 * DEG, "max shear")]:
        t = t % np.pi
        ax.axvline(t / DEG, color=P.GRID, lw=1.2)
        ax.text(t / DEG + 1.5, -75, lab, fontsize=8, color=P.TEXT_MUTED)
    ax.axhline(0, color=P.TEXT_MUTED, lw=0.8)
    ax.set_xlabel("Rotation of the element θ [deg]")
    ax.set_ylabel("Stress [MPa]")
    ax.set_title("Rotating the element: σx = −20, σy = 90, τxy = 60 MPa")
    ax.legend(loc="upper right")
    return save(fig, "l1_transformation.png")


def l2_macaulay():
    x = np.linspace(0, 4, 801)
    a = 1.5
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for n, lab in [
        (0, "⟨x − a⟩⁰ (step: point force in V, couple in M)"),
        (1, "⟨x − a⟩¹ (ramp: point force in M)"),
        (2, "⟨x − a⟩² / 2 (uniform load in M)"),
    ]:
        y = np.where(x >= a, (x - a) ** n if n else 1.0, 0.0) / (2 if n == 2 else 1)
        ax.plot(x, y, color=P.SERIES[n], label=lab)
    ax.axvline(a, color=P.GRID, lw=1)
    ax.text(a + 0.03, 3.3, "x = a", fontsize=8, color=P.TEXT_MUTED)
    ax.set_ylim(-0.2, 3.6)
    ax.set_xlabel("x")
    ax.set_title("Singularity functions switch on at x = a")
    ax.legend(loc="upper left", fontsize=8)
    return save(fig, "l2_macaulay.png")


def l3_hermite():
    xi = np.linspace(0, 1, 201)
    H = hermite(xi)
    labels = ["N1 (v1)", "N2 (L θ1)", "N3 (v2)", "N4 (L θ2)"]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for k in range(4):
        ax.plot(xi, H[k], color=P.SERIES[k], label=labels[k])
    ax.axhline(0, color=P.TEXT_MUTED, lw=0.8)
    ax.set_xlabel("ξ = x / L")
    ax.set_title("Cubic Hermite shape functions of the beam element")
    ax.legend(fontsize=8)
    return save(fig, "l3_hermite.png")


def l3_sparsity():
    m = Model()
    panels, a, h = 8, 2.0, 2.0
    bottom = [m.add_node(i * a, 0) for i in range(panels + 1)]
    top = [m.add_node(i * a, h) for i in range(panels + 1)]
    for i in range(panels):
        m.add_frame(bottom[i], bottom[i + 1], 200e9, 1e-3, 1e-6)
        m.add_frame(top[i], top[i + 1], 200e9, 1e-3, 1e-6)
        m.add_bar(bottom[i], top[i + 1], 200e9, 5e-4)
    for i in range(panels + 1):
        m.add_frame(bottom[i], top[i], 200e9, 1e-3, 1e-6)
    K = m.stiffness()
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    ax.spy(K, markersize=1.6, color=P.SERIES[0])
    ax.set_title(
        f"Global K: {K.shape[0]} DOFs, {K.nnz} non-zeros ({K.nnz / K.shape[0] ** 2:.1%})",
        fontsize=10,
    )
    ax.grid(False)
    return save(fig, "l3_sparsity.png")


def l4_imperfection():
    p = np.linspace(0, 0.98, 300)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for k, d0 in enumerate([0.001, 0.005, 0.02]):
        ax.plot(d0 / (1 - p), p, color=P.SERIES[k], label=f"initial bow δ₀ = {d0:g} L")
    ax.axhline(1, color=P.REFERENCE, ls="--", lw=1, label="Euler load (perfect column)")
    ax.set_xlim(0, 0.12)
    ax.set_xlabel("Midspan deflection δ / L")
    ax.set_ylabel("P / P_cr")
    ax.set_title("Imperfect columns: δ = δ₀ / (1 − P/P_cr)")
    ax.legend(fontsize=8, loc="lower right")
    return save(fig, "l4_imperfection.png")


def l5_qbar():
    ge = ply("T300/5208")
    th = np.linspace(-90, 90, 361)
    Q = np.array([qbar(ge, t * DEG) for t in th]) / 1e9
    fig, ax = plt.subplots(figsize=(7, 4))
    for k, (i, j, lab) in enumerate(
        [(0, 0, "Q̄11"), (1, 1, "Q̄22"), (0, 1, "Q̄12"), (2, 2, "Q̄66"), (0, 2, "Q̄16"), (1, 2, "Q̄26")]
    ):
        ax.plot(th, Q[:, i, j], color=P.SERIES[k], label=lab)
    ax.axhline(0, color=P.TEXT_MUTED, lw=0.8)
    ax.set_xlabel("Ply angle θ [deg]")
    ax.set_ylabel("GPa")
    ax.set_title("Transformed stiffness of a T300/5208 ply")
    ax.legend(ncol=3, fontsize=8)
    return save(fig, "l5_qbar.png")


def l5_off_axis():
    ge = ply("T300/5208")
    th = np.linspace(0.5, 90, 300)
    fig, ax = plt.subplots(figsize=(6.6, 4))
    for k, (lab, f) in enumerate(
        [
            ("Maximum stress", lambda s: failure.max_stress(s, ge)[0]),
            ("Tsai–Hill", lambda s: failure.tsai_hill(s, ge)),
            ("Tsai–Wu", lambda s: failure.tsai_wu(s, ge)),
        ]
    ):
        strength = [f(stress_transformation(t * DEG) @ np.array([1.0, 0, 0])) for t in th]
        ax.semilogy(th, np.array(strength) / MPA, color=P.SERIES[k], label=lab)
    ax.set_xlabel("Fibre angle to the load θ [deg]")
    ax.set_ylabel("Uniaxial tensile strength σx [MPa]")
    ax.set_title("Off-axis strength of a unidirectional lamina")
    ax.legend()
    return save(fig, "l5_off_axis.png")


def l6_rainflow():
    h = [-2, 1, -3, 5, -1, 3, -4, 4, -2]
    cycles = count_cycles(h)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(range(len(h)), h, color=P.TEXT_MUTED, lw=1, marker="o", ms=4, zorder=3)
    for c in cycles:
        col = P.SERIES[0] if c.count == 1 else P.SERIES[1]
        ax.plot([c.start, c.end], [h[c.start], h[c.end]], color=col, lw=3, alpha=0.6)
        f = 0.5 if c.count == 1 else 0.3
        ax.text(
            c.start + f * (c.end - c.start) - 0.12,
            h[c.start] + f * (h[c.end] - h[c.start]),
            f"{c.range:g}",
            fontsize=9,
            color=col,
            ha="right",
            fontweight="bold",
        )
    ax.plot([], [], color=P.SERIES[0], lw=3, alpha=0.6, label="full cycle")
    ax.plot([], [], color=P.SERIES[1], lw=3, alpha=0.6, label="half cycle")
    ax.set_xlabel("Reversal")
    ax.set_ylabel("Load")
    ax.set_title("ASTM E1049 example: the counted ranges")
    ax.legend(loc="lower left", fontsize=8)
    return save(fig, "l6_rainflow.png")


def l6_haigh():
    m = isotropic("2024-T3")
    su, sy = m.Ftu / KSI, m.Fty / KSI
    sar = 20.0  # fully reversed amplitude for the target life [ksi]
    sm = np.linspace(-10, su, 400)
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(sm, sar * (1 - sm / su), color=P.SERIES[0], label="Goodman")
    ax.plot(sm[sm >= 0], sar * (1 - (sm[sm >= 0] / su) ** 2), color=P.SERIES[1], label="Gerber")
    ax.plot(sm[sm <= sy], sar * (1 - sm[sm <= sy] / sy), color=P.SERIES[2], label="Soderberg")
    # SWT: sqrt(Smax Sa) = Sar -> Sa^2 + Sm Sa - Sar^2 = 0
    sa_swt = (-sm + np.sqrt(sm**2 + 4 * sar**2)) / 2
    ax.plot(sm, sa_swt, color=P.SERIES[3], label="Smith–Watson–Topper")
    c = handbook_curve("2024-T3_Kt1")
    N = c.life(sar * KSI, -1.0)
    R = np.linspace(-1, 0.9, 200)
    smax = c.max_stress_for_life(N, R) / KSI
    ok = smax <= su  # beyond Ftu the fitted equation is meaningless
    smax, R = smax[ok], R[ok]
    ax.plot(
        smax * (1 + R) / 2,
        smax * (1 - R) / 2,
        color=P.REFERENCE,
        ls="--",
        lw=1.2,
        label=f"MIL-HDBK-5J 2024-T3 curve, N = {N:.1e}",
    )
    ax.set_ylim(0, 26)
    ax.set_xlabel("Mean stress Sm [ksi]")
    ax.set_ylabel("Alternating stress Sa [ksi]")
    ax.set_title("Constant-life (Haigh) diagram, Sar = 20 ksi")
    ax.legend(fontsize=8)
    return save(fig, "l6_haigh.png")


def l7_orders():
    E, A, I, L = 70e9, 4e-4, 2e-6, 3.0
    ref_v = 5 * 2e3 * L**4 / (384 * E * I)
    ns = np.array([2, 4, 8, 16, 32])
    lumped, buck = [], []
    for n in ns:
        m = beam_model(L, n, E, A, I)
        m.fix(0, "xy").fix(n, "y")
        for e in range(n):
            m.distributed(e, -2e3, lumped=True)
        lumped.append(abs(-m.solve().displacement(n // 2)[1] - ref_v) / ref_v)
        P_fe = linear_buckling(column_model(1.2, n, E, 3e-4, 1e-8)).load_factors[0]
        buck.append(P_fe / euler_load(E, 1e-8, 1.2) - 1)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(1 / ns, lumped, "s-", color=P.SERIES[1], label="beam deflection, lumped loads")
    ax.loglog(1 / ns, buck, "o-", color=P.SERIES[0], label="buckling load, consistent K_G")
    h = 1 / ns
    ax.loglog(h, lumped[0] * (h / h[0]) ** 2, "--", color=P.REFERENCE, lw=1, label="h²")
    ax.loglog(h, buck[0] * (h / h[0]) ** 4, ":", color=P.REFERENCE, lw=1, label="h⁴")
    ax.set_xlabel("Element size h / L")
    ax.set_ylabel("Relative error")
    ax.set_title("Observed order of accuracy is a test, not a hope")
    ax.legend(fontsize=8)
    return save(fig, "l7_orders.png")


def main() -> None:
    P.use_style()
    for f in [
        l1_transformation,
        l2_macaulay,
        l3_hermite,
        l3_sparsity,
        l4_imperfection,
        l5_qbar,
        l5_off_axis,
        l6_rainflow,
        l6_haigh,
        l7_orders,
    ]:
        print("wrote", Path(f()).name)


if __name__ == "__main__":
    sys.exit(main())
