"""Phase 5 validation: Classical Lamination Theory against Kaw's worked examples, plus
ply-by-ply stresses and failure envelopes for a MIL-HDBK-17-2F carbon/epoxy laminate.

Writes ``docs/figures/laminate_ply_stresses.png``, ``docs/figures/laminate_failure_envelope.png``,
``docs/figures/laminate_angle_ply.png`` and ``validation/results/composite_laminates.md``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from report import rel_err, write_table

from structures import plotting as P
from structures.composite_laminates import Laminate, failure, reduced_compliance
from structures.composite_laminates import stress_transformation as T
from structures.materials import ply

DEG = np.pi / 180
GE = ply("T300/5208")
HT = ply("T300/976")
T_PLY_976 = 0.0053 * 0.0254  # MIL-HDBK-17-2F Table 4.2.17(a): cured ply thickness 0.0053 in


def kaw_rows():
    rows = []

    def add(name, computed, ref, fmt="{:.4g}"):
        rows.append([name, fmt.format(computed), fmt.format(ref), f"{rel_err(computed, ref):+.1e}"])

    lam = Laminate.from_angles(GE, [0, 90, 0], 0.005)
    A, _, D = lam.stiffness()
    add("[0/90/0] A11 [Pa·m]", A[0, 0], 1.870e9)
    add("[0/90/0] A22 [Pa·m]", A[1, 1], 1.013e9)
    add("[0/90/0] D11 [Pa·m³]", D[0, 0], 4.935e4)
    add("[0/90/0] D22 [Pa·m³]", D[1, 1], 4.696e3)
    c = lam.engineering_constants()
    add("[0/90/0] in-plane Ex [GPa]", c["Ex"] / 1e9, 124.5)
    add("[0/90/0] in-plane Ey [GPa]", c["Ey"] / 1e9, 67.43)
    add("[0/90/0] in-plane νxy", c["nuxy"], 0.04292)
    add("[0/90/0] flexural Ex [GPa]", c["Ex_flex"] / 1e9, 175.0)
    add("[0/90/0] flexural Ey [GPa]", c["Ey_flex"] / 1e9, 16.65)
    pts = {pt.ply: pt for pt in lam.ply_points(N=(1, 0, 0))}
    add("[0/90/0], Nx = 1 N/m: σ1 in 0° ply [Pa]", pts[0].stress_12[0], 97.26)
    add("[0/90/0], Nx = 1 N/m: σ2 in 90° ply [Pa]", pts[1].stress_12[1], 5.472)
    s12 = T(60 * DEG) @ np.array([2.0, -3.0, 4.0])
    add("60° lamina, max stress S [MPa]", failure.max_stress(s12, GE)[0] / 1e6, 16.33)
    add(
        "60° lamina, max strain S [MPa]",
        failure.max_strain(reduced_compliance(GE) @ s12, GE)[0] / 1e6,
        16.33,
    )
    add("60° lamina, Tsai-Hill S [MPa]", failure.tsai_hill(s12, GE) / 1e6, 10.94)
    add("60° lamina, modified Tsai-Hill S [MPa]", failure.tsai_hill(s12, GE, True) / 1e6, 16.06)
    add("60° lamina, Tsai-Wu S [MPa]", failure.tsai_wu(s12, GE) / 1e6, 22.39)
    hist = lam.ply_by_ply_failure(N=(1, 0, 0))
    add("[0/90/0] first ply failure Nx [N/m] (90°, Tsai-Wu)", hist[0][0], 7.277e6)
    add("[0/90/0] last ply failure Nx [N/m] (0°, discounted)", hist[1][0], 1.5e7)
    return rows


def quasi_iso():
    return Laminate.from_angles(HT, [0, 45, -45, 90], T_PLY_976, symmetric=True)


def figure_ply_stresses():
    lam = quasi_iso()
    N, M = (150e3, 0, 0), (8.0, 0, 0)
    pts = lam.ply_points(N, M)
    z = np.array([p.z for p in pts]) * 1e3
    fig, axes = plt.subplots(1, 3, figsize=(10, 4.6), sharey=True)
    labels = [
        ("σx (global)", lambda p: p.stress_xy[0]),
        ("σ1 (fibre)", lambda p: p.stress_12[0]),
        ("σ2 (transverse)", lambda p: p.stress_12[1]),
    ]
    angles = [round(p.angle / DEG) for p in lam.plies]
    colors = {0: P.SERIES[0], 45: P.SERIES[1], -45: P.SERIES[2], 90: P.SERIES[3]}
    for ax, (lab, f) in zip(axes, labels, strict=True):
        vals = np.array([f(p) for p in pts]) / 1e6
        for k in range(len(lam.plies)):
            sl = slice(3 * k, 3 * k + 3)
            ax.plot(
                vals[sl],
                z[sl],
                color=colors[angles[k]],
                lw=2.5,
                label=(f"{angles[k]:+d}°" if angles[k] % 90 else f"{angles[k]}°")
                if k < 4
                else None,
            )
        for zi in lam.z * 1e3:
            ax.axhline(zi, color=P.GRID, lw=0.8)
        ax.axvline(0, color=P.TEXT_MUTED, lw=0.8)
        ax.set_xlabel(f"{lab} [MPa]")
    axes[0].set_ylabel("z [mm]")
    axes[0].legend(title="ply", fontsize=8, loc="lower left")
    fig.suptitle(
        "T300/976 [0/±45/90]s under Nx = 150 kN/m and Mx = 8 N·m/m: ply-by-ply stress",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    P.save(fig, "laminate_ply_stresses.png")


def figure_envelope():
    lam = quasi_iso()
    phi = np.linspace(0, 2 * np.pi, 361)
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    rows = []
    for k, (crit, label) in enumerate(
        [("tsai_wu", "Tsai-Wu"), ("max_stress", "Maximum stress"), ("tsai_hill", "Tsai-Hill")]
    ):
        r = np.array(
            [lam.first_ply_failure(N=(np.cos(p), np.sin(p), 0), criterion=crit)[0] for p in phi]
        )
        ax.plot(r * np.cos(phi) / 1e3, r * np.sin(phi) / 1e3, color=P.SERIES[k], label=label)
        rows.append(r)
    ax.axhline(0, color=P.TEXT_MUTED, lw=0.8)
    ax.axvline(0, color=P.TEXT_MUTED, lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("Nx [kN/m]")
    ax.set_ylabel("Ny [kN/m]")
    ax.set_title("First-ply-failure envelope, T300/976 [0/±45/90]s")
    ax.legend(loc="upper left", fontsize=8)
    P.save(fig, "laminate_failure_envelope.png")
    tw = rows[0]
    h = lam.thickness
    return [
        [
            "T300/976 [0/±45/90]s FPF, Nx tension (Tsai-Wu)",
            f"{tw[0] / 1e3:.0f} kN/m",
            f"{tw[0] / h / 1e6:.0f} MPa average",
            "MIL-HDBK-17-2F ply data",
        ],
        [
            "T300/976 [0/±45/90]s FPF, Nx compression (Tsai-Wu)",
            f"{tw[180] / 1e3:.0f} kN/m",
            f"{tw[180] / h / 1e6:.0f} MPa average",
            "MIL-HDBK-17-2F ply data",
        ],
    ]


def figure_angle_ply():
    th = np.linspace(0, 90, 91)
    Ex, Gxy, nu = [], [], []
    for t in th:
        c = Laminate.from_angles(HT, [t, -t], T_PLY_976, symmetric=True).engineering_constants()
        Ex.append(c["Ex"])
        Gxy.append(c["Gxy"])
        nu.append(c["nuxy"])
    fig, a1 = plt.subplots(figsize=(6.6, 4.2))
    a1.plot(th, np.array(Ex) / 1e9, color=P.SERIES[0], label="Ex")
    a1.plot(th, np.array(Gxy) / 1e9, color=P.SERIES[1], label="Gxy")
    a1.set_xlabel("Ply angle θ in [±θ]s [deg]")
    a1.set_ylabel("Modulus [GPa]")
    a2 = a1.twinx()
    a2.plot(th, nu, color=P.SERIES[2], label="νxy")
    a2.set_ylabel("Poisson's ratio νxy")
    a2.grid(False)
    a2.spines["right"].set_visible(True)
    lines = a1.get_lines() + a2.get_lines()
    a1.legend(lines, [ln.get_label() for ln in lines], loc="upper right")
    a1.set_title("Angle-ply laminate constants, T300/976")
    P.save(fig, "laminate_angle_ply.png")
    i = int(np.argmax(nu))
    return [
        [
            "[±θ]s peak νxy",
            f"{nu[i]:.2f} at θ = {th[i]:.0f}°",
            "νxy > 1 is possible for angle plies",
            "CLT",
        ],
        ["[±45]s Gxy / G12", f"{Gxy[45] / HT.G12:.2f}", "max of Gxy over θ", "CLT"],
    ]


if __name__ == "__main__":
    P.use_style()
    write_table(
        "composite_laminates",
        "Composite laminates: CLT vs Kaw's worked examples",
        ["Quantity", "Computed", "Kaw", "Rel. error"],
        kaw_rows(),
        notes="T300/5208 graphite/epoxy (Kaw Table 2.1), 5 mm plies as in Kaw's examples. "
        "Kaw prints four significant figures, so ~1e-3 agreement is the best "
        "possible.",
    )
    figure_ply_stresses()
    write_table(
        "composite_laminates_demo",
        "T300/976 laminate studies",
        ["Quantity", "Computed", "Note", "Data"],
        figure_envelope() + figure_angle_ply(),
    )
