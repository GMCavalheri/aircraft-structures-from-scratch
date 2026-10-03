"""Phase 1 validation: stress transformation, principal stresses and yield criteria.

Writes ``docs/figures/mohr_circle.png``, ``docs/figures/yield_envelopes.png`` and
``validation/results/stress_strain.md``. Run with
``uv run python validation/stress_strain_validation.py``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from report import write_table

from structures import plotting as P
from structures.materials import KSI, isotropic
from structures.stress_strain import (
    MohrCircle,
    mohr_circles_3d,
    plot_mohr,
    principal_2d,
    principal_3d,
    stress_invariants,
    stress_tensor,
    transform_stress_2d,
    transform_tensor,
    tresca,
    von_mises,
)

DEG = np.pi / 180
MPA = 1e6


def checks():
    rows = []
    s = [-20 * MPA, 90 * MPA, 60 * MPA]
    s1, s2, tp = principal_2d(s)
    rows.append(
        [
            "Principal stresses, (−20, 90, 60) MPa",
            f"{s1 / MPA:.2f} / {s2 / MPA:.2f} MPa",
            f"{(35 + np.hypot(55, 60)):.2f} / {(35 - np.hypot(55, 60)):.2f} MPa",
            "closed form",
        ]
    )
    shear = abs(transform_stress_2d(s, tp)[2])
    rows.append(["Shear on principal plane", f"{shear:.1e} Pa", "0", "exact"])
    tau = 40 * MPA
    p1, p2, ang = principal_2d([0, 0, tau])
    rows.append(
        [
            "Pure shear: principal values, angle",
            f"±{p1 / MPA:.1f} MPa at {ang / DEG:.1f}°",
            "±40.0 MPa at 45°",
            "exact",
        ]
    )
    rng = np.random.default_rng(0)
    sigma = stress_tensor(*rng.normal(size=6) * 100 * MPA)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    a, b = (
        np.array(stress_invariants(sigma)),
        np.array(stress_invariants(transform_tensor(sigma, q))),
    )
    rows.append(
        [
            "Invariants under random rotation (max rel. change)",
            f"{np.max(np.abs((b - a) / a)):.1e}",
            "0",
            "round-off",
        ]
    )
    rows.append(
        [
            "von Mises, uniaxial / pure shear",
            f"{von_mises([1, 0, 0]):.4f} σ / {von_mises([0, 0, 1]):.4f} τ",
            "1 σ / √3 τ = 1.7321 τ",
            "exact",
        ]
    )
    rows.append(["Tresca, pure shear", f"{tresca([0, 0, 1]):.4f} τ", "2 τ", "exact"])
    return rows


def figure_mohr():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.4))
    c = MohrCircle(-20.0, 90.0, 60.0)
    plot_mohr(a1, c, color=P.SERIES[0], unit="MPa")
    a1.set_title("Plane stress: σx = −20, σy = 90, τxy = 60 MPa")
    sig = stress_tensor(80, 20, -40, 10, 0, 30)
    p, _ = principal_3d(sig)
    for i, (ctr, r) in enumerate(mohr_circles_3d(sig)):
        phi = np.linspace(0, np.pi, 181)
        a2.plot(ctr + r * np.cos(phi), r * np.sin(phi), color=P.SERIES[i], lw=1.5)
    a2.plot(p, np.zeros(3), "s", color=P.TEXT, ms=4)
    for k, v in enumerate(p):
        a2.annotate(
            f"σ{k + 1} = {v:.1f}", (v, 0), textcoords="offset points", xytext=(-12, -14), fontsize=8
        )
    a2.set_aspect("equal")
    a2.set_xlabel("Normal stress σ [MPa]")
    a2.set_ylabel("Shear stress |τ| [MPa]")
    a2.set_title("3D state: the three Mohr circles")
    fig.tight_layout()
    P.save(fig, "mohr_circle.png")


def figure_envelopes():
    m = isotropic("2024-T3", "B")
    fy = m.Fty / KSI
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    t = np.linspace(0, 2 * np.pi, 721)
    # von Mises ellipse s1^2 - s1 s2 + s2^2 = Fy^2, parametrised along its principal axes
    u, v = fy * np.sqrt(2) * np.cos(t), fy * np.sqrt(2 / 3) * np.sin(t)
    ax.plot((u - v) / np.sqrt(2), (u + v) / np.sqrt(2), color=P.SERIES[0], label="von Mises")
    hexagon = np.array([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1], [1, 0]]) * fy
    ax.plot(hexagon[:, 0], hexagon[:, 1], color=P.SERIES[1], label="Tresca")
    ax.plot([-fy, fy], [fy, -fy], color=P.GRID, lw=0.8)
    ax.annotate(
        "pure shear",
        (fy / np.sqrt(3), -fy / np.sqrt(3)),
        textcoords="offset points",
        xytext=(6, -2),
        fontsize=8,
        color=P.TEXT_MUTED,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("σ1 [ksi]")
    ax.set_ylabel("σ2 [ksi]")
    ax.set_title(f"Plane-stress yield envelopes, 2024-T3 (Fty = {fy:.0f} ksi)")
    ax.legend(loc="lower right")
    P.save(fig, "yield_envelopes.png")


if __name__ == "__main__":
    P.use_style()
    write_table(
        "stress_strain",
        "Stress/strain validation",
        ["Check", "Computed", "Reference", "Source"],
        checks(),
    )
    figure_mohr()
    figure_envelopes()
