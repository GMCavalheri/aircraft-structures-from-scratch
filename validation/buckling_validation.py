"""Phase 4 validation: FEA eigen-buckling vs Euler, column curves, plate buckling.

Writes ``docs/figures/buckling_modes.png``, ``docs/figures/column_curves.png``,
``docs/figures/plate_buckling_k.png`` and ``validation/results/buckling.md``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from reference import ramberg_osgood_n
from report import rel_err, write_table

from structures import plotting as P
from structures.buckling import (
    EFFECTIVE_LENGTH,
    column_model,
    compression_k,
    euler_load,
    euler_stress,
    johnson_stress,
    johnson_transition,
    linear_buckling,
    tangent_modulus_stress,
)
from structures.materials import KSI, isotropic

E, A, I, L = 70e9, 3e-4, 1e-8, 1.2


def _fixed_pinned(x):
    # v = A sin kx + B cos kx + C x + D with v(0) = v'(0) = 0 and v(1) = v''(1) = 0 gives
    # v = sin kx - kx + k (1 - cos kx), with tan k = k (k = 4.4934).
    k = np.pi / EFFECTIVE_LENGTH["fixed-pinned"]
    return np.sin(k * x) - k * x + k * (1 - np.cos(k * x))


# Exact first buckling modes on x in [0, 1] (Timoshenko & Gere, Theory of Elastic Stability, Ch. 2)
EXACT_MODE = {
    "pinned-pinned": lambda x: np.sin(np.pi * x),
    "fixed-free": lambda x: 1 - np.cos(np.pi * x / 2),
    "fixed-fixed": lambda x: (1 - np.cos(2 * np.pi * x)) / 2,
    "fixed-pinned": _fixed_pinned,
}


def _shape_error(end, n):
    """Max nodal deviation of the FEA mode from the exact mode, both scaled to max |v| = 1."""
    v = linear_buckling(column_model(L, n, E, A, I, end)).modes[1::3, 0]
    ref = EXACT_MODE[end](np.linspace(0, 1, n + 1))
    ref = ref / ref[np.argmax(np.abs(ref))]
    return float(np.max(np.abs(v - np.sign(np.dot(v, ref)) * ref)))


def modes_table_and_figure():
    rows = []
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.4), sharey=True)
    for ax, (k, end) in zip(axes, enumerate(EFFECTIVE_LENGTH), strict=True):
        ref = euler_load(E, I, L, end)
        errs = []
        for n in (4, 8, 16):
            lam = linear_buckling(column_model(L, n, E, A, I, end)).load_factors[0]
            errs.append(rel_err(lam, ref))
        rows.append(
            [end, f"{EFFECTIVE_LENGTH[end]:.4f}", f"{ref:.2f} N", *[f"{e:+.1e}" for e in errs]]
        )
        rows[-1].append(f"{_shape_error(end, 8):.1e}")
        res = linear_buckling(column_model(L, 40, E, A, I, end))
        x = np.linspace(0, 1, 41)
        v = res.modes[1::3, 0]
        ref_shape = EXACT_MODE[end](x)
        ref_shape = ref_shape / ref_shape[np.argmax(np.abs(ref_shape))]
        if np.dot(ref_shape, v) < 0:
            v = -v
        ax.plot(np.zeros_like(x), x, color=P.GRID, lw=1.0)
        ax.plot(v, x, color=P.SERIES[k], lw=2.5, label="FEA, 40 elements")
        ax.plot(ref_shape, x, color=P.REFERENCE, lw=1, ls="--", label="exact")
        ax.set_title(f"{end}\nK = {EFFECTIVE_LENGTH[end]:.3f}", fontsize=10)
        ax.set_xlabel("lateral deflection (normalised)")
        ax.set_xlim(-1.2, 1.2)
    axes[0].set_ylabel("x / L")
    axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle(
        "Buckling mode shapes from the FEA eigenvalue solver (40 frame elements)",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout()
    P.save(fig, "buckling_modes.png")
    return rows


def column_curves():
    fig, ax = plt.subplots(figsize=(7, 4.5))
    lam = np.linspace(1, 160, 600)
    rows = []
    for k, key in enumerate(["2024-T3", "7075-T6"]):
        m = isotropic(key)
        ax.plot(
            lam,
            johnson_stress(m.Ec, m.Fcy, lam) / KSI,
            color=P.SERIES[k],
            label=f"{key}: Johnson–Euler (Fcy = {m.Fcy / KSI:.0f} ksi)",
        )
        lt = johnson_transition(m.Ec, m.Fcy)
        ax.plot(lt, m.Fcy / 2 / KSI, "o", color=P.SERIES[k], ms=4)
        rows.append(
            [
                f"{key} Johnson–Euler transition",
                f"L'/ρ = {lt:.1f}",
                f"σ = Fcy/2 = {m.Fcy / 2 / KSI:.1f} ksi",
                "tangent",
            ]
        )
    m = isotropic("2024-T3")
    n, src = ramberg_osgood_n("2024-T3", "L", "compression")
    # MIL-HDBK-5J 1.6.2.2: the column stress may not exceed the compressive strength, so the
    # tangent-modulus curve (which keeps rising as L'/rho -> 0) is capped at Fcy.
    f_t = np.minimum(tangent_modulus_stress(m.Ec, m.Fcy, n, lam), m.Fcy)
    ax.plot(
        lam,
        f_t / KSI,
        color=P.SERIES[2],
        ls="-.",
        label=f"2024-T3: tangent modulus, n = {n:.0f}, capped at Fcy",
    )
    with np.errstate(divide="ignore"):
        ax.plot(
            lam,
            euler_stress(m.Ec, lam) / KSI,
            color=P.REFERENCE,
            ls="--",
            lw=1,
            label="2024-T3: Euler, π²E/(L'/ρ)²",
        )
    ax.set_ylim(0, 85)
    ax.set_xlabel("Effective slenderness L'/ρ")
    ax.set_ylabel("Column stress [ksi]")
    ax.set_title("Column curves with MIL-HDBK-5J B-basis allowables")
    ax.legend(fontsize=8)
    P.save(fig, "column_curves.png")
    return rows


def plate_k():
    r = np.linspace(0.3, 5, 800)
    fig, ax = plt.subplots(figsize=(7, 4))
    for m in range(1, 6):
        ax.plot(r, (m / r + r / m) ** 2, color=P.GRID, lw=1)
        ax.text(m * 1.0, 4.25, f"m = {m}", fontsize=8, ha="center", color=P.TEXT_MUTED)
    k, _ = compression_k(r)
    ax.plot(r, k, color=P.SERIES[0], lw=2, label="governing k (minimum over m)")
    ax.axhline(4, color=P.REFERENCE, ls="--", lw=1, label="k = 4 (long plate)")
    ax.set_ylim(0, 12)
    ax.set_xlabel("Aspect ratio a / b")
    ax.set_ylabel("Buckling coefficient k")
    ax.set_title("Simply supported plate, uniaxial compression")
    ax.legend(loc="upper right")
    P.save(fig, "plate_buckling_k.png")
    return [
        [
            "Plate k at a/b = 1, 2, 3",
            ", ".join(f"{compression_k(x)[0]:.4f}" for x in (1, 2, 3)),
            "4, 4, 4",
            "Timoshenko & Gere 9.2",
        ],
        [
            "Plate k at a/b = √2 (m = 1/2 crossover)",
            f"{compression_k(np.sqrt(2))[0]:.4f}",
            "4.5",
            "Timoshenko & Gere 9.2",
        ],
    ]


if __name__ == "__main__":
    P.use_style()
    write_table(
        "buckling",
        "Buckling validation: FEA eigenvalue vs Euler",
        [
            "End condition",
            "K",
            "Euler P_cr",
            "error, 4 el.",
            "8 el.",
            "16 el.",
            "mode-shape error, 8 el.",
        ],
        modes_table_and_figure(),
        notes="E = 70 GPa, I = 1e-8 m⁴, L = 1.2 m. Consistent geometric stiffness: the "
        "FEA load is an upper bound and the error falls as h⁴. On a uniform mesh the nodal "
        "deflections of the three trigonometric modes are exact to round-off; the fixed-pinned "
        "mode converges at about h⁶.",
    )
    write_table(
        "buckling_columns_plates",
        "Column curves and plate buckling",
        ["Check", "Computed", "Reference", "Source"],
        column_curves() + plate_k(),
    )
