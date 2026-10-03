"""Phase 6 validation: rainflow counting, MIL-HDBK-5J S-N equations, curve fitting and
Miner's-rule life of a synthetic flight spectrum.

Writes ``docs/figures/sn_curves.png``, ``docs/figures/sn_fit.png``,
``docs/figures/spectrum_rainflow.png``, ``docs/figures/miner_life.png`` and
``validation/results/fatigue.md``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from report import write_table

from structures import plotting as P
from structures.fatigue import (
    blocks_to_failure,
    count_cycles,
    fit_basquin,
    fit_equivalent_stress,
    flight_history,
    handbook_curve,
    histogram,
    mean_stress,
    miner_damage,
    range_mean_matrix,
)
from structures.materials import KSI, isotropic

ASTM = [-2, 1, -3, 5, -1, 3, -4, 4, -2]
SPECTRUM = {"s_1g": 12 * KSI, "s_ground": -4 * KSI, "beta": 2.5 * KSI, "cycles_per_flight": 40}
N_FLIGHTS = 500


def rainflow_rows():
    got = histogram(count_cycles(ASTM))
    ref = [(3, 0.5), (4, 1.5), (6, 0.5), (8, 1.0), (9, 0.5)]
    return [
        [
            "ASTM E1049 example, (range, count)",
            "; ".join(f"{r:g}: {c:g}" for r, c in got),
            "; ".join(f"{r:g}: {c:g}" for r, c in ref),
            "✔ exact" if got == ref else "✘",
        ]
    ]


def figure_sn():
    # Stop at 1e7 cycles: beyond the tested range the fitted equations are extrapolations (the
    # handbook's own caution) and the 7075-T6 curves start to cross.
    N = np.geomspace(1e3, 1e7, 300)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, mat in zip(axes, ("2024-T3", "7075-T6"), strict=True):
        for kt, ls in (("Kt1", "-"), ("Kt2", "--")):
            c = handbook_curve(f"{mat}_{kt}")
            for k, R in enumerate((-1.0, 0.0, 0.5)):
                s = c.max_stress_for_life(N, R) / KSI
                ax.semilogx(
                    N, s, color=P.SERIES[k], ls=ls, label=f"R = {R:g}" if kt == "Kt1" else None
                )
        ax.set_title(f"{mat} sheet (MIL-HDBK-5J)")
        ax.set_xlabel("Cycles to failure N")
        ax.set_ylim(0, 80)
    axes[0].set_ylabel("Maximum stress Smax [ksi]")
    axes[0].plot([], [], color=P.TEXT_MUTED, ls="-", label="unnotched, Kt = 1")
    axes[0].plot([], [], color=P.TEXT_MUTED, ls="--", label="notched, Kt = 2")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    P.save(fig, "sn_curves.png")
    c = handbook_curve("2024-T3_Kt1")
    computed = float(np.log10(c.life(40 * KSI, 0.0)))
    return [
        [
            "2024-T3 Kt = 1, Smax = 40 ksi, R = 0: log Nf",
            f"{computed:.4f}",
            f"{11.1 - 3.97 * np.log10(24.2):.4f}",
            "hand calc. of Fig. 3.2.3.1.8(e)",
        ]
    ]


def figure_fit():
    true = handbook_curve("2024-T3_Kt2")
    rng = np.random.default_rng(7)
    R = np.repeat([-1.0, 0.0, 0.5], 40)
    s_max = true.max_stress_for_life(10 ** rng.uniform(4, 7, R.size), R)
    N = true.life(s_max, R) * 10 ** rng.normal(0, true.std_error, R.size)
    fit = fit_equivalent_stress(s_max, R, N)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    seq = true.equivalent_stress(s_max, R) / KSI
    for k, r in enumerate((-1.0, 0.0, 0.5)):
        sel = R == r
        ax.semilogx(
            N[sel],
            seq[sel],
            "o",
            ms=4,
            color=P.SERIES[k],
            alpha=0.8,
            label=f"synthetic tests, R = {r:g}",
        )
    Ns = np.geomspace(1e3, 1e8, 300)
    ax.semilogx(
        Ns,
        true.max_stress_for_life(Ns, 0.0) / KSI,
        color=P.REFERENCE,
        ls="--",
        label="MIL-HDBK-5J curve (Seq form)",
    )
    ax.semilogx(
        Ns,
        (10 ** ((np.log10(Ns) - fit.A1) / fit.A2) + fit.A4),
        color=P.SERIES[3],
        label="refitted from the synthetic data",
    )
    ax.set_xlabel("Cycles to failure N")
    ax.set_ylabel("Equivalent stress Seq = Smax(1−R)^A3 [ksi]")
    ax.set_title("2024-T3 Kt = 2: all stress ratios collapse onto one Seq curve")
    ax.legend(fontsize=8)
    P.save(fig, "sn_fit.png")
    return [
        [
            "Refit of 120 synthetic tests (handbook scatter): A1, A2, A3, A4",
            f"{fit.A1:.2f}, {fit.A2:.2f}, {fit.A3:.3f}, {fit.A4:.1f}",
            f"{true.A1:.2f}, {true.A2:.2f}, {true.A3:.3f}, {true.A4:.1f}",
            "Fig. 3.2.3.1.8(g)",
        ],
        [
            "Refit standard error of log life",
            f"{fit.std_error:.3f}",
            f"{true.std_error:.3f}",
            "Fig. 3.2.3.1.8(g)",
        ],
    ]


def spectrum_study():
    hist = flight_history(N_FLIGHTS, seed=11, **SPECTRUM)
    cycles = count_cycles(hist)
    hb = handbook_curve("2024-T3_Kt2")
    m = isotropic("2024-T3")
    # Basquin fitted to the handbook curve at R = -1 between 1e4 and 1e7 cycles (amplitude form)
    Nfit = np.geomspace(1e4, 1e7, 30)
    bas = fit_basquin(hb.max_stress_for_life(Nfit, -1.0), Nfit)

    def corrected(correction):
        return lambda a, s: _pos(bas.life, correction(a, s))

    models = {
        "MIL-HDBK-5J equivalent stress (R-dependent curve)": hb.life_from_amplitude,
        "Basquin (R = −1 fit) + Goodman": corrected(lambda a, s: mean_stress.goodman(a, s, m.Ftu)),
        "Basquin (R = −1 fit) + SWT": corrected(mean_stress.smith_watson_topper),
        "Basquin (R = −1 fit) + Walker, γ = A3": corrected(
            lambda a, s: mean_stress.walker(a, s, hb.A3)
        ),
    }
    rows = []
    lives = {}
    for name, life in models.items():
        D = miner_damage(cycles, life)
        flights = blocks_to_failure(D) * N_FLIGHTS
        lives[name] = flights
        rows.append([name, f"{D:.3e}", f"{flights:,.0f}"])
    figure_spectrum(hist, cycles)
    figure_lives(lives)
    return rows, len(cycles)


def _pos(fn, s):
    s = np.asarray(s, float)
    with np.errstate(divide="ignore"):
        return np.where(s > 0, fn(np.maximum(s, 1e-30)), np.inf)


def figure_spectrum(hist, cycles):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [1.4, 1]})
    n = 3 * (2 + 2 * SPECTRUM["cycles_per_flight"])
    a1.plot(np.arange(n), hist[:n] / KSI, color=P.SERIES[0], lw=0.9)
    a1.axhline(SPECTRUM["s_1g"] / KSI, color=P.TEXT_MUTED, lw=0.8, ls="--")
    a1.set_xlabel("Reversal number")
    a1.set_ylabel("Net-section stress [ksi]")
    a1.set_title("Synthetic flight history (first three flights)")
    r_bins = np.linspace(0, 32, 17)
    m_bins = np.linspace(-2, 22, 13)
    H = range_mean_matrix(cycles, np.array(r_bins) * KSI, np.array(m_bins) * KSI)
    im = a2.pcolormesh(
        m_bins,
        r_bins,
        np.where(H > 0, H, np.nan),
        cmap="Blues",
        norm=LogNorm(),
    )
    fig.colorbar(im, ax=a2, label="cycles per 500 flights")
    a2.set_xlabel("Cycle mean [ksi]")
    a2.set_ylabel("Cycle range [ksi]")
    a2.set_title("Rainflow range–mean matrix")
    a2.grid(False)
    fig.tight_layout()
    P.save(fig, "spectrum_rainflow.png")


def figure_lives(lives):
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    names = list(lives)
    vals = [lives[k] for k in names]
    ax.barh(range(len(names)), vals, color=[P.SERIES[i] for i in range(len(names))])
    ax.set_yticks(range(len(names)), names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Predicted median life [flights] (Miner, D = 1)")
    ax.set_title("2024-T3 Kt = 2 under the synthetic spectrum")
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:,.0f}", va="center", fontsize=8)
    ax.set_xlim(0, max(vals) * 1.25)
    P.save(fig, "miner_life.png")


if __name__ == "__main__":
    P.use_style()
    rows = rainflow_rows() + figure_sn() + figure_fit()
    write_table(
        "fatigue",
        "Fatigue validation",
        ["Check", "Computed", "Reference", "Source"],
        rows,
        notes="A1, A2 and A4 are strongly correlated (a steeper slope with a higher "
        "fatigue limit describes nearly the same curve), so the individual "
        "coefficients wander while the refitted curve stays within a factor of 1.5 "
        "in life of the handbook curve over 10⁴–10⁶ cycles (tested).",
    )
    life_rows, n_cycles = spectrum_study()
    write_table(
        "fatigue_spectrum",
        "Miner's-rule life of a synthetic flight spectrum",
        ["S-N / mean-stress model", "Damage per 500 flights", "Life [flights]"],
        life_rows,
        notes=f"{n_cycles} rainflow cycles per 500 flights; 1 g stress 12 ksi, ground "
        "−4 ksi, exponential gust increments with mean 2.5 ksi, 40 gust cycles per "
        "flight. Median lives (no scatter factor).",
    )
