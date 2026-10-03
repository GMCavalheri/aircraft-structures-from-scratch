import numpy as np
import pytest

from structures.fatigue import (
    BasquinCurve,
    Cycle,
    blocks_to_failure,
    count_cycles,
    fit_basquin,
    fit_equivalent_stress,
    flight_history,
    handbook_curve,
    histogram,
    mean_stress,
    miner_damage,
    reversals,
)
from structures.materials import KSI

# ASTM E1049-85 (2017), Sec. 5.4.4 / Fig. 6 rainflow example (also the first test case of the
# open-source `rainflow` Python package).
ASTM = [-2, 1, -3, 5, -1, 3, -4, 4, -2]


def test_rainflow_astm_example_histogram():
    assert histogram(count_cycles(ASTM)) == [(3, 0.5), (4, 1.5), (6, 0.5), (8, 1.0), (9, 0.5)]


def test_rainflow_astm_example_cycles():
    expected = [
        (3, -0.5, 0.5, 0, 1),
        (4, -1.0, 0.5, 1, 2),
        (4, 1.0, 1.0, 4, 5),
        (8, 1.0, 0.5, 2, 3),
        (9, 0.5, 0.5, 3, 6),
        (8, 0.0, 0.5, 6, 7),
        (6, 1.0, 0.5, 7, 8),
    ]
    got = [(c.range, c.mean, c.count, c.start, c.end) for c in count_cycles(ASTM)]
    assert got == expected


def test_reversals_and_conservation():
    x = [0, 1, 2, 2, 1, -1, -1, 3, 3, 0]
    v, i = reversals(x)
    assert v.tolist() == [0, 2, -1, 3, 0]
    assert i.tolist() == [0, 2, 5, 7, 9]
    # every range between successive reversals is counted exactly once (as two halves or a pair)
    rng = np.random.default_rng(4)
    series = rng.normal(size=500)
    v, _ = reversals(series)
    cycles = count_cycles(series)
    assert sum(2 * c.count for c in cycles) == pytest.approx(len(v) - 1)


def test_rainflow_constant_amplitude():
    t = np.linspace(0, 20 * np.pi, 2001)  # 10 full sine cycles
    cycles = count_cycles(3 + 2 * np.sin(t))
    full = [c for c in cycles if c.count == 1.0]
    assert sum(c.count for c in cycles) == pytest.approx(10.0, abs=0.5)
    assert all(c.range == pytest.approx(4, rel=1e-4) for c in full)
    assert all(c.mean == pytest.approx(3, abs=1e-4) for c in full)


def test_miner_constant_amplitude():
    curve = BasquinCurve(C=1e12, m=4)
    cycles = [Cycle(range=200.0, mean=0.0, count=1.0, start=0, end=1)] * 5000
    D = miner_damage(cycles, lambda a, m: curve.life(a))
    assert D == pytest.approx(5000 / curve.life(100.0))
    assert blocks_to_failure(D) == pytest.approx(curve.life(100.0) / 5000)
    assert blocks_to_failure(0.0) == np.inf


def test_basquin_fit_recovers_parameters():
    S = np.geomspace(100e6, 400e6, 12)
    true = BasquinCurve(C=3e40, m=4.2)
    fit = fit_basquin(S, true.life(S))
    assert fit.m == pytest.approx(4.2)
    assert fit.C == pytest.approx(3e40, rel=1e-8)
    assert fit.strength(fit.life(250e6)) == pytest.approx(250e6)


def test_handbook_equation_hand_calculation():
    # MIL-HDBK-5J Fig. 3.2.3.1.8(e), 2024-T3 unnotched: log Nf = 11.1 - 3.97 log(Seq - 15.8),
    # Seq = Smax (1 - R)^0.56. Smax = 40 ksi, R = 0: log Nf = 11.1 - 3.97 log10(24.2) = 5.6062.
    c = handbook_curve("2024-T3_Kt1")
    assert np.log10(c.life(40 * KSI, 0.0)) == pytest.approx(11.1 - 3.97 * np.log10(24.2))
    # fully reversed: Seq = 2^0.56 Smax
    assert c.equivalent_stress(30 * KSI, -1.0) / KSI == pytest.approx(30 * 2**0.56)
    # at or below the fatigue limit the life is infinite
    assert c.life(15.8 * KSI, 0.0) == np.inf
    assert c.life(10 * KSI, 0.0) == np.inf
    # inverse
    assert c.max_stress_for_life(1e6, 0.1) == pytest.approx(
        ((10 ** ((6 - 11.1) / -3.97) + 15.8) / 0.9**0.56) * KSI
    )


def test_handbook_curves_ordering():
    # A notch (Kt = 2) shortens life at the same net-section stress; higher R (more mean) at
    # the same Smax lengthens it.
    for mat in ("2024-T3", "7075-T6"):
        un, no = handbook_curve(f"{mat}_Kt1"), handbook_curve(f"{mat}_Kt2")
        assert no.life(35 * KSI, 0.1) < un.life(35 * KSI, 0.1)
        assert un.life(40 * KSI, 0.5) > un.life(40 * KSI, 0.0)


def test_life_from_amplitude_consistent_with_R_form():
    c = handbook_curve("7075-T6_Kt2")
    s_a, s_m = 12 * KSI, 18 * KSI
    R = (s_m - s_a) / (s_m + s_a)
    assert c.life_from_amplitude(s_a, s_m) == pytest.approx(c.life(s_m + s_a, R))
    assert c.life_from_amplitude(10 * KSI, -20 * KSI) == np.inf  # wholly compressive


def test_fit_equivalent_stress_recovers_handbook_coefficients():
    # Synthetic test programme from the handbook curve, with the handbook's own scatter
    # (standard error of log life), refitted: coefficients come back within their scatter.
    true = handbook_curve("2024-T3_Kt2")
    rng = np.random.default_rng(7)
    R = np.repeat([-1.0, 0.0, 0.5], 40)
    life_target = 10 ** rng.uniform(4, 7, R.size)
    s_max = true.max_stress_for_life(life_target, R)
    N = true.life(s_max, R) * 10 ** rng.normal(0, true.std_error, R.size)
    fit = fit_equivalent_stress(s_max, R, N)
    assert fit.A3 == pytest.approx(true.A3, abs=0.05)
    assert fit.std_error == pytest.approx(true.std_error, rel=0.25)
    # predicted lives agree within a factor of 1.5 over the fitted range
    for r in (-1.0, 0.0, 0.5):
        for n in (1e4, 1e5, 1e6):
            s = true.max_stress_for_life(n, r)
            assert abs(np.log10(fit.life(s, r) / n)) < np.log10(1.5)


def test_mean_stress_corrections():
    s_a, s_u, s_y = 100e6, 480e6, 330e6
    for f in (
        lambda a, m: mean_stress.goodman(a, m, s_u),
        lambda a, m: mean_stress.gerber(a, m, s_u),
        lambda a, m: mean_stress.soderberg(a, m, s_y),
        mean_stress.smith_watson_topper,
        lambda a, m: mean_stress.walker(a, m, 0.6),
    ):
        assert f(s_a, 0.0) == pytest.approx(s_a)  # all reduce to S_a at zero mean
    assert mean_stress.goodman(s_a, s_u / 2, s_u) == pytest.approx(2 * s_a)
    assert mean_stress.walker(s_a, 50e6, 0.5) == pytest.approx(
        mean_stress.smith_watson_topper(s_a, 50e6)
    )
    # Gerber is the least conservative for tensile means, Soderberg the most
    sm = 150e6
    assert (
        mean_stress.gerber(s_a, sm, s_u)
        < mean_stress.goodman(s_a, sm, s_u)
        < mean_stress.soderberg(s_a, sm, s_y)
    )


def test_walker_matches_handbook_equivalent_stress():
    # Seq = Smax (1 - R)^A3 = 2^A3 Smax^(1 - A3) Sa^A3 = 2^A3 * Walker(Sa, Sm, A3)
    c = handbook_curve("2024-T3_Kt1")
    s_a, s_m = 15 * KSI, 10 * KSI
    R = (s_m - s_a) / (s_m + s_a)
    assert c.equivalent_stress(s_a + s_m, R) == pytest.approx(
        2**c.A3 * mean_stress.walker(s_a, s_m, c.A3)
    )


def test_flight_history_reproducible_and_bounded():
    a = flight_history(5, 80e6, -20e6, 10e6, seed=3)
    b = flight_history(5, 80e6, -20e6, 10e6, seed=3)
    assert np.array_equal(a, b)
    assert a.size == 5 * (2 + 2 * 40)
    # each flight contributes one ground-air-ground cycle from -20 MPa up to the flight maximum
    cycles = count_cycles(a)
    assert max(c.range for c in cycles) >= 100e6
