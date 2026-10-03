import numpy as np
import pytest

from structures.buckling import (
    EFFECTIVE_LENGTH,
    column_model,
    compression_k,
    euler_load,
    euler_stress,
    johnson_stress,
    johnson_transition,
    linear_buckling,
    plate_critical_stress,
    shear_k,
    tangent_modulus,
    tangent_modulus_stress,
)
from structures.materials import isotropic

E, A, I, L = 70e9, 3e-4, 1e-8, 1.2


def test_fixed_pinned_factor_from_transcendental_root():
    # tan(kL) = kL first root 4.4934 -> K = pi / 4.4934
    x = np.pi / EFFECTIVE_LENGTH["fixed-pinned"]
    assert np.tan(x) == pytest.approx(x, rel=1e-9)
    assert EFFECTIVE_LENGTH["fixed-pinned"] == pytest.approx(0.6992, abs=1e-4)


@pytest.mark.parametrize("end", list(EFFECTIVE_LENGTH))
def test_fea_buckling_converges_to_euler(end):
    ref = euler_load(E, I, L, end)
    lam = linear_buckling(column_model(L, 16, E, A, I, end)).load_factors[0]
    # 16 elements: the error falls as h^4 (see the convergence-order test); fixed-fixed has the
    # shortest buckled half-wave (L/2) and is the least resolved, at 3e-5.
    assert lam == pytest.approx(ref, rel=1e-4)


def test_fea_buckling_convergence_order():
    ref = euler_load(E, I, L, "pinned-pinned")
    errs = [
        linear_buckling(column_model(L, n, E, A, I)).load_factors[0] / ref - 1 for n in (2, 4, 8)
    ]
    assert all(e > 0 for e in errs)  # consistent K_G: upper bound
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert rates == pytest.approx([4, 4], abs=0.15)  # cubic Hermite: error ~ h^4


def test_higher_modes_and_shapes_pinned():
    # Pinned-pinned: P_n = n^2 P_1 and mode n is sin(n pi x / L).
    res = linear_buckling(column_model(L, 40, E, A, I), n_modes=3)
    assert res.load_factors / res.load_factors[0] == pytest.approx([1, 4, 9], rel=1e-3)
    x = np.linspace(0, L, 41)
    for j, n in enumerate((1, 2, 3)):
        v = res.modes[1::3, j]
        ref = np.sin(n * np.pi * x / L)
        ref /= ref[np.argmax(np.abs(ref))]
        assert np.allclose(v, ref, atol=1e-3) or np.allclose(v, -ref, atol=1e-3)


def test_fixed_free_mode_shape():
    # Cantilever column: v = 1 - cos(pi x / 2L)
    res = linear_buckling(column_model(L, 30, E, A, I, "fixed-free"))
    x = np.linspace(0, L, 31)
    assert res.modes[1::3, 0] == pytest.approx(1 - np.cos(np.pi * x / (2 * L)), abs=1e-4)


def test_johnson_tangent_to_euler():
    m = isotropic("2024-T3")
    lam_t = johnson_transition(m.Ec, m.Fcy)
    j = johnson_stress(m.Ec, m.Fcy, [lam_t * (1 - 1e-9)])[0]
    assert j == pytest.approx(euler_stress(m.Ec, lam_t), rel=1e-6)
    assert j == pytest.approx(m.Fcy / 2, rel=1e-6)
    h = 1e-4 * lam_t
    dj = (johnson_stress(m.Ec, m.Fcy, [lam_t - h])[0] - j) / -h
    de = (euler_stress(m.Ec, lam_t + h) - euler_stress(m.Ec, lam_t)) / h
    assert dj == pytest.approx(de, rel=1e-3)
    assert johnson_stress(m.Ec, m.Fcy, [0.0])[0] == pytest.approx(m.Fcy)


def test_tangent_modulus_column():
    m = isotropic("2024-T3")
    n = 15  # MIL-HDBK-5J Fig. 3.2.3.1.6(a), 2024-T3 sheet, L-compression
    assert tangent_modulus(1e-3 * m.Fcy, m.Ec, m.Fcy, n) == pytest.approx(m.Ec, rel=1e-9)
    # Slender columns stay elastic: tangent-modulus stress -> Euler stress
    assert tangent_modulus_stress(m.Ec, m.Fcy, n, 200.0) == pytest.approx(
        euler_stress(m.Ec, 200.0), rel=1e-3
    )
    # Stocky columns: the tangent-modulus stress stays below Euler and satisfies Engesser
    f = tangent_modulus_stress(m.Ec, m.Fcy, n, 30.0)
    assert f < euler_stress(m.Ec, 30.0)
    assert f == pytest.approx(np.pi**2 * tangent_modulus(f, m.Ec, m.Fcy, n) / 30.0**2)


def test_plate_compression_k():
    k, mm = compression_k(1.0)
    assert (k, mm) == (pytest.approx(4.0), 1)
    k, mm = compression_k(3.0)
    assert (k, mm) == (pytest.approx(4.0), 3)
    # m = 1 and m = 2 curves cross at a/b = sqrt(2) with k = 4.5
    assert compression_k(np.sqrt(2))[0] == pytest.approx(4.5)
    ks, _ = compression_k(np.linspace(0.5, 10, 500))
    assert ks.min() == pytest.approx(4.0, abs=1e-3)
    assert compression_k(0.5)[0] == pytest.approx((2 + 0.5) ** 2)


def test_plate_shear_and_stress():
    assert shear_k(1.0) == pytest.approx(9.35)
    assert shear_k(1e6) == pytest.approx(5.35)
    m = isotropic("2024-T3")
    s = plate_critical_stress(m.Ec, m.nu, 0.001, 0.15, 4.0)
    assert s == pytest.approx(4 * np.pi**2 * m.Ec / (12 * (1 - m.nu**2)) * (1 / 150) ** 2)
