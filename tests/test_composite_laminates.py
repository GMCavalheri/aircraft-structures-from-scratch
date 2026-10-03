"""CLT checks. [Kaw] = Kaw, Mechanics of Composite Materials, 2nd ed. (CRC, 2006), worked
examples as reproduced in the author's open courseware slides (USF compositesOCW):
Example 2.x (60 deg lamina failure), Ex. 4.x ([0/90/0] laminate moduli), Ex. 5.x (ply-by-ply
failure of [0/90/0])."""

import numpy as np
import pytest

from structures.composite_laminates import (
    Laminate,
    Ply,
    failure,
    invariants,
    qbar,
    qbar_from_invariants,
    reduced_compliance,
    reduced_stiffness,
    strain_transformation,
    stress_transformation,
)
from structures.materials import OrthotropicPly, ply

GE = ply("T300/5208")  # Kaw Table 2.1 graphite/epoxy
DEG = np.pi / 180


def kaw_090():
    return Laminate.from_angles(GE, [0, 90, 0], 0.005)


def test_reduced_stiffness_kaw():
    # [Kaw] [Q] = [181.8 2.897 0; 2.897 10.35 0; 0 0 7.17] GPa
    Q = reduced_stiffness(GE) / 1e9
    assert Q[0, 0] == pytest.approx(181.8, abs=0.05)
    assert Q[0, 1] == pytest.approx(2.897, abs=5e-4)
    assert Q[1, 1] == pytest.approx(10.35, abs=5e-3)
    assert Q[2, 2] == pytest.approx(7.17)
    assert reduced_compliance(GE) @ reduced_stiffness(GE) == pytest.approx(np.eye(3))


def test_transformations():
    # [Kaw] 60 deg lamina, global (2, -3, 4) -> local (1.714, -2.714, -4.165)
    s12 = stress_transformation(60 * DEG) @ np.array([2.0, -3.0, 4.0])
    assert s12 == pytest.approx([1.714, -2.714, -4.165], abs=1e-3)
    # 0 deg: identity; 90 deg: swap 1 and 2 and flip shear
    assert qbar(GE, 0) == pytest.approx(reduced_stiffness(GE))
    Q, Q90 = reduced_stiffness(GE), qbar(GE, 90 * DEG)
    assert Q90[0, 0] == pytest.approx(Q[1, 1])
    assert Q90[1, 1] == pytest.approx(Q[0, 0])
    assert Q90[0, 2] == pytest.approx(0, abs=1e-3)
    # work invariance s12 . e12 = sxy . exy requires T_sigma^T T_eps = I
    t = 0.37
    assert stress_transformation(t).T @ strain_transformation(t) == pytest.approx(np.eye(3))


def test_qbar_invariants_and_symmetry():
    U = invariants(GE)
    for th in np.radians([0, 15, 30, 45, 60, 75, 90, -30]):
        assert qbar(GE, th) == pytest.approx(qbar_from_invariants(U, th), rel=1e-12, abs=1)
        assert np.allclose(qbar(GE, th), qbar(GE, th).T)
    # +theta and -theta differ only in the sign of the coupling terms Q16, Q26
    a, b = qbar(GE, 30 * DEG), qbar(GE, -30 * DEG)
    assert a[0, 2] == pytest.approx(-b[0, 2])
    assert a[0, 0] == pytest.approx(b[0, 0])


def test_kaw_090_extensional_and_bending():
    lam = kaw_090()
    A, B, D = lam.stiffness()
    # [Kaw] A = [1.870e9 4.345e7 0; 4.345e7 1.013e9 0; 0 0 1.076e8] Pa m
    assert A[0, 0] == pytest.approx(1.870e9, rel=1e-3)
    assert A[0, 1] == pytest.approx(4.345e7, rel=1e-3)
    assert A[1, 1] == pytest.approx(1.013e9, rel=1e-3)
    assert A[2, 2] == pytest.approx(1.076e8, rel=1e-3)
    assert B == pytest.approx(np.zeros((3, 3)), abs=1e-6)
    # [Kaw] D = [4.935e4 8.148e2 0; 8.148e2 4.696e3 0; 0 0 2.017e3] Pa m^3
    assert D[0, 0] == pytest.approx(4.935e4, rel=1e-3)
    assert D[0, 1] == pytest.approx(8.148e2, rel=1e-3)
    assert D[1, 1] == pytest.approx(4.696e3, rel=1e-3)
    assert D[2, 2] == pytest.approx(2.017e3, rel=1e-3)


def test_kaw_090_engineering_constants():
    c = kaw_090().engineering_constants()
    # [Kaw] in-plane: Ex 124.5 GPa, Ey 67.43, Gxy 7.17, nuxy 0.04292, nuyx 0.02323
    assert c["Ex"] / 1e9 == pytest.approx(124.5, rel=1e-3)
    assert c["Ey"] / 1e9 == pytest.approx(67.43, rel=1e-3)
    assert c["Gxy"] / 1e9 == pytest.approx(7.17, rel=1e-3)
    assert c["nuxy"] == pytest.approx(0.04292, rel=1e-3)
    assert c["nuyx"] == pytest.approx(0.02323, rel=1e-3)
    # [Kaw] flexural: Ex 175.0 GPa, Ey 16.65, nuxy 0.1735, nuyx 0.01651
    assert c["Ex_flex"] / 1e9 == pytest.approx(175.0, rel=1e-3)
    assert c["Ey_flex"] / 1e9 == pytest.approx(16.65, rel=1e-3)
    assert c["nuxy_flex"] == pytest.approx(0.1735, rel=1e-3)
    assert c["nuyx_flex"] == pytest.approx(0.01651, rel=1e-3)


def test_kaw_090_ply_stresses_unit_load():
    # [Kaw] Nx = 1 N/m: eps0_x = 5.353e-10, eps0_y = -2.297e-11; local stresses (Pa)
    # 0 deg plies (97.26, 1.313, 0), 90 deg ply (-2.626, 5.472, 0), uniform through each ply.
    lam = kaw_090()
    eps0, kap = lam.response(N=(1, 0, 0))
    assert eps0 == pytest.approx([5.353e-10, -2.297e-11, 0], rel=1e-3, abs=1e-20)
    assert kap == pytest.approx([0, 0, 0], abs=1e-18)
    for pt in lam.ply_points(N=(1, 0, 0)):
        ref = [-2.626, 5.472, 0] if pt.ply == 1 else [97.26, 1.313, 0]
        assert pt.stress_12 == pytest.approx(ref, rel=1e-3, abs=1e-9)


def test_kaw_lamina_failure_criteria():
    # [Kaw] 60 deg lamina under sx = 2S, sy = -3S, txy = 4S: maximum S is
    # 16.33 MPa (maximum stress and maximum strain), 10.94 MPa (Tsai-Hill), 16.06 MPa
    # (modified Tsai-Hill), 22.39 MPa (Tsai-Wu with H12 = -1/2 sqrt(H11 H22)).
    s12 = stress_transformation(60 * DEG) @ np.array([2.0, -3.0, 4.0])
    e12 = reduced_compliance(GE) @ s12
    assert failure.max_stress(s12, GE)[0] / 1e6 == pytest.approx(16.33, abs=0.01)
    assert failure.max_stress(s12, GE)[1] == "12S"
    assert failure.max_strain(e12, GE)[0] / 1e6 == pytest.approx(16.33, abs=0.01)
    assert failure.tsai_hill(s12, GE) / 1e6 == pytest.approx(10.94, abs=0.01)
    assert failure.tsai_hill(s12, GE, modified=True) / 1e6 == pytest.approx(16.06, abs=0.01)
    assert failure.tsai_wu(s12, GE) / 1e6 == pytest.approx(22.39, abs=0.01)
    # the Tsai-Wu index equals 1 at the strength ratio
    assert failure.tsai_wu_index(s12 * failure.tsai_wu(s12, GE), GE) == pytest.approx(1.0)


def test_kaw_090_ply_by_ply_failure():
    # [Kaw] [0/90/0] under Nx: strength ratios at Nx = 1 N/m are 1.339e7 (0 deg) and 7.277e6
    # (90 deg) by Tsai-Wu; 1.548e7 (1T) and 7.254e6 (2T) by maximum strain. First ply failure
    # is the 90 deg ply at Nx = 7.277e6 N/m; with it discounted the 0 deg plies carry 100 Pa
    # per N/m and fail at Nx = 1.5e7 N/m.
    lam = kaw_090()
    by_ply = {pt.ply: sr for pt, sr in lam.strength_ratios(N=(1, 0, 0))}
    assert by_ply[0] == pytest.approx(1.339e7, rel=1e-3)
    assert by_ply[1] == pytest.approx(7.277e6, rel=1e-3)
    by_ply = {pt.ply: sr for pt, sr in lam.strength_ratios(N=(1, 0, 0), criterion="max_strain")}
    assert by_ply[0] == pytest.approx(1.548e7, rel=1e-3)
    assert by_ply[1] == pytest.approx(7.254e6, rel=1e-3)
    history = lam.ply_by_ply_failure(N=(1, 0, 0))
    assert history[0][0] == pytest.approx(7.277e6, rel=1e-3)
    assert history[0][1] == [1]
    assert history[1][0] == pytest.approx(1.5e7, rel=1e-6)
    assert history[1][1] == [0, 2]
    assert lam.active.all()  # progressive analysis leaves the laminate intact


def test_isotropic_stack_reproduces_plate_theory():
    E, nu, t = 70e9, 0.3, 0.001
    iso = OrthotropicPly("iso", E, E, E / (2 * (1 + nu)), nu)
    lam = Laminate.from_angles(iso, [0, 37, -80, 12], t)
    h = 4 * t
    A, B, D = lam.stiffness()
    assert A[0, 0] == pytest.approx(E * h / (1 - nu**2))
    assert A[0, 1] == pytest.approx(nu * E * h / (1 - nu**2))
    assert D[0, 0] == pytest.approx(E * h**3 / (12 * (1 - nu**2)))
    assert B == pytest.approx(np.zeros((3, 3)), abs=1e-3)


def test_symmetric_and_antisymmetric_coupling():
    t = 0.000125
    sym = Laminate.from_angles(GE, [0, 45, -45, 90], t, symmetric=True)
    assert sym.B == pytest.approx(np.zeros((3, 3)), abs=1e-6)
    anti = Laminate.from_angles(GE, [30, -30], t)
    A, B, _ = anti.stiffness()
    assert A[0, 2] == pytest.approx(0, abs=1e-3) and A[1, 2] == pytest.approx(0, abs=1e-3)
    assert abs(B[0, 2]) > 0 and abs(B[1, 2]) > 0
    assert B[0, 0] == pytest.approx(0, abs=1e-6) and B[2, 2] == pytest.approx(0, abs=1e-6)
    cross = Laminate.from_angles(GE, [0, 90], t)
    B = cross.B
    assert B[0, 0] == pytest.approx(-B[1, 1])
    assert B[0, 1] == pytest.approx(0, abs=1e-6) and B[2, 2] == pytest.approx(0, abs=1e-6)


def test_quasi_isotropic_in_plane():
    # [0/45/-45/90]s: A11 = A22 = U1 h, A12 = U4 h, A66 = U5 h = (A11 - A12) / 2, A16 = A26 = 0
    t = 0.000125
    lam = Laminate.from_angles(GE, [0, 45, -45, 90], t, symmetric=True)
    U, h, A = invariants(GE), lam.thickness, lam.A
    assert A[0, 0] == pytest.approx(U[0] * h)
    assert A[1, 1] == pytest.approx(U[0] * h)
    assert A[0, 1] == pytest.approx(U[3] * h)
    assert A[2, 2] == pytest.approx((A[0, 0] - A[0, 1]) / 2)
    assert A[0, 2] == pytest.approx(0, abs=1e-3)
    # rotating the laminate does not change its in-plane stiffness
    rot = Laminate([Ply(p.material, p.angle + 0.3, p.thickness) for p in lam.plies])
    assert rot.engineering_constants()["Ex"] == pytest.approx(lam.engineering_constants()["Ex"])


def test_unidirectional_laminate_constants():
    lam = Laminate.from_angles(GE, [0] * 8, 0.000125)
    c = lam.engineering_constants()
    assert c["Ex"] == pytest.approx(GE.E1)
    assert c["Ey"] == pytest.approx(GE.E2)
    assert c["nuxy"] == pytest.approx(GE.nu12)
    assert c["Gxy"] == pytest.approx(GE.G12)
    # uniaxial tension on UD: first-ply failure at Nx = Xt h for every criterion (Tsai-Wu
    # reduces to F1 s + F11 s^2 = 1, whose positive root is s = Xt exactly)
    for crit in ("tsai_wu", "tsai_hill", "max_stress", "max_strain"):
        sr, _ = lam.first_ply_failure(N=(1, 0, 0), criterion=crit)
        assert sr == pytest.approx(GE.Xt * lam.thickness, rel=1e-9)


def test_pure_bending_strain_linear_in_z():
    lam = Laminate.from_angles(GE, [0, 90, 90, 0], 0.000125)
    pts = lam.ply_points(M=(10, 0, 0))
    z = np.array([p.z for p in pts])
    ex = np.array([p.strain_xy[0] for p in pts])
    k = np.polyfit(z, ex, 1)
    assert np.allclose(np.polyval(k, z), ex, atol=1e-15)
    assert k[1] == pytest.approx(0, abs=1e-15)  # symmetric: no mid-plane strain
