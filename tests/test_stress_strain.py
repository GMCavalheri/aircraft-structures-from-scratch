import numpy as np
import pytest

from structures.stress_strain import (
    MohrCircle,
    compliance_matrix,
    max_shear_2d,
    mohr_circles_3d,
    plane_strain_stiffness,
    plane_stress_stiffness,
    principal_2d,
    principal_3d,
    safety_factor,
    stiffness_matrix,
    stress_invariants,
    stress_tensor,
    transform_strain_2d,
    transform_stress_2d,
    transform_tensor,
    tresca,
    von_mises,
)

DEG = np.pi / 180


def test_transform_90_degrees_swaps_normals():
    s = [80.0, -30.0, 25.0]
    sx, sy, txy = transform_stress_2d(s, 90 * DEG)
    assert (sx, sy, txy) == pytest.approx((-30, 80, -25))


def test_transform_matches_tensor_rotation():
    # The closed-form plane transformation equals R sigma R^T with R about z.
    s = [50e6, -20e6, 35e6]
    t = 27 * DEG
    c, n = np.cos(t), np.sin(t)
    R = np.array([[c, n, 0], [-n, c, 0], [0, 0, 1]])
    T = transform_tensor(stress_tensor(sx=s[0], sy=s[1], txy=s[2]), R)
    assert transform_stress_2d(s, t) == pytest.approx([T[0, 0], T[1, 1], T[0, 1]])


def test_invariants_preserved_under_rotation():
    rng = np.random.default_rng(0)
    sigma = stress_tensor(*rng.normal(size=6))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    assert stress_invariants(transform_tensor(sigma, q)) == pytest.approx(stress_invariants(sigma))


def test_pure_shear_principal_at_45_degrees():
    tau = 40e6
    s1, s2, tp = principal_2d([0, 0, tau])
    assert (s1, s2) == pytest.approx((tau, -tau))
    assert tp == pytest.approx(45 * DEG)
    assert transform_stress_2d([0, 0, tau], tp) == pytest.approx([tau, -tau, 0], abs=1e-6)


def test_principal_hand_calculation():
    # sx = -20, sy = 90, txy = 60 MPa: centre 35, R = sqrt(55^2 + 60^2) = 81.39 MPa,
    # s1 = 116.39, s2 = -46.39, theta_p = atan2(120, -110)/2 = 132.51/2 = 66.26 deg (s1 direction).
    s1, s2, tp = principal_2d([-20, 90, 60])
    assert s1 == pytest.approx(35 + np.hypot(55, 60))
    assert s2 == pytest.approx(35 - np.hypot(55, 60))
    assert tp / DEG == pytest.approx(66.26, abs=0.01)
    # shear vanishes on principal planes and is extreme 45 deg away
    assert transform_stress_2d([-20, 90, 60], tp)[2] == pytest.approx(0, abs=1e-12)
    tmax, ts = max_shear_2d([-20, 90, 60])
    assert abs(transform_stress_2d([-20, 90, 60], ts)[2]) == pytest.approx(tmax)


def test_principal_3d_matches_2d_and_sorted():
    p, v = principal_3d([30.0, -10.0, 20.0])
    s1, s2, _ = principal_2d([30.0, -10.0, 20.0])
    assert sorted([s1, s2, 0.0], reverse=True) == pytest.approx(p)
    assert np.allclose(v.T @ v, np.eye(3))


def test_von_mises_special_cases():
    s = 250e6
    assert von_mises([s, 0, 0]) == pytest.approx(s)  # uniaxial
    assert von_mises([0, 0, s]) == pytest.approx(np.sqrt(3) * s)  # pure shear
    assert von_mises(stress_tensor(s, s, s)) == pytest.approx(0, abs=1e-6)  # hydrostatic
    assert von_mises([s, s, 0]) == pytest.approx(s)  # equibiaxial
    # principal-stress form agrees with the component form
    sig = stress_tensor(10, -4, 7, 3, -2, 5)
    p, _ = principal_3d(sig)
    ref = np.sqrt(((p[0] - p[1]) ** 2 + (p[1] - p[2]) ** 2 + (p[2] - p[0]) ** 2) / 2)
    assert von_mises(sig) == pytest.approx(ref)


def test_tresca_and_safety_factor():
    s = 100e6
    assert tresca([0, 0, s]) == pytest.approx(2 * s)  # pure shear: s1 - s3 = 2 tau
    assert tresca([s, s, 0]) == pytest.approx(s)  # equibiaxial: out-of-plane zero governs
    # Tresca is never less than von Mises, and at most 2/sqrt(3) larger.
    rng = np.random.default_rng(1)
    for _ in range(50):
        sig = stress_tensor(*rng.normal(size=6))
        assert von_mises(sig) <= tresca(sig) + 1e-12 <= 2 / np.sqrt(3) * von_mises(sig) + 1e-9
    assert safety_factor([s, 0, 0], 300e6) == pytest.approx(3.0)
    assert safety_factor([0, 0, s], 300e6, "tresca") == pytest.approx(1.5)


def test_strain_transformation_engineering_shear():
    # Pure shear strain gxy rotates into principal strains +-gxy/2 at 45 deg.
    g = 1e-3
    assert transform_strain_2d([0, 0, g], 45 * DEG) == pytest.approx([g / 2, -g / 2, 0], abs=1e-15)


def test_hooke_uniaxial_and_inverse():
    E, nu = 70e9, 0.33
    S, C = compliance_matrix(E, nu), stiffness_matrix(E, nu)
    assert S @ C == pytest.approx(np.eye(6), abs=1e-12)
    eps = S @ np.array([100e6, 0, 0, 0, 0, 0])
    assert eps[:3] == pytest.approx([100e6 / E, -nu * 100e6 / E, -nu * 100e6 / E])
    gamma = S @ np.array([0, 0, 0, 0, 0, 50e6])
    assert gamma[5] == pytest.approx(50e6 * 2 * (1 + nu) / E)


def test_plane_stress_and_strain_reduce_from_3d():
    E, nu = 200e9, 0.3
    S = compliance_matrix(E, nu)
    idx = [0, 1, 5]
    assert np.linalg.inv(S[np.ix_(idx, idx)]) == pytest.approx(plane_stress_stiffness(E, nu))
    C = stiffness_matrix(E, nu)
    assert C[np.ix_(idx, idx)] == pytest.approx(plane_strain_stiffness(E, nu))


def test_mohr_circle_geometry():
    c = MohrCircle(-20.0, 90.0, 60.0)
    assert c.principal == pytest.approx(principal_2d([-20, 90, 60])[:2])
    for t in np.linspace(0, np.pi, 7):
        x, y = c.point(t)
        assert np.hypot(x - c.center, y) == pytest.approx(c.radius)
    circles = mohr_circles_3d([-20.0, 90.0, 60.0])
    assert max(r for _, r in circles) == pytest.approx(tresca([-20.0, 90.0, 60.0]) / 2)
