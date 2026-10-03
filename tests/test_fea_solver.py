import numpy as np
import pytest

from structures.beam_theory import closed_form
from structures.fea_solver import Model, beam_model, consistent_load, frame_stiffness, rotation

E, A, I = 70e9, 4e-4, 2e-6  # aluminium, 400 mm^2, 2e-6 m^4
P, w = 5e3, 2e3


def test_two_bar_truss():
    # Symmetric two-bar truss, apex load P down: N = -P / (2 sin a) in each bar,
    # apex deflection = P Lb / (2 E A sin^2 a) (virtual work / method of joints).
    span, h = 3.0, 1.0
    m = Model()
    a, b, c = m.add_node(0, 0), m.add_node(span, 0), m.add_node(span / 2, h)
    m.add_bar(a, c, E, A)
    m.add_bar(b, c, E, A)
    m.fix(a, "xy").fix(b, "xy").load(c, fy=-P)
    s = m.solve()
    Lb = np.hypot(span / 2, h)
    sin_a = h / Lb
    assert s.axial_force(0) == pytest.approx(-P / (2 * sin_a))
    assert s.axial_force(1) == pytest.approx(-P / (2 * sin_a))
    assert s.displacement(c)[1] == pytest.approx(-P * Lb / (2 * E * A * sin_a**2))
    assert s.displacement(c)[0] == pytest.approx(0, abs=1e-15)
    assert s.reaction(a)[1] + s.reaction(b)[1] == pytest.approx(P)


def test_three_bar_indeterminate_truss():
    # Vertical bar of length L flanked by two bars at angle th from vertical, all meeting at the
    # loaded node: N_mid = P / (1 + 2 cos^3 th), N_side = P cos^2 th / (1 + 2 cos^3 th)
    # (compatibility of the three elongations; Timoshenko, Strength of Materials, Part I).
    Lv, th = 1.5, np.radians(35)
    m = Model()
    tip = m.add_node(0, 0)
    for x in (-Lv * np.tan(th), 0.0, Lv * np.tan(th)):
        top = m.add_node(x, Lv)
        m.add_bar(top, tip, E, A)
        m.fix(top, "xy")
    m.load(tip, fy=-P)
    s = m.solve()
    d = 1 + 2 * np.cos(th) ** 3
    assert s.axial_force(1) == pytest.approx(P / d)
    assert s.axial_force(0) == pytest.approx(P * np.cos(th) ** 2 / d)
    assert s.displacement(tip)[1] == pytest.approx(-P / d * Lv / (E * A))


def test_determinate_truss_method_of_joints():
    # Two-panel truss A(0,0) B(a,0) C(2a,0) D(a,h), load P down at B, pin A, roller C:
    # N_BD = P, N_AD = N_CD = -P/(2 sin al), N_AB = N_BC = P/(2 tan al).
    a, h = 2.0, 1.5
    m = Model()
    A_, B, C, D = (m.add_node(*p) for p in [(0, 0), (a, 0), (2 * a, 0), (a, h)])
    bars = {
        k: m.add_bar(i, j, E, A)
        for k, (i, j) in {
            "AB": (A_, B),
            "BC": (B, C),
            "AD": (A_, D),
            "DC": (D, C),
            "BD": (B, D),
        }.items()
    }
    m.fix(A_, "xy").fix(C, "y").load(B, fy=-P)
    s = m.solve()
    al = np.arctan2(h, a)
    assert s.axial_force(bars["BD"]) == pytest.approx(P)
    assert s.axial_force(bars["AD"]) == pytest.approx(-P / (2 * np.sin(al)))
    assert s.axial_force(bars["DC"]) == pytest.approx(-P / (2 * np.sin(al)))
    assert s.axial_force(bars["AB"]) == pytest.approx(P / (2 * np.tan(al)))
    assert s.reaction(A_)[0] == pytest.approx(0, abs=1e-9)


@pytest.mark.parametrize("n", [1, 2, 5])
def test_cantilever_tip_load_nodally_exact(n):
    L = 2.0
    m = beam_model(L, n, E, A, I)
    m.fix(0).load(n, fy=-P)
    s = m.solve()
    ref = closed_form.cantilever_tip_load(P, L, E * I)
    assert s.displacement(n)[1] == pytest.approx(-ref["deflection"])
    assert s.displacement(n)[2] == pytest.approx(-ref["slope"])
    assert s.reaction(0)[2] == pytest.approx(P * L)


@pytest.mark.parametrize("n", [2, 4, 8])
def test_simply_supported_udl_consistent_loads_exact(n):
    # Hermite elements with consistent loads reproduce the exact nodal solution of
    # EI v'''' = q for any mesh (the exact solution of a uniform load is a quartic, but the
    # nodal values are interpolated exactly), so midspan is exact for even n.
    L = 3.0
    m = beam_model(L, n, E, A, I)
    m.fix(0, "xy").fix(n, "y")
    for e in range(n):
        m.distributed(e, -w)
    s = m.solve()
    ref = closed_form.simply_supported_udl(w, L, E * I)
    assert s.displacement(n // 2)[1] == pytest.approx(-ref["deflection"])
    assert s.displacement(0)[2] == pytest.approx(-ref["slope"])
    # internal moment from end forces + element load is exact everywhere
    assert s.moment(n // 2 - 1, 1.0) == pytest.approx(ref["moment"])
    assert s.moment(0, 0.5) == pytest.approx(w * (L / n / 2) * (L - L / n / 2) / 2)


def test_lumped_loads_converge_quadratically():
    L = 3.0
    ref = closed_form.simply_supported_udl(w, L, E * I)["deflection"]
    errs = []
    for n in (4, 8, 16, 32):
        m = beam_model(L, n, E, A, I)
        m.fix(0, "xy").fix(n, "y")
        for e in range(n):
            m.distributed(e, -w, lumped=True)
        errs.append(abs(-m.solve().displacement(n // 2)[1] - ref) / ref)
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert rates == pytest.approx([2, 2, 2], abs=0.05)


def test_propped_cantilever_and_fixed_fixed():
    L, n = 2.5, 6
    m = beam_model(L, n, E, A, I)
    m.fix(0).fix(n, "y")
    for e in range(n):
        m.distributed(e, -w)
    s = m.solve()
    ref = closed_form.propped_cantilever_udl(w, L, E * I)
    assert s.reaction(n)[1] == pytest.approx(ref["R_prop"])
    assert s.reaction(0)[2] == pytest.approx(ref["moment"])
    m = beam_model(L, 4, E, A, I)
    m.fix(0).fix(4)
    m.load(2, fy=-P)
    ref = closed_form.fixed_fixed_center_load(P, L, E * I)
    assert m.solve().displacement(2)[1] == pytest.approx(-ref["deflection"])


def test_linearly_varying_load():
    # Triangular load on a simply supported beam: reactions w0L/6 and w0L/3.
    L, n = 2.0, 4
    m = beam_model(L, n, E, A, I)
    m.fix(0, "xy").fix(n, "y")
    for e in range(n):
        m.distributed(e, -w * e / n, -w * (e + 1) / n)
    s = m.solve()
    ref = closed_form.simply_supported_triangular(w, L, E * I)
    assert s.reaction(0)[1] == pytest.approx(ref["R_left"])
    assert s.reaction(n)[1] == pytest.approx(ref["R_right"])


def test_l_frame_rigid_knee():
    # Column (height h, fixed base) + beam (length b), tip load P down. Tip deflection:
    # P b^3 / 3EI (beam) + P b^2 h / EI (column rotation) + P h / EA (column shortening);
    # the knee sways towards the load by P b h^2 / 2EI (Castigliano / unit-load method).
    h, b = 2.0, 1.2
    m = Model()
    n0, n1, n2 = m.add_node(0, 0), m.add_node(0, h), m.add_node(b, h)
    m.add_frame(n0, n1, E, A, I)
    m.add_frame(n1, n2, E, A, I)
    m.fix(n0).load(n2, fy=-P)
    s = m.solve()
    EI = E * I
    expected = P * b**3 / (3 * EI) + P * b**2 * h / EI + P * h / (E * A)
    assert s.displacement(n2)[1] == pytest.approx(-expected)
    assert s.displacement(n1)[0] == pytest.approx(P * b * h**2 / (2 * EI))


def test_rotated_model_is_invariant():
    # The same cantilever built at 37 degrees must give the same displacement magnitude.
    L, n, ang = 2.0, 3, np.radians(37)
    flat, tilted = beam_model(L, n, E, A, I), beam_model(L, n, E, A, I, angle=ang)
    flat.fix(0).load(n, fy=-P)
    tilted.fix(0).load(n, fx=P * np.sin(ang), fy=-P * np.cos(ang))
    u1, u2 = flat.solve().displacement(n), tilted.solve().displacement(n)
    assert np.hypot(*u2[:2]) == pytest.approx(np.hypot(*u1[:2]))
    assert u2[2] == pytest.approx(u1[2])


def test_element_stiffness_properties():
    k = frame_stiffness(E, A, I, 1.3)
    assert np.allclose(k, k.T)
    eig = np.linalg.eigvalsh(k)
    assert np.sum(np.abs(eig) < 1e-6 * eig.max()) == 3  # three rigid-body modes
    T = rotation(np.cos(0.4), np.sin(0.4))
    assert T @ T.T == pytest.approx(np.eye(6))
    f = consistent_load(2.0, 1.0, 1.0)
    assert f == pytest.approx([0, 1, 1 / 3, 0, 1, -1 / 3])


def test_mechanism_raises():
    m = Model()
    a, b = m.add_node(0, 0), m.add_node(1, 0)
    m.add_bar(a, b, E, A)
    m.fix(a, "x").load(b, fy=-1)
    with pytest.raises(np.linalg.LinAlgError):
        m.solve()


def test_prescribed_settlement():
    # Propped cantilever with the prop settled by d: prop reaction 3 EI d / L^3.
    L, d = 2.0, 1e-3
    m = beam_model(L, 2, E, A, I)
    m.fix(0).fix(2, "y", value=-d)
    assert m.solve().reaction(2)[1] == pytest.approx(-3 * E * I * d / L**3)
