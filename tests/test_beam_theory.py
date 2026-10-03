import numpy as np
import pytest

from structures.beam_theory import (
    Beam,
    Couple,
    DistributedLoad,
    PointLoad,
    Support,
    bredt_batho,
    circle,
    circular_shaft,
    closed_form,
    i_section,
    open_thin_walled,
    rectangle,
    thin_walled,
    tube,
)

L, EI = 2.0, 3.0e5  # m, N m^2
P, w = 1.2e3, 800.0  # N, N/m (magnitudes; applied downward below)


def reactions(sol):
    return [f for _, f, _ in sol.reactions], [m for _, _, m in sol.reactions]


def test_cantilever_tip_load():
    s = Beam(L, EI).add(Support(0, "fixed"), PointLoad(L, -P)).solve()
    ref = closed_form.cantilever_tip_load(P, L, EI)
    assert s.deflection(L) == pytest.approx(-ref["deflection"])
    assert s.slope(L) == pytest.approx(-ref["slope"])
    assert s.moment(0) == pytest.approx(-ref["moment"])  # hogging at the root
    F, M = reactions(s)
    assert F == pytest.approx([P])
    assert M == pytest.approx([P * L])  # CCW reaction couple


def test_cantilever_udl_fixed_at_right():
    # Mirror image: fixed at x = L, free at x = 0.
    s = Beam(L, EI).add(Support(L, "fixed"), DistributedLoad(0, L, -w)).solve()
    ref = closed_form.cantilever_udl(w, L, EI)
    assert s.deflection(0) == pytest.approx(-ref["deflection"])
    assert s.moment(L) == pytest.approx(-ref["moment"])


def test_simply_supported_udl():
    s = Beam(L, EI).add(Support(0), Support(L, "roller"), DistributedLoad(0, L, -w)).solve()
    ref = closed_form.simply_supported_udl(w, L, EI)
    assert s.deflection(L / 2) == pytest.approx(-ref["deflection"])
    assert s.moment(L / 2) == pytest.approx(ref["moment"])
    assert s.slope(0) == pytest.approx(-ref["slope"])
    assert s.shear(1e-9) == pytest.approx(w * L / 2)
    assert s.deflection(0) == pytest.approx(0, abs=1e-15)


@pytest.mark.parametrize("a", [1.4, 0.5])  # load right and left of midspan
def test_simply_supported_offset_point_load(a):
    s = Beam(L, EI).add(Support(0), Support(L), PointLoad(a, -P)).solve()
    ref = closed_form.simply_supported_point_load(P, a, L, EI)
    x, v = s.extreme("deflection", n=20001)
    assert v == pytest.approx(-ref["deflection"], rel=1e-6)
    assert x == pytest.approx(ref["x_max"], abs=2e-4)
    assert reactions(s)[0] == pytest.approx([ref["R_left"], ref["R_right"]])
    assert s.moment(a) == pytest.approx(ref["moment"])


def test_triangular_load():
    s = Beam(L, EI).add(Support(0), Support(L), DistributedLoad(0, L, 0.0, -w)).solve()
    ref = closed_form.simply_supported_triangular(w, L, EI)
    F, _ = reactions(s)
    assert F == pytest.approx([ref["R_left"], ref["R_right"]])
    x, v = s.extreme("deflection", n=20001)
    assert v == pytest.approx(-ref["deflection"], rel=1e-8)
    assert x == pytest.approx(ref["x_max"], abs=1e-3)
    assert s.extreme("moment", n=20001)[1] == pytest.approx(ref["moment"], rel=1e-6)


def test_propped_cantilever_udl_indeterminate():
    s = Beam(L, EI).add(Support(0, "fixed"), Support(L), DistributedLoad(0, L, -w)).solve()
    ref = closed_form.propped_cantilever_udl(w, L, EI)
    F, M = reactions(s)
    assert F == pytest.approx([ref["R_fixed"], ref["R_prop"]])
    assert s.moment(0) == pytest.approx(-ref["moment"])
    x, v = s.extreme("deflection", n=20001)
    assert v == pytest.approx(-ref["deflection"], rel=1e-8)
    assert x == pytest.approx(ref["x_max"], abs=1e-3)


def test_fixed_fixed():
    s = Beam(L, EI).add(Support(0, "fixed"), Support(L, "fixed"), DistributedLoad(0, L, -w))
    s = s.solve()
    ref = closed_form.fixed_fixed_udl(w, L, EI)
    assert s.moment(0) == pytest.approx(-ref["moment"])
    assert s.moment(L / 2) == pytest.approx(ref["moment_mid"])
    assert s.deflection(L / 2) == pytest.approx(-ref["deflection"])
    s = Beam(L, EI).add(Support(0, "fixed"), Support(L, "fixed"), PointLoad(L / 2, -P)).solve()
    ref = closed_form.fixed_fixed_center_load(P, L, EI)
    assert s.deflection(L / 2) == pytest.approx(-ref["deflection"])
    assert s.moment(L) == pytest.approx(-ref["moment"])


def test_two_span_continuous():
    s = Beam(2 * L, EI).add(Support(0), Support(L), Support(2 * L), DistributedLoad(0, 2 * L, -w))
    s = s.solve()
    ref = closed_form.two_span_continuous_udl(w, L, EI)
    F, _ = reactions(s)
    assert F == pytest.approx([ref["R_end"], ref["R_mid"], ref["R_end"]])
    assert s.moment(L) == pytest.approx(-ref["moment"])


def test_end_couple():
    M0 = 500.0
    s = Beam(L, EI).add(Support(0), Support(L), Couple(L, M0)).solve()
    ref = closed_form.simply_supported_end_couple(M0, L, EI)
    assert abs(s.slope(L)) == pytest.approx(ref["slope_loaded"])
    assert abs(s.slope(0)) == pytest.approx(ref["slope_far"])
    F, _ = reactions(s)
    assert F == pytest.approx([ref["R"], -ref["R"]])
    # A CCW couple at the right end bends the beam concave-up there: sagging M = M0 x / L
    assert s.moment(L / 2) == pytest.approx(M0 / 2)
    assert s.moment(L) == pytest.approx(M0)


def test_equilibrium_random_loading():
    rng = np.random.default_rng(3)
    b = Beam(5.0, EI).add(Support(0, "fixed"), Support(2.0), Support(5.0))
    for _ in range(6):
        b.add(PointLoad(rng.uniform(0, 5), rng.normal(0, 1e3)))
    b.add(DistributedLoad(1.0, 4.0, -300.0, 200.0), Couple(3.3, 750.0))
    s = b.solve()
    F, M = reactions(s)
    applied = sum(ld.P for ld in b.loads if isinstance(ld, PointLoad)) + (-300 + 200) / 2 * 3
    assert sum(F) == pytest.approx(-applied)
    # just left of the right-hand roller: V balances its reaction, M vanishes
    assert s.shear(5.0) == pytest.approx(-F[-1])
    assert s.moment(5.0) == pytest.approx(0, abs=1e-6)
    assert [float(s.deflection(x)) for x in (0, 2, 5)] == pytest.approx([0, 0, 0], abs=1e-12)


def test_mechanism_rejected():
    with pytest.raises(ValueError):
        Beam(L, EI).add(Support(0), PointLoad(1, -1)).solve()


def test_bending_stress_sign():
    # Sagging moment: compression on top (y > 0), tension at the bottom.
    sec = rectangle(0.05, 0.1)
    s = Beam(L, EI).add(Support(0), Support(L), PointLoad(L / 2, -P)).solve()
    top = s.bending_stress(L / 2, sec.y_top, sec.I)
    assert top == pytest.approx(-(P * L / 4) / sec.S_top)
    assert s.bending_stress(L / 2, -sec.y_bottom, sec.I) == pytest.approx(-top)


def test_section_properties():
    r = rectangle(0.04, 0.12)
    assert r.I == pytest.approx(0.04 * 0.12**3 / 12)
    # square of side a: exact Saint-Venant J = 0.1406 a^4 (Roark, Table 10.7 case 4)
    assert rectangle(1, 1).J == pytest.approx(0.1406, rel=2e-3)
    c = circle(0.05)
    assert c.J == pytest.approx(np.pi * 0.05**4 / 32)
    t = tube(0.1, 0.09)
    assert t.I == pytest.approx(np.pi * (0.1**4 - 0.09**4) / 64)


def test_thin_walled_matches_i_section():
    bf, tf, h, tw = 0.1, 0.008, 0.2, 0.005
    exact = i_section(bf, tf, h, tw)
    hm = h - tf  # flange mid-line separation
    segs = [
        ((-bf / 2, hm / 2), (bf / 2, hm / 2), tf),
        ((-bf / 2, -hm / 2), (bf / 2, -hm / 2), tf),
        ((0, -hm / 2 + tf / 2), (0, hm / 2 - tf / 2), tw),
    ]
    tw_sec = thin_walled(segs)
    assert tw_sec.A == pytest.approx(exact.A)
    assert tw_sec.I == pytest.approx(exact.I, rel=1e-12)
    assert tw_sec.I_product == pytest.approx(0, abs=1e-15)
    assert tw_sec.J == pytest.approx(exact.J, rel=1e-12)


def test_thin_walled_channel_centroid_and_product():
    # Unequal-flange angle has a non-zero product of inertia; check against direct integration.
    t = 0.002
    s = thin_walled([((0, 0), (0.06, 0), t), ((0, 0), (0, 0.04), t)])
    A1, A2 = 0.06 * t, 0.04 * t
    zc = (A1 * 0.03) / (A1 + A2)
    yc = (A2 * 0.02) / (A1 + A2)
    Izy = A1 * (0.03 - zc) * (0 - yc) + A2 * (0 - zc) * (0.02 - yc)
    assert s.I_product == pytest.approx(Izy, rel=1e-12)


def test_torsion_thin_tube_three_ways():
    # Thin circular tube r = 50 mm, t = 1 mm: polar J ~ 2 pi r^3 t from exact, Bredt-Batho and
    # the thin_walled builder (closed polygon approximating the circle).
    r, t, G, T = 0.05, 0.001, 27e9, 400.0
    exact = circular_shaft(T, G, 2 * r + t, 2 * r - t)
    bb = bredt_batho(T, G, np.pi * r**2, [(2 * np.pi * r, t)])
    assert bb.J == pytest.approx(exact.J, rel=1e-3)
    assert bb.tau_max == pytest.approx(T / (2 * np.pi * r**2 * t))
    n = 720
    pts = [(r * np.cos(a), r * np.sin(a)) for a in np.linspace(0, 2 * np.pi, n, endpoint=False)]
    segs = [(pts[i], pts[(i + 1) % n], t) for i in range(n)]
    assert thin_walled(segs, closed=True).J == pytest.approx(bb.J, rel=1e-4)
    # Open (slit) tube: J = 2 pi r t^3 / 3, i.e. (3/t^2) r^2 times weaker
    op = open_thin_walled(T, G, [(2 * np.pi * r, t)])
    assert op.J == pytest.approx(2 * np.pi * r * t**3 / 3)
    assert bb.J / op.J == pytest.approx(3 * r**2 / t**2, rel=1e-12)


def test_circular_shaft():
    d, T, G, Lsh = 0.03, 150.0, 80e9, 1.5
    res = circular_shaft(T, G, d)
    assert res.tau_max == pytest.approx(16 * T / (np.pi * d**3))
    assert res.twist(Lsh) == pytest.approx(T * Lsh / (G * np.pi * d**4 / 32))
