import pytest

from structures.materials import KSI, MSI, isotropic, ply


def test_ksi_conversion():
    # 1 ksi = 1000 lbf / in^2 = 4448.2216 N / 6.4516e-4 m^2 (exact lbf and inch definitions)
    assert pytest.approx(4448.2216152605 / 6.4516e-4, rel=1e-6) == KSI
    assert MSI == pytest.approx(1000 * KSI)


def test_2024_t3_sheet_mil_hdbk_5j():
    # MIL-HDBK-5J Table 3.2.3.0(b1), T3 sheet 0.010-0.128 in, L direction.
    a, b = isotropic("2024-T3", "A"), isotropic("2024-T3", "B")
    assert a.Ftu / KSI == pytest.approx(64)
    assert b.Ftu / KSI == pytest.approx(65)
    assert b.Fty / KSI == pytest.approx(48)
    assert b.E / MSI == pytest.approx(10.5)
    assert b.nu == 0.33
    assert b.density == pytest.approx(2768, rel=1e-3)  # 0.100 lb/in^3


def test_7075_t6_and_others():
    # MIL-HDBK-5J Tables 3.7.6.0(b1), 5.4.1.0(b), 2.3.1.0(f1).
    assert isotropic("7075-T6").Ftu / KSI == pytest.approx(80)
    assert isotropic("Ti-6Al-4V", "A").Fty / KSI == pytest.approx(126)
    steel = isotropic("4340")  # S-basis only
    assert steel.basis == "S"
    assert steel.Ftu / KSI == pytest.approx(260)


def test_shear_modulus_close_to_isotropic_relation():
    # Tabulated G should be consistent with G = E / 2(1 + nu) to the handbook's rounding.
    for key in ["2024-T3", "7075-T6", "Ti-6Al-4V", "4340"]:
        m = isotropic(key)
        assert m.G == pytest.approx(m.E / (2 * (1 + m.nu)), rel=0.03)


def test_plies():
    t = ply("T300/5208")
    assert t.E1 == pytest.approx(181e9)
    assert t.nu21 == pytest.approx(0.28 * 10.3 / 181)
    assert t.Yt == pytest.approx(40e6)
    # MIL-HDBK-17-2F Table 4.2.17(a): F1tu = 211 ksi, E1t = 19.6 Msi.
    h = ply("T300/976")
    assert h.Xt / KSI == pytest.approx(211)
    assert h.E1 / MSI == pytest.approx(19.6)


def test_unknown_key():
    with pytest.raises(KeyError):
        isotropic("6061-T6")
    with pytest.raises(KeyError):
        ply("IM7/8552")
