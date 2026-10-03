"""Generate and execute the demo notebooks from the cell sources below.

Run with ``uv run python notebooks/build_notebooks.py``. Keeping the sources here (instead of
hand-editing .ipynb JSON) makes the notebooks reviewable in diffs and reproducible.
"""

from __future__ import annotations

from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent

SETUP = """\
import numpy as np
import matplotlib.pyplot as plt
from structures import plotting as P

P.use_style()
DEG = np.pi / 180
MPA = 1e6"""

NOTEBOOKS: dict[str, list[tuple[str, str]]] = {
    "01_stress_strain.ipynb": [
        (
            "md",
            """\
# 1 · Stress and strain at a point

Rotate a plane stress state, find its principal stresses, draw Mohr's circle, and check yield
against a MIL-HDBK-5J allowable. Theory: [docs/stress_strain.md](../docs/stress_strain.md).""",
        ),
        (
            "code",
            SETUP
            + """
from structures.stress_strain import (MohrCircle, plot_mohr, principal_2d, max_shear_2d,
                                      transform_stress_2d, von_mises, tresca, safety_factor)
from structures.materials import isotropic, KSI""",
        ),
        (
            "code",
            """\
s = np.array([-20, 90, 60]) * MPA          # sx, sy, txy
s1, s2, tp = principal_2d(s)
tmax, ts = max_shear_2d(s)
print(f"s1 = {s1/MPA:.2f} MPa, s2 = {s2/MPA:.2f} MPa at theta_p = {tp/DEG:.2f} deg")
print(f"max in-plane shear = {tmax/MPA:.2f} MPa at {ts/DEG:.2f} deg")
print("check: shear on the principal plane =", transform_stress_2d(s, tp)[2])""",
        ),
        (
            "code",
            """\
theta = np.linspace(0, np.pi, 181)
sx, sy, txy = transform_stress_2d(s, theta)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
a1.plot(theta / DEG, sx / MPA, label="σx'"); a1.plot(theta / DEG, txy / MPA, label="τx'y'")
a1.axvline(tp / DEG, color=P.GRID); a1.set_xlabel("rotation θ [deg]"); a1.set_ylabel("MPa"); a1.legend()
plot_mohr(a2, MohrCircle(*s), scale=MPA, unit="MPa");""",
        ),
        (
            "md",
            "As the element rotates through 180°, the point on Mohr's circle goes once around: "
            "a rotation θ of the element is 2θ on the circle.",
        ),
        (
            "code",
            """\
m = isotropic("2024-T3", "B")
state = np.array([180, 60, 70]) * MPA
print(f"{m.name}: Fty = {m.Fty/KSI:.0f} ksi = {m.Fty/MPA:.0f} MPa ({m.source})")
print(f"von Mises = {von_mises(state)/MPA:.1f} MPa -> SF = {safety_factor(state, m.Fty):.3f}")
print(f"Tresca    = {tresca(state)/MPA:.1f} MPa -> SF = {safety_factor(state, m.Fty, 'tresca'):.3f}")""",
        ),
    ],
    "02_beam_theory.ipynb": [
        (
            "md",
            """\
# 2 · Beams and torsion

The Macaulay solver handles determinate and indeterminate beams with the same code. Theory:
[docs/beam_theory.md](../docs/beam_theory.md).""",
        ),
        (
            "code",
            SETUP
            + """
from structures.beam_theory import (Beam, Support, PointLoad, DistributedLoad, Couple,
                                    closed_form, i_section, thin_walled, bredt_batho,
                                    open_thin_walled)""",
        ),
        (
            "code",
            """\
L, EI, w = 2.0, 3.0e5, 800.0
s = Beam(L, EI).add(Support(0, "fixed"), Support(L), DistributedLoad(0, L, -w)).solve()
ref = closed_form.propped_cantilever_udl(w, L, EI)
print("reactions:", [(sup.kind, round(F, 3), round(M, 3)) for sup, F, M in s.reactions])
print(f"prop reaction {s.reactions[1][1]:.1f} N vs 3wL/8 = {ref['R_prop']:.1f} N")
x, v = s.extreme("deflection")
print(f"max deflection {-v*1e3:.4f} mm (down) at x = {x:.4f} m (closed form {ref['deflection']*1e3:.4f} mm at {ref['x_max']:.4f} m)")""",
        ),
        (
            "code",
            """\
beam = Beam(4.0, 2.0e6).add(Support(0), Support(3.0), PointLoad(1.0, -5e3),
                            DistributedLoad(1.5, 4.0, -2e3, -1e3), Couple(2.0, 3e3))
sol = beam.solve()
x = np.linspace(0, 4, 801)
fig, axes = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
for ax, f, lab in zip(axes, [sol.shear(x)/1e3, sol.moment(x)/1e3, sol.deflection(x)*1e3],
                      ["V [kN]", "M [kN m]", "v [mm]"]):
    ax.plot(x, f); ax.axhline(0, color=P.TEXT_MUTED, lw=0.8); ax.set_ylabel(lab)
axes[-1].set_xlabel("x [m]");""",
        ),
        (
            "md",
            "The couple at x = 2 m makes the bending moment jump, and the overhang beyond the "
            "support at 3 m hogs. Now torsion: a closed tube against the same tube slit open.",
        ),
        (
            "code",
            """\
r, t, G, T = 0.05, 0.001, 27e9, 400.0
closed = bredt_batho(T, G, np.pi * r**2, [(2 * np.pi * r, t)])
opened = open_thin_walled(T, G, [(2 * np.pi * r, t)])
print(f"closed: tau = {closed.tau_max/MPA:.2f} MPa, twist = {closed.twist(1.0)/DEG:.4f} deg/m")
print(f"open:   tau = {opened.tau_max/MPA:.1f} MPa, twist = {opened.twist(1.0)/DEG:.1f} deg/m"
      "  <- far beyond yield: the slit tube cannot carry this torque")
print(f"stiffness ratio = {closed.J/opened.J:.0f} (3 r^2 / t^2 = {3*r**2/t**2:.0f})")""",
        ),
    ],
    "03_fea_solver.ipynb": [
        (
            "md",
            """\
# 3 · The finite element solver

Bars and Euler–Bernoulli frames, sparse assembly, and a convergence study. Theory:
[docs/fea_solver.md](../docs/fea_solver.md).""",
        ),
        (
            "code",
            SETUP
            + """
from structures.fea_solver import Model, beam_model
from structures.beam_theory import closed_form""",
        ),
        (
            "code",
            """\
# Two-panel truss: A(0,0) B(2,0) C(4,0) D(2,1.5), 10 kN down at B
m = Model()
A, B, C, D = (m.add_node(*p) for p in [(0, 0), (2, 0), (4, 0), (2, 1.5)])
names = {}
for k, (i, j) in {"AB": (A, B), "BC": (B, C), "AD": (A, D), "DC": (D, C), "BD": (B, D)}.items():
    names[k] = m.add_bar(i, j, 200e9, 1e-3)
m.fix(A, "xy").fix(C, "y").load(B, fy=-10e3)
s = m.solve()
for k, e in names.items():
    print(f"{k}: {s.axial_force(e)/1e3:+8.3f} kN")
print("reactions (kN):", (s.reaction(A)[:2] / 1e3).round(3) + 0, (s.reaction(C)[:2] / 1e3).round(3) + 0)""",
        ),
        (
            "code",
            """\
L, E, A_, I = 3.0, 70e9, 4e-4, 2e-6
ref = closed_form.simply_supported_udl(2e3, L, E * I)["deflection"]
for n in (2, 4, 8):
    for lumped in (False, True):
        mdl = beam_model(L, n, E, A_, I)
        mdl.fix(0, "xy").fix(n, "y")
        for e in range(n):
            mdl.distributed(e, -2e3, lumped=lumped)
        v = -mdl.solve().displacement(n // 2)[1]
        print(f"n = {n}, {'lumped    ' if lumped else 'consistent'}: error {abs(v - ref)/ref:.2e}")""",
        ),
        (
            "md",
            "Consistent loads give exact nodal deflections on every mesh (Tong's theorem); "
            "lumped loads converge as h².",
        ),
        (
            "code",
            """\
# Portal frame: deflected shape
m = Model()
n = [m.add_node(0, 0), m.add_node(0, 3), m.add_node(4, 3), m.add_node(4, 0)]
for i in range(3):
    m.add_frame(n[i], n[i + 1], 200e9, 5e-3, 8e-5)
m.fix(n[0]).fix(n[3]).distributed(1, -20e3).load(n[1], fx=40e3)
s = m.solve()
fig, ax = plt.subplots(figsize=(5, 4))
for e in range(3):
    d0 = s.deflected_shape(e, scale=0); d = s.deflected_shape(e, scale=50)
    ax.plot(*d0.T, color=P.TEXT_MUTED, ls="--", lw=0.8); ax.plot(*d.T, color=P.SERIES[0])
ax.set_aspect("equal"); ax.set_title("portal frame, displacements x50");""",
        ),
    ],
    "04_buckling.ipynb": [
        (
            "md",
            """\
# 4 · Buckling

Euler columns, inelastic column curves, plate buckling, and FEA eigenvalue buckling. Theory:
[docs/buckling.md](../docs/buckling.md).""",
        ),
        (
            "code",
            SETUP
            + """
from structures.buckling import (EFFECTIVE_LENGTH, euler_load, column_model, linear_buckling,
                                 johnson_stress, tangent_modulus_stress, compression_k,
                                 plate_critical_stress)
from structures.materials import isotropic, KSI""",
        ),
        (
            "code",
            """\
E, A, I, L = 70e9, 3e-4, 1e-8, 1.2
for end, K in EFFECTIVE_LENGTH.items():
    P_e = euler_load(E, I, L, end)
    P_fe = linear_buckling(column_model(L, 8, E, A, I, end)).load_factors[0]
    print(f"{end:14s} K = {K:.4f}  Euler {P_e:9.2f} N   FEA(8 el.) {P_fe:9.2f} N   error {P_fe/P_e - 1:.1e}")""",
        ),
        (
            "code",
            """\
res = linear_buckling(column_model(L, 40, E, A, I, "pinned-pinned"), n_modes=3)
x = np.linspace(0, 1, 41)
fig, ax = plt.subplots(figsize=(6, 3.5))
for j in range(3):
    ax.plot(x, res.modes[1::3, j], label=f"mode {j+1}: P/P1 = {res.load_factors[j]/res.load_factors[0]:.3f}")
ax.set_xlabel("x / L"); ax.legend();""",
        ),
        (
            "code",
            """\
m = isotropic("2024-T3")
lam = np.array([10, 30, 50, 70, 90, 120])
print("L'/rho  Johnson[ksi]  tangent-modulus(n=15)[ksi]")
for l, j, t in zip(lam, johnson_stress(m.Ec, m.Fcy, lam), tangent_modulus_stress(m.Ec, m.Fcy, 15, lam)):
    print(f"{l:6.0f}  {j/KSI:12.1f}  {min(t, m.Fcy)/KSI:12.1f}")""",
        ),
        (
            "code",
            """\
# A 2024-T3 skin bay between stringers: 150 mm wide, 450 mm long, 1.6 mm thick
k, half_waves = compression_k(450 / 150)
s_cr = plate_critical_stress(m.Ec, m.nu, 1.6e-3, 0.150, k)
print(f"k = {k:.3f} with {half_waves} half-waves; sigma_cr = {s_cr/MPA:.1f} MPa "
      f"({s_cr/m.Fcy:.0%} of Fcy): thin skins buckle long before they yield")""",
        ),
    ],
    "05_composite_laminates.ipynb": [
        (
            "md",
            """\
# 5 · Composite laminates

Classical Lamination Theory: ABD matrices, ply stresses, and failure. Theory:
[docs/composite_laminates.md](../docs/composite_laminates.md).""",
        ),
        (
            "code",
            SETUP
            + """
from structures.composite_laminates import Laminate, reduced_stiffness
from structures.materials import ply""",
        ),
        (
            "code",
            """\
ge = ply("T300/5208")
print("Q [GPa] =\\n", np.round(reduced_stiffness(ge) / 1e9, 3))
lam = Laminate.from_angles(ge, [0, 90, 0], 0.005)   # Kaw's [0/90/0] example, 5 mm plies
A, B, D = lam.stiffness()
print("A [Pa m] =\\n", A.round(-4)); print("max |B| =", np.abs(B).max())
c = lam.engineering_constants()
print(f"Ex = {c['Ex']/1e9:.1f} GPa, Ey = {c['Ey']/1e9:.2f} GPa, nuxy = {c['nuxy']:.5f}")""",
        ),
        (
            "code",
            """\
for pt in lam.ply_points(N=(1, 0, 0)):
    if pt.position == "middle":
        print(f"ply {pt.ply} ({round(lam.plies[pt.ply].angle/DEG)} deg): sigma_12 = {pt.stress_12.round(3)} Pa per N/m")
for load, failed in lam.ply_by_ply_failure(N=(1, 0, 0)):
    print(f"Nx = {load:.4e} N/m: plies {failed} fail")""",
        ),
        (
            "md",
            "The 90° ply cracks first, at Nx = 7.28 MN/m. The 0° plies then carry everything "
            "until fibre failure at 15 MN/m, as in Kaw's ply-by-ply example.",
        ),
        (
            "code",
            """\
ht = ply("T300/976")
t = 0.0053 * 0.0254
qi = Laminate.from_angles(ht, [0, 45, -45, 90], t, symmetric=True)
phi = np.linspace(0, 2 * np.pi, 181)
fig, ax = plt.subplots(figsize=(5, 5))
for k, crit in enumerate(["tsai_wu", "max_stress"]):
    r = np.array([qi.first_ply_failure(N=(np.cos(p), np.sin(p), 0), criterion=crit)[0] for p in phi])
    ax.plot(r * np.cos(phi) / 1e3, r * np.sin(phi) / 1e3, label=crit, color=P.SERIES[k])
ax.set_aspect("equal"); ax.set_xlabel("Nx [kN/m]"); ax.set_ylabel("Ny [kN/m]"); ax.legend();""",
        ),
    ],
    "06_fatigue.ipynb": [
        (
            "md",
            """\
# 6 · Fatigue

MIL-HDBK-5J S-N curves, rainflow counting, and Miner's rule on a synthetic flight spectrum.
Theory: [docs/fatigue.md](../docs/fatigue.md).""",
        ),
        (
            "code",
            SETUP
            + """
from structures.fatigue import (handbook_curve, count_cycles, histogram, flight_history,
                                miner_damage, blocks_to_failure)
from structures.materials import KSI""",
        ),
        (
            "code",
            """\
astm = [-2, 1, -3, 5, -1, 3, -4, 4, -2]          # ASTM E1049 example
for c in count_cycles(astm):
    print(f"range {c.range:g}, mean {c.mean:+.1f}, count {c.count}")
print("histogram:", histogram(count_cycles(astm)))""",
        ),
        (
            "code",
            """\
c = handbook_curve("2024-T3_Kt2")
print(c.name, c.source, f"A1..A4 = {c.A1}, {c.A2}, {c.A3}, {c.A4} ksi")
for R in (-1, 0, 0.5):
    print(f"R = {R:+.1f}: Smax for 1e5 cycles = {c.max_stress_for_life(1e5, R)/KSI:.1f} ksi")""",
        ),
        (
            "code",
            """\
hist = flight_history(500, s_1g=12 * KSI, s_ground=-4 * KSI, beta=2.5 * KSI, seed=11)
cycles = count_cycles(hist)
D = miner_damage(cycles, c.life_from_amplitude)
print(f"{len(cycles)} cycles, damage per 500 flights = {D:.3e}, "
      f"median life = {blocks_to_failure(D) * 500:,.0f} flights")
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(hist[:250] / KSI, lw=0.8); ax.set_ylabel("stress [ksi]"); ax.set_xlabel("reversal");""",
        ),
        (
            "md",
            "That is a median life with no scatter factor. Certification divides it by a scatter "
            "factor and backs it up with crack-growth-based inspections.",
        ),
    ],
}


def _merge_streams(outputs):
    """Join consecutive stdout/stderr chunks; how the kernel splits them depends on timing."""
    merged = []
    for out in outputs:
        prev = merged[-1] if merged else None
        if (
            prev is not None
            and out.get("output_type") == "stream"
            and prev.get("output_type") == "stream"
            and prev.get("name") == out.get("name")
        ):
            prev["text"] += out["text"]
        else:
            merged.append(out)
    return merged


def build(name: str, cells: list[tuple[str, str]]) -> Path:
    nb = nbformat.v4.new_notebook()
    nb.metadata["kernelspec"] = {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    }
    nb.cells = [
        nbformat.v4.new_markdown_cell(src) if kind == "md" else nbformat.v4.new_code_cell(src)
        for kind, src in cells
    ]
    for i, cell in enumerate(nb.cells):
        cell.id = f"cell-{i:02d}"  # stable ids keep rebuilds diff-free
    NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        record_timing=False,
        resources={"metadata": {"path": str(HERE)}},
    ).execute()
    for cell in nb.cells:
        cell.get("outputs", [])[:] = _merge_streams(cell.get("outputs", []))
    path = HERE / name
    nbformat.write(nb, path)
    return path


def main() -> None:
    for name, cells in NOTEBOOKS.items():
        print("built", build(name, cells).name)


if __name__ == "__main__":
    main()
