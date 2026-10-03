# Beam theory validation

| Case | Computed | Closed form | Formula | Rel. error |
|---|---|---|---|---|
| Cantilever, tip load — tip deflection | 0.0106667 | 0.0106667 | PL³/3EI | +0.0e+00 |
| Cantilever, UDL — tip deflection | 0.00533333 | 0.00533333 | wL⁴/8EI | +0.0e+00 |
| Simply supported, centre load — midspan deflection | 0.000666667 | 0.000666667 | PL³/48EI | +0.0e+00 |
| Simply supported, UDL — midspan deflection | 0.000555556 | 0.000555556 | 5wL⁴/384EI | -5.9e-16 |
| Simply supported, triangular load — max deflection | 0.00027828 | 0.00027828 | 0.006522 w₀L⁴/EI | -2.0e-09 |
| Simply supported, end couple — far-end slope | 0.000555556 | 0.000555556 | M₀L/6EI | +0.0e+00 |
| Propped cantilever, UDL — prop reaction (indeterminate) | 600 | 600 | 3wL/8 | +1.9e-16 |
| Fixed-fixed, UDL — end moment (indeterminate) | 266.667 | 266.667 | wL²/12 | -4.3e-16 |
| Two-span continuous, UDL — middle reaction (indeterminate) | 2000 | 2000 | 10wL/8 | +0.0e+00 |
| Thin tube r = 50 mm, t = 1 mm: J, Bredt–Batho vs exact annulus | 7.8540e-07 m⁴ | 7.8548e-07 m⁴ | π(dₒ⁴−dᵢ⁴)/32 | -1.0e-04 |
| Same tube slit open: stiffness ratio closed/open | 7500 | 7500 | 3r²/t² | -2.4e-16 |

EI = 3×10⁵ N·m², L = 2 m, P = 1.2 kN, w = 800 N/m; SI units.
