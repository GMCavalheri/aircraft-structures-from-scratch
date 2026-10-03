# Composite laminates: CLT vs Kaw's worked examples

| Quantity | Computed | Kaw | Rel. error |
|---|---|---|---|
| [0/90/0] A11 [Pa·m] | 1.87e+09 | 1.87e+09 | -8.4e-05 |
| [0/90/0] A22 [Pa·m] | 1.013e+09 | 1.013e+09 | -4.8e-04 |
| [0/90/0] D11 [Pa·m³] | 4.935e+04 | 4.935e+04 | -3.5e-05 |
| [0/90/0] D22 [Pa·m³] | 4696 | 4696 | -1.1e-05 |
| [0/90/0] in-plane Ex [GPa] | 124.5 | 124.5 | +2.6e-04 |
| [0/90/0] in-plane Ey [GPa] | 67.43 | 67.43 | +5.7e-05 |
| [0/90/0] in-plane νxy | 0.04292 | 0.04292 | -7.8e-05 |
| [0/90/0] flexural Ex [GPa] | 175 | 175 | -2.4e-04 |
| [0/90/0] flexural Ey [GPa] | 16.65 | 16.65 | -6.7e-05 |
| [0/90/0], Nx = 1 N/m: σ1 in 0° ply [Pa] | 97.26 | 97.26 | +4.0e-05 |
| [0/90/0], Nx = 1 N/m: σ2 in 90° ply [Pa] | 5.472 | 5.472 | +2.6e-05 |
| 60° lamina, max stress S [MPa] | 16.33 | 16.33 | -2.3e-04 |
| 60° lamina, max strain S [MPa] | 16.33 | 16.33 | -2.3e-04 |
| 60° lamina, Tsai-Hill S [MPa] | 10.94 | 10.94 | -2.2e-04 |
| 60° lamina, modified Tsai-Hill S [MPa] | 16.06 | 16.06 | +4.4e-05 |
| 60° lamina, Tsai-Wu S [MPa] | 22.39 | 22.39 | -5.7e-05 |
| [0/90/0] first ply failure Nx [N/m] (90°, Tsai-Wu) | 7.277e+06 | 7.277e+06 | -6.0e-05 |
| [0/90/0] last ply failure Nx [N/m] (0°, discounted) | 1.5e+07 | 1.5e+07 | +0.0e+00 |

T300/5208 graphite/epoxy (Kaw Table 2.1), 5 mm plies as in Kaw's examples. Kaw prints four significant figures, so ~1e-3 agreement is the best possible.
