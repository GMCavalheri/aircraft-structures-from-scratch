# Buckling validation: FEA eigenvalue vs Euler

| End condition | K | Euler P_cr | error, 4 el. | 8 el. | 16 el. | mode-shape error, 8 el. |
|---|---|---|---|---|---|---|
| pinned-pinned | 1.0000 | 4797.72 N | +5.1e-04 | +3.3e-05 | +2.1e-06 | 7.8e-16 |
| fixed-free | 2.0000 | 1199.43 N | +3.3e-05 | +2.1e-06 | +1.3e-07 | 1.3e-14 |
| fixed-fixed | 0.5000 | 19190.90 N | +7.5e-03 | +5.1e-04 | +3.3e-05 | 3.1e-15 |
| fixed-pinned | 0.6992 | 9814.94 N | +2.1e-03 | +1.4e-04 | +8.6e-06 | 3.2e-07 |

E = 70 GPa, I = 1e-8 m⁴, L = 1.2 m. Consistent geometric stiffness: the FEA load is an upper bound and the error falls as h⁴. On a uniform mesh the nodal deflections of the three trigonometric modes are exact to round-off; the fixed-pinned mode converges at about h⁶.
