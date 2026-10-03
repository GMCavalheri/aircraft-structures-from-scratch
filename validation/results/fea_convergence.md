# FEA mesh convergence

| Study | Result | Expected |  |
|---|---|---|---|
| Consistent loads, nodal midspan deflection (n = 2…64) | max error 2.8e-11 | exact nodal values | round-off |
| Lumped loads, nodal midspan deflection | order 2.00 | 2 | ✔ |
| Hermite interpolation between nodes | order 4.00 | 4 | ✔ |
