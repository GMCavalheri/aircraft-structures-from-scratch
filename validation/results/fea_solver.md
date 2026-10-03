# FEA solver validation

| Case | FEA | Reference | Source | Rel. error |
|---|---|---|---|---|
| Cantilever tip deflection (4 el.) | 0.321429 | 0.321429 | PL³/3EI | -2.3e-14 |
| Simply supported UDL, midspan (6 el.) | 0.015067 | 0.015067 | 5wL⁴/384EI | -3.2e-15 |
| Propped cantilever, prop reaction (8 el.) | 2250 | 2250 | 3wL/8 | +0.0e+00 |
| Overhang beam, tip deflection vs Macaulay solver | -0.00342229 | -0.00342229 | Phase 2 | -3.5e-15 |
| Overhang beam, moment over roller vs Macaulay | -912 | -912 | Phase 2 | -1.7e-15 |
| L-frame tip deflection (2 el.) | 0.123786 | 0.123786 | Pb³/3EI + Pb²h/EI + Ph/EA | +1.9e-14 |
| Pratt truss: support reaction | 25.000 kN | 25.000 kN | statics | -8.7e-16 |
| Pratt truss: midspan bottom-chord force | 40.000 kN | 40.000 kN | method of sections | -2.7e-15 |
| Portal frame: ΣRx, ΣRy | -40.000, 80.000 kN | −40.000, 80.000 kN | equilibrium | 2.1e-14 |
