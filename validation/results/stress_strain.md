# Stress/strain validation

| Check | Computed | Reference | Source |
|---|---|---|---|
| Principal stresses, (−20, 90, 60) MPa | 116.39 / -46.39 MPa | 116.39 / -46.39 MPa | closed form |
| Shear on principal plane | 7.5e-09 Pa | 0 | exact |
| Pure shear: principal values, angle | ±40.0 MPa at 45.0° | ±40.0 MPa at 45° | exact |
| Invariants under random rotation (max rel. change) | 4.5e-16 | 0 | round-off |
| von Mises, uniaxial / pure shear | 1.0000 σ / 1.7321 τ | 1 σ / √3 τ = 1.7321 τ | exact |
| Tresca, pure shear | 2.0000 τ | 2 τ | exact |
