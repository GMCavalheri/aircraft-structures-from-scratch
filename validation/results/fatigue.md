# Fatigue validation

| Check | Computed | Reference | Source |
|---|---|---|---|
| ASTM E1049 example, (range, count) | 3: 0.5; 4: 1.5; 6: 0.5; 8: 1; 9: 0.5 | 3: 0.5; 4: 1.5; 6: 0.5; 8: 1; 9: 0.5 | ✔ exact |
| 2024-T3 Kt = 1, Smax = 40 ksi, R = 0: log Nf | 5.6063 | 5.6063 | hand calc. of Fig. 3.2.3.1.8(e) |
| Refit of 120 synthetic tests (handbook scatter): A1, A2, A3, A4 | 8.68, -3.04, 0.680, 13.5 | 9.20, -3.33, 0.680, 12.3 | Fig. 3.2.3.1.8(g) |
| Refit standard error of log life | 0.252 | 0.270 | Fig. 3.2.3.1.8(g) |

A1, A2 and A4 are strongly correlated (a steeper slope with a higher fatigue limit describes nearly the same curve), so the individual coefficients wander while the refitted curve stays within a factor of 1.5 in life of the handbook curve over 10⁴–10⁶ cycles (tested).
