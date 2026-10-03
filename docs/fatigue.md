# Fatigue

Code: [`src/structures/fatigue/`](../src/structures/fatigue) ·
Notebook: [`06_fatigue.ipynb`](../notebooks/06_fatigue.ipynb)

## Stress-life curves

**Basquin**: $N = C\,S^{-m}$, a straight line in log–log axes. `fit_basquin` regresses log N on
log S, with life as the dependent variable as in ASTM E739.

**MIL-HDBK-5J equivalent stress.** The handbook folds every stress ratio $R = S_{min}/S_{max}$
onto one curve:

```math
\log_{10}N_f = A_1 + A_2\log_{10}(S_{eq} - A_4), \qquad S_{eq} = S_{max}(1-R)^{A_3}\quad\text{(ksi)}
```

$A_4$ behaves as a fatigue limit, so life is infinite when $S_{eq} \le A_4$. Since
$S_a = S_{max}(1-R)/2$, the equivalent stress is $S_{eq} = 2^{A_3}S_{max}^{1-A_3}S_a^{A_3}$, a Walker
mean-stress correction with exponent $A_3$.

| Curve | A1 | A2 | A3 | A4 (ksi) | std. error log N | Source |
|---|---|---|---|---|---|---|
| 2024-T3 sheet, unnotched | 11.1 | −3.97 | 0.56 | 15.8 | 0.38 | Fig. 3.2.3.1.8(e) |
| 2024-T3 sheet, Kt = 2.0 | 9.2 | −3.33 | 0.68 | 12.3 | 0.27 | Fig. 3.2.3.1.8(g) |
| 7075-T6 sheet, unnotched | 14.86 | −5.80 | 0.49 | — | 0.41 | Fig. 3.7.6.1.8(d) |
| 7075-T6 sheet, Kt = 2.0 | 7.50 | −2.46 | 0.54 | 18.6 | 0.31 | Fig. 3.7.6.1.8(f) |

![S-N curves](figures/sn_curves.png)

`fit_equivalent_stress` refits $A_1 \ldots A_4$ by nonlinear least squares on log life.
$A_1$, $A_2$ and $A_4$ are strongly correlated, so a refit of scattered data can return quite
different coefficients that describe nearly the same curve.

![Equivalent-stress fit](figures/sn_fit.png)

## Mean-stress corrections

Each correction gives the fully reversed amplitude $S_{ar}$ that is as damaging as amplitude
$S_a$ at mean $S_m$:

| Model | $S_{ar}$ |
|---|---|
| Goodman | $S_a/(1 - S_m/S_u)$ |
| Gerber | $S_a/(1 - (S_m/S_u)^2)$ |
| Soderberg | $S_a/(1 - S_m/S_y)$ |
| Smith–Watson–Topper | $\sqrt{S_{max}S_a}$ |
| Walker | $S_{max}^{1-\gamma}S_a^{\gamma}$ (γ = 0.5 is SWT) |

## Rainflow counting

ASTM E1049-85 Sec. 5.4.4:
1. Reduce the history to its reversals.
2. Push each reversal on a stack.
3. While the newest range X is at least the previous range Y, count Y and remove its points.
   Count it as a half cycle if Y contains the starting point, otherwise as a full cycle.
4. Count whatever is left on the stack at the end as half cycles.

The implementation reproduces the standard's example history (−2, 1, −3, 5, −1, 3, −4, 4, −2)
cycle by cycle, including the means and start/end indices.

## Cumulative damage

Palmgren–Miner: $D = \sum n_i/N_i$, failure at D = 1. Wholly compressive cycles
($S_{max} \le 0$) are taken as non-damaging.

## A synthetic flight spectrum

`flight_history` generates a seeded lower-wing-skin-like history:
- a ground–air–ground cycle every flight, from −4 ksi on the ground to the 12 ksi 1 g mean;
- 40 gust cycles per flight about the 1 g mean, with increments drawn from an exponential
  distribution (mean 2.5 ksi), i.e. an exponential exceedance curve.

It is not TWIST or FALSTAFF. The study applies it to the 2024-T3 Kt = 2 curve:

![Spectrum and rainflow matrix](figures/spectrum_rainflow.png)
![Predicted lives](figures/miner_life.png)

The handbook's R-dependent curve predicts about 2.1 × 10⁵ flights. A Basquin line fitted at
R = −1 gives 2.4 × 10⁵ flights with Goodman or Walker (γ = A3) corrections, but SWT halves the
life to 1.2 × 10⁵. The choice of mean-stress model matters as much as the S-N curve. This
spectrum is dominated by tensile means, where SWT is the most severe of the three.

## Validation

| Check | Result |
|---|---|
| ASTM E1049 rainflow example | exact, cycle by cycle |
| Range conservation on random histories | every reversal-to-reversal range counted once |
| Handbook equation | hand calculation, fatigue limit, inverse |
| Refit of 120 synthetic tests with the handbook's scatter | lives within a factor of 1.5; A3 within 0.05; std. error within 25 % |
| Mean-stress corrections | all reduce to $S_a$ at $S_m = 0$; Walker(γ = A3) = handbook $S_{eq}/2^{A_3}$ |
| Miner, constant amplitude | D = n/N |

## Limits

Stress-life only: no strain-life (low-cycle) or crack-growth (damage-tolerance) analysis.
Miner's rule ignores load-sequence effects such as overload retardation. Lives are medians
with no scatter factor. Real aircraft use reliability-based scatter factors and inspection
intervals from crack growth. The handbook equations should not be extrapolated beyond the
tested stress ratios or lives.
