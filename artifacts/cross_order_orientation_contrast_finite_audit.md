# Exact finite audit of orientation selection versus rare bridges

Status: **reproducible finite computation, not an asymptotic theorem**.  The
audit uses exact pressure-minimizing children at each child's own scale,
enumerates both orientations and the complete bridge cube, and evaluates
quantities which directly upper-bound the cross-order defect.

For `N=m+n`, let `A,D` minimize `P_m(beta),P_n(beta)` and put

```math
L_\epsilon(B)=\log\overline Z_N
 \left(\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix},
 {\beta\over\sqrt N}\right),
\qquad T=P_m(\beta)+P_n(\beta).
```

The three reported certificates are

```math
U={1\over2}\sum_\epsilon\mathbb E_BL_\epsilon(B)-T,
\quad
O=\mathbb E_B\min_\epsilon L_\epsilon(B)-T,
\quad
R=\min_{\epsilon,B}L_\epsilon(B)-T.                   \tag{1}
```

They satisfy

```math
\boxed{E_{m,n}(\beta)\le R\le O\le U,}
\qquad
U-O={1\over2}\mathbb E_B|L_+(B)-L_-(B)|.             \tag{2}
```

Thus every number below is already in the cross-order normalization; no
surrogate cap, conference child, or target-order minimizer is used.

## Balanced and nearest-balanced data

The columns are `(U, U-O, O, R)`.

| `N` | split | `beta` | uniform | orientation gain | selected mean | best rare bridge |
|---:|:---:|---:|---:|---:|---:|---:|
| 4 | 2+2 | 0.5 | 0.063080 | 0.003585 | 0.059495 | 0.059495 |
| 6 | 3+3 | 0.5 | 0.063690 | 0.005563 | 0.058127 | 0.041401 |
| 8 | 4+4 | 0.5 | 0.069346 | 0.005883 | 0.063462 | 0.045942 |
| 9 | 4+5 | 0.5 | 0.074550 | 0.005976 | 0.068574 | 0.043996 |
| 10 | 5+5 | 0.5 | 0.080711 | 0.006208 | 0.074503 | 0.040956 |
| 4 | 2+2 | 1 | 0.254574 | 0.043727 | 0.210847 | 0.210847 |
| 6 | 3+3 | 1 | 0.258160 | 0.065087 | 0.193073 | 0.015572 |
| 8 | 4+4 | 1 | 0.328465 | 0.068751 | 0.259714 | 0.074212 |
| 9 | 4+5 | 1 | 0.388315 | 0.068209 | 0.320106 | 0.075026 |
| 10 | 5+5 | 1 | 0.459155 | 0.069645 | 0.389510 | 0.078234 |
| 4 | 2+2 | 2 | 0.912557 | 0.276975 | 0.635582 | 0.635582 |
| 6 | 3+3 | 2 | 0.797231 | 0.364482 | 0.432749 | -0.645080 |
| 8 | 4+4 | 2 | 1.318982 | 0.385598 | 0.933384 | -0.159123 |
| 9 | 4+5 | 2 | 1.761474 | 0.372630 | 1.388844 | 0.055127 |
| 10 | 5+5 | 2 | 2.264344 | 0.386573 | 1.877771 | 0.451428 |
| 10 | 5+5 | 4 | 7.396641 | 1.076901 | 6.319740 | 1.362881 |

At `beta=1`, the joint uniform defect grows from `0.258` at `N=6` to
`0.459` at `N=10`, while the free-orientation gain stays near `0.07` and
the globally optimized bridge certificate stays below `0.079`.  The
moderate-temperature data therefore separate two mechanisms cleanly:

1. orientation selection alone does not track the growing averaged defect;
2. a very small subset of correlated bridges can remove almost all of it.

This is consistent with an `O_beta(1)` or another sublinear rare-bridge
certificate, but five small orders cannot distinguish that from delayed
linear growth.  It is evidence for attacking `R` directly, not evidence for
the desired recurrence.

The order-ten `beta=1` bridge basin also separates ordinary fixed fractional
moments from the true minimum.  Among the `2^(25+1)` joint
orientation/bridge choices, the exact masses below `T+s` are

| `s` | 0.05 | 0.10 | 0.25 | 0.50 | 1.00 |
|---:|---:|---:|---:|---:|---:|
| mass | 0 | `5.9605e-7` | 0.015091 | 0.692006 | 0.995344 |

The corresponding fractional certificates
`R_lambda-T=-(1/lambda)log E exp(-lambda L)-T` are

| `lambda` | 0.25 | 0.5 | 1 | 2 | 4 | 8 |
|---:|---:|---:|---:|---:|---:|---:|
| defect | 0.456838 | 0.454592 | 0.450291 | 0.442354 | 0.428534 | 0.406359 |

Thus the exact best defect `0.078234` is carried by an exponentially small
finite basin and is invisible to these fixed fractional parameters.  This
is precisely the quantitative basin obligation in the selector-transport
audit; it is not a proof that the exponential rate persists.

At `beta=0.25`, the uniform defects remain between `0.0156` and `0.0169`
through order ten, so this range is too small to detect the conditional
`Theta(beta^4 N)` high-temperature obstruction.  At `beta=2`, the best
certificate changes sign, is nonmonotone, and rises to `0.451` at order ten;
at `beta=4` it rises to `1.363`.  The apparent boundedness at `beta=1` is
therefore not a temperature-uniform phenomenon.  None of these observations
supports an asymptotic claim.

## Verification and reproduction

The signing cube, bridge cube, and all integer energies are enumerated
exactly.  Absolute energy histograms are compared at 80 decimal digits to
select every child thermal-minimizer class.  The pressure and convolution
evaluations are floating point, with three independent direct bridge checks
per cube and observed log-pressure error below `2e-9`.  The identity in
(2) is checked to absolute error below `2e-10` in every record.

```bash
.venv/bin/python computations/cross_order_orientation_contrast_exact.py \
  --max-total-n 9 \
  --output computations/results/cross_order_orientation_contrast_exact.json

.venv/bin/python computations/cross_order_orientation_contrast_exact.py \
  --min-total-n 10 --max-total-n 10 --balanced-only --betas 1 \
  --output computations/results/cross_order_orientation_contrast_exact_n10_beta1.json

.venv/bin/python computations/cross_order_orientation_contrast_exact.py \
  --min-total-n 10 --max-total-n 10 --balanced-only --betas 0.25 0.5 2 4 \
  --output computations/results/cross_order_orientation_contrast_exact_n10_other_betas.json
```

The program deliberately labels `R` as a direct construction certificate,
not a value of the globally minimized order-`N` pressure.
