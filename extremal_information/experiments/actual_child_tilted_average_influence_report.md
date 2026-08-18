# Tilted average influence on the actual-child IC path

Classification: **complete finite bridge enumeration and numerical
quadrature; reproducible evidence, not an interval certificate and not an
asymptotic theorem**.

> **Precision correction (2026-08-18).**  This report inherits the archived
> child list from `actual_child_projective_synchronization.py`.  That run used
> default `mpmath` precision with an inappropriately loose minimizer tie
> tolerance: only `36` of the `102` archived records survive corrected
> minimizing-child classification.  Per-record numerical identities remain
> implementation checks, but the quoted all-record ranges and optimizer-level
> empirical conclusions are withdrawn pending a deliberately regenerated run.
> The frozen JSON is retained as provenance rather than rewritten.

The experiment extends
[`actual_child_projective_synchronization.py`](actual_child_projective_synchronization.py)
to the exact hybrid law

```math
 {dq_s\over dr}\propto e^{-s h},\qquad 0\le s\le1,
```

for every balanced contracted-temperature minimizer class at
`N=4,...,9`, both orientations, and `beta=1,2,4`.  All bridges are
enumerated.  Twelve-node Gauss--Legendre quadrature evaluates

```math
 \int_0^1(1-s)\operatorname{Var}_{q_s}(h)\,ds
 \quad\hbox{and}\quad
 \int_0^1\mathcal A_s\,ds,
```

where `A_s` is exactly IC.13, using
`bar h_i=E_(r_i)h(R_i,B_(-i))` rather than a sampled approximation.
It also evaluates the optimized conditional-entropic influence `E_s` from
IC.4 and the exact forward path split

```math
D(q_s\Vert r)
=\operatorname {TC}_{\rm row}(q_s)
 +\sum_iD(q_{s,i}\Vert r_i),                         \tag{TI.1}
```

whose two terms integrate separately to the canonical reverse error:

```math
\mathcal J=D(r\Vert q_1)
=\int_0^1 {\operatorname {TC}_{\rm row}(q_s)\over s^2}\,ds
 +\int_0^1 {\sum_iD(q_{s,i}\Vert r_i)\over s^2}\,ds.\tag{TI.2}
```

## Numerical identity and inequality checks

Across all `102` actual-child laws, the first quadrature divided by the
directly evaluated canonical `J` differs from one by at most `4.37e-11`.
This independently verifies the IC.7 curvature identity and every path
normalization.

The independently evaluated two terms in (TI.2) sum to `J` with maximum
relative error `5.29e-11`.  This checks the KL directions and row-marginal
normalizations in the total-correlation decomposition.  Notice that the
integrand uses the forward divergence `D(q_s||r)`, even though its weighted
path integral recovers the reverse endpoint divergence `D(r||q_1)`.

The IC.15 right side divided by `J` lies in

| `beta` | range of quadrature `int A_s ds / J` |
|---:|---:|
| `1` | `[1.9664,2.0708]` |
| `2` | `[1.6973,2.2374]` |
| `4` | `[1.7348,2.6334]` |

Thus the tilted influence theorem is a constant-factor diagnostic on every
enumerated actual law, unlike the struck worst-context projective bound.
The factor near two at `beta=1` is consistent with the quadratic/small-
interaction regime; it is observed, not asserted asymptotically.

The sharper IC.23 conditional-entropic right side divided by `J` lies in

| `beta` | range of quadrature `int E_s ds / J` |
|---:|---:|
| `1` | `[1.9599,2.0672]` |
| `2` | `[1.6585,2.1655]` |
| `4` | `[1.4864,2.2731]` |

Thus optimizing each one-row conditional comparison improves the coarse
influence certificate, but it does not close the finite interaction.
By ES.23, the excess over one in this ratio is exactly the integrated dual
total correlation divided by `J`.  It ranges over `[.9599,1.0672]`,
`[.6585,1.1655]`, and `[.4864,1.2731]` for `beta=1,2,4`, respectively.
Thus every enumerated actual-child path has a dual-total-correlation mass
comparable to its canonical error even though every conditional row has
uniformly bounded Renyi-two complexity.  At `beta=2,4`, where the observed
`J/N` does not trend toward zero, this is direct finite evidence for the
"tight components but collective dependence" side of the structural
dichotomy; it remains numerical rather than an asymptotic lower bound.

## Which row component carries the canonical mismatch?

The exact integrated shares in (TI.2), across all orders `4,...,9`, are:

| `beta` | row total-correlation share | row-marginal-drift share |
|---:|---:|---:|
| `1` | `[.9072,1.0000]` | `[0,.0928]` |
| `2` | `[.5403,1.0000]` | `[0,.4597]` |
| `4` | `[.3462,1.0000]` | `[0,.6538]` |

At the largest enumerated order `N=9`, the split is:

| `beta` | row total-correlation share | row-marginal-drift share |
|---:|---:|---:|
| `1` | `[.9072,.9972]` | `[.0028,.0928]` |
| `2` | `[.5883,.8771]` | `[.1229,.4117]` |
| `4` | `[.3462,.9234]` | `[.0766,.6538]` |

This gives a sharp finite falsifier to either one-component explanation.
At `beta=1`, essentially all canonical error is genuine dependence among
rows while the one-row marginals remain near their canonical factors.  At
`beta=4`, the most strongly interacting `N=9` orientation has
`J/N=.43119` but only `.34618` of it in total correlation; `.65382` is
ordinary marginal drift.  Therefore a uniform theorem that attributes the
actual-child canonical mismatch solely to irreducible row dependence is
false already at the enumerated orders.  Conversely, a theorem controlling
only row marginals misses over ninety percent of some actual-child errors.

## Scaling data

The main quantities at the largest enumerated order are:

| `N=9` | range of `J/N` | range of `sup_s A_s/N` | IC.15/J range |
|---|---:|---:|---:|
| `beta=1` | `[.0008662,.0021510]` | `[.0018923,.0044019]` | `[1.9664,2.0663]` |
| `beta=2` | `[.0253895,.0458262]` | `[.0506598,.0860465]` | `[1.7955,2.1140]` |
| `beta=4` | `[.1263989,.4311872]` | `[.4012773,1.1670505]` | `[1.7348,2.6334]` |

For comparison across orders:

| `beta` | `N` | range of `J/N` | range of `sup_s A_s/N` |
|---:|---:|---:|---:|
| 1 | 4 | `.000696` | `.001434` |
| 1 | 6 | `[.001024,.001132]` | `[.002193,.002363]` |
| 1 | 8 | `[.000996,.002069]` | `[.002157,.004361]` |
| 1 | 9 | `[.000866,.002151]` | `[.001892,.004402]` |
| 2 | 4 | `.023432` | `.057537` |
| 2 | 6 | `[.029373,.035317]` | `[.063823,.076616]` |
| 2 | 8 | `[.030862,.054262]` | `[.070456,.103959]` |
| 2 | 9 | `[.025389,.045826]` | `[.050660,.086047]` |
| 4 | 4 | `.135677` | `.501620` |
| 4 | 6 | `[.209098,.239638]` | `[.511874,.698902]` |
| 4 | 8 | `[.272582,.364710]` | `[.683318,.846890]` |
| 4 | 9 | `[.126399,.431187]` | `[.401277,1.167051]` |

The maximizing grid point is an endpoint (`s=0` or `s=1`) in every case.
At `beta=1` it is always `s=0`; at the larger temperatures the endpoint
depends on the child class/orientation.

## Judgment

This audit touches the actual optimized-child law and survives the finite
falsification test: the IC.3 observable tracks the canonical interaction
within a factor below `2.64` on every enumerated instance.  It does **not**
show that either quantity is sublinear.  In fact the `beta=2,4` data retain
nonzero normalized influence throughout the observed orders, so finite
evidence does not support a no-gain closure there.

The narrower next question must retain both mechanisms: prove that the
weighted path mass of **both** row total correlation and marginal drift is
`o(N)` for contracted-temperature minimizing children, or prove that one
has a fixed positive density along a subsequence.  Conditional Renyi
tightness alone does not distinguish these alternatives.  The complete
records, including every `s` profile and both exact integrated components,
are in
[`../../computations/results/actual_child_projective_synchronization.json`](../../computations/results/actual_child_projective_synchronization.json).
