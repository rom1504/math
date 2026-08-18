# Tilted average influence on the actual-child IC path

Classification: **complete finite bridge enumeration and numerical
quadrature; reproducible evidence, not an interval certificate and not an
asymptotic theorem**.

The experiment extends
[`actual_child_projective_synchronization.py`](actual_child_projective_synchronization.py)
to the exact hybrid law

```math
 {dq_s\over dr}\propto e^{-s h},\qquad 0\le s\le1,
```

for every balanced contracted-temperature minimizer class at
`N=4,...,9`, both orientations, and `beta=1,2,4`.  All bridges are
enumerated.  Twelve-node Gauss--Legendre quadrature evaluates both

```math
 \int_0^1(1-s)\operatorname{Var}_{q_s}(h)\,ds
 \quad\hbox{and}\quad
 \int_0^1\mathcal A_s\,ds,                           \tag{TI.1}
```

where `A_s` is exactly IC.13, using
`bar h_i=E_(r_i)h(R_i,B_(-i))` rather than a sampled approximation.

## Numerical identity and inequality checks

Across all `102` actual-child laws, the first quadrature divided by the
directly evaluated canonical `J` differs from one by at most `4.37e-11`.
This independently verifies the IC.7 curvature identity and every path
normalization.

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

The narrower next question is no longer a worst-context synchronization
bound.  It is an optimizer-specific theorem controlling the averaged
one-row collision--cavity increments under the tilted path, or a converse
showing that their observed extensive mass persists.  The complete records,
including every `s` profile, are in
[`../../computations/results/actual_child_projective_synchronization.json`](../../computations/results/actual_child_projective_synchronization.json).
