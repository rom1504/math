# Stable audit and stress test of the greedy Gibbs bridge

Status: exact finite-state normalization and increment audit; reproducible
finite scaling evidence. No uniform gradient theorem or scalable obstruction
is proved.

## Independent normalization audit

For child matrices of orders `m,n`, the implementation fixes one global spin
in the first child and retains both orientations of the second child. Hence it
enumerates exactly

```math
2^{m-1}2^n=2^{m+n-1}
```

projective parent spin states. With current energy `E`, its stable
log-sum-exp calculation evaluates

```math
Z=\mathbb E\cosh(\gamma E),\qquad
r_{ij}=\frac{\mathbb E[x_i y_j\sinh(\gamma E)]}{Z}.
```

Assigning an unrevealed edge the sign `c` therefore has the exact ratio

```math
\frac{Z_{\rm new}}{Z}
=\cosh\gamma\left(1+c r_{ij}\tanh\gamma\right).
```

The independent program replayed the saved greedy trajectories for equal
exact children at every order `4` through `10`. It checked that each saved
edge maximized the absolute response among unrevealed edges, that every
nonzero response received the opposite sign, that every one-edge increment
agreed with the direct partition-function difference, and that the final
integer cap was unchanged. Maximum errors were `7.8e-16` in a response,
`8.8e-15` in an increment, and `3.4e-16` in greedy maximality. At exact-zero
responses the two floating summation orders can choose opposite signs; both
are mathematically minimizing choices.

## Calibrated threshold behavior

Let `L=mn` and

```math
\Delta=(M_m^{2/3}+M_n^{2/3})^{3/2}-M_m-M_n.
```

The exact tested greedy cost is

```math
K=\sum_e\left[\log\cosh\gamma
 +\log(1-|r_e|\tanh\gamma)\right].
```

Thus the exact finite-temperature comparison is `gamma Delta-K`. The simpler
linear sufficient condition obtained from
`-log(1-z) >= z` has threshold

```math
\sum_e|r_e|\ge
\frac{L\log\cosh\gamma-\gamma\Delta}{\tanh\gamma}.
```

For scaled temperature `t=4`, `gamma=4/sqrt(2n)`, the saved equal-child
trajectories give:

| `n` | `sum abs(r) / n^(3/2)` | linear surplus `/ n^(3/2)` | exact `gamma Delta-K` | final cap |
|---:|---:|---:|---:|---:|
| 4 | 0.8471 | -0.2461 | 1.254 | 10 |
| 5 | 0.8508 | -0.4104 | -1.349 | 13 |
| 6 | 0.7949 | -0.4708 | -3.749 | 20 |
| 7 | 0.8896 | -0.1989 | -0.088 | 27 |
| 8 | 0.8538 | -0.2764 | -2.464 | 34 |
| 9 | 0.9093 | -0.2117 | -1.589 | 39 |
| 10 | 0.8745 | -0.2756 | -3.794 | 48 |

There is no improving normalized linear-threshold surplus: it remains
negative and of roughly constant scale. A log-log fit over these seven small
orders gives `sum abs(r) ~ n^1.57`; the normalized values themselves are the
more transparent evidence for `Theta(n^(3/2))` behavior at fixed `t`.
This is finite evidence, not an asymptotic law.

Failure of the linear sufficient condition must not be confused with failure
of the exact `K` comparison. At `n=4,t=4`, for example, the linear surplus is
negative while the exact margin is positive. Exact margins are also positive
at some tested orders at `t=2`. Nevertheless, at `t=4` the exact margin is
negative at six of the seven orders (slightly negative at `n=7`) and shows no
improvement with order. The temperatures required for a summable
soft-to-ground error must grow, so fixed-`t` success would not by itself prove
the desired recurrence.

## Random orders and held-out minimizers

For orders `4` through `8`, five deterministic random edge orders were tested
at each of `t=2,4`. The greedy order had a larger cumulative normalized
response than the random-order mean in every order-temperature pair. At
`t=4`, the greedy values ranged from `0.795` to `0.890`, while random-order
means ranged from `0.645` to `0.778`. Random final caps were generally equal
or worse. This supports the local greedy rule but supplies no objective-
independent reveal order or uniform lower bound.

All known exact-minimizer switching/permutation classes at orders `7` and `8`
were also held out and tested at `t=4`. The normalized cumulative responses
were `0.890, 0.973, 0.907` for the three order-seven classes and `0.867,
0.855` for the two order-eight classes. Their exact margins were respectively
`-0.088, 1.254, -0.066` and `-2.254, -2.468`. All members at a fixed order had
the same final cap (`27` or `34`), but the cost margins are not class-uniform.

## Stopping judgment

The data rule out neither the exact cumulative-gradient lemma nor the greedy
construction. They do show that the normalized first-moment shortfall is not
shrinking through order `10`, that random reveal orders do not repair it, and
that finite exact margins depend on the minimizer class. The first bridge
response is exactly zero by block spin-flip symmetry, but only one or a few
machine-zero steps occur afterward; there is no observed extensive zero-
gradient obstruction.

Computing the next choice still requires the full Gibbs response matrix of
the exponentially large parent state law. The experiment found no bounded
state, invariant, or scalable obstruction. Additional instances of the same
full-response heuristic would therefore not produce theorem-level leverage;
the route should resume only with a uniform cumulative-response potential or
a compressed state that certifies the exact `K` margin.

## Reproduction

```bash
PYTHONPATH=computations .venv/bin/python \
  computations/phase2l_gibbs_bridge_audit_stress.py \
  --output computations/results/phase2l_gibbs_bridge_audit_stress.json
```

Result SHA-256:
`302deb64aaefde1ac92d0f352a12549ab933b34bb258271e0ad7a93e2f3d0131`.
