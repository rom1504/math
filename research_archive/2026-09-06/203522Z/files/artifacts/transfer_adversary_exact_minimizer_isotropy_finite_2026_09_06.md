# Exact minimizers need not have isotropic exact-ground laws: a finite counterexample

Date: 2026-09-06. This is a finite theorem and a bounded computational
screen, not an asymptotic exact-minimizer obstruction. The separate
scalable near-minimizer counterexample is in
`transfer_adversary_nearmin_clique_orientation_gap_2026_09_06.md`.

Write `H_A(x)=sum_(i<j) A_ij x_i x_j`, `Q(A)=max_x |H_A(x)|` and
`g_A(x)=Q(A)-|H_A(x)|`. An isotropic spin law means `E xx^T=I`.
Its mean is immaterial here; symmetrizing by an independent global sign
can also make the mean zero without changing any statement.

## 1. Elementary order-four counterexample

Consider the hollow signing

```text
A = [ 0  1  1  1 ]
    [ 1  0  1  1 ]
    [ 1  1  0 -1 ]
    [ 1  1 -1  0 ].
```

For every spin `x in {−1,1}^4`, direct conditioning on the two pair
products gives the exact identity

```math
|H_A(x)|=2+x_0x_1+x_2x_3.                             (1)
```

Indeed, if both pair products are `+1`, the two internal-edge
contributions cancel and the four cross edges have absolute sum four.
If both are `−1`, the internal and cross sums vanish. If just one is
`−1`, the cross sum vanishes and the internal sum has absolute value two.
Thus `Q(A)=4`.

For ANY order-four signing, the uniform-cube second moment is
`E H_A(x)^2=6`, by orthogonality of distinct pair characters. All its
energies are even integers. Hence `Q(A)<4` would imply `Q(A)<=2` and
`E H_A^2<=4`, a contradiction. Therefore `M_4=4`, and the displayed
matrix is an exact minimizer with a self-contained lower certificate.

Its projective absolute grounds, using `x_0=1`, are exactly

```text
(1,1,1,1),  (1,1,-1,-1).
```

Both have `x_0x_1=1`. No probability law supported on their global-sign
orbits can therefore be isotropic. This is not a failure of a numerical
solver or merely a numerical rank test.

More sharply, (1) implies for EVERY isotropic law

```math
E g_A(x)=E[2-x_0x_1-x_2x_3]=2.                       (2)
```

The uniform law on the following four mutually orthogonal spins has
`E xx^T=I`; each spin has absolute energy two:

```text
(1, 1, 1,-1),  (1, 1,-1, 1),
(1,-1, 1, 1),  (1,-1,-1,-1).
```

Consequently the smallest possible maximum support slack is exactly
two, and the smallest possible mean slack is exactly two. In particular
no exact-ground isotropic law exists, even though `A` is globally
minimizing among the discrete signings.

The accompanying checker additionally enumerates all 64 labeled
order-four signings. Exactly 48 have cap four, and every one has exactly
two projective absolute ground states. Each such support spans at most
two dimensions, so NONE of the order-four exact minimizers has an
isotropic exact-ground law. This stronger finite classification is
computational; the displayed single counterexample above is elementary.

## 2. What discrete global minimality does force

There is a useful finite one-edge condition, weaker than a common law.
If `A` is an exact minimizer and `e={i,j}`, let `A^e` reverse just that
edge. Since `Q(A^e)>=Q(A)`, some oriented spin `(sigma,x)`, with
`sigma in {−1,1}`, satisfies

```math
\sigma H_A(x)-2\sigma A_{ij}x_ix_j\ge Q(A).
```

The slack `Q(A)-sigma H_A(x)` is nonnegative. The last display therefore
forces

```math
\sigma A_{ij}x_ix_j=-1,
\qquad 0\le Q(A)-\sigma H_A(x)\le2.                  (3)
```

The witness in (3) may depend on the edge. No common distribution, no
zero pair correlations, and no continuous stationarity conclusion has
been established by this argument. The order-four counterexample shows
that exact-ground isotropy cannot simply be substituted for (3).

## 3. Preserved bounded LP screen and exact replay

Source: `computations/transfer_adversary_minimizer_isotropy_2026_09_06.py`.
Data and full rational certificates:
`computations/results/transfer_adversary_minimizer_isotropy_2026_09_06.json`.
Run from the repository root with

```text
.venv/bin/python computations/transfer_adversary_minimizer_isotropy_2026_09_06.py --output computations/results/transfer_adversary_minimizer_isotropy_2026_09_06.json
```

The final run passed all 15 cases. This is a diagnostic sample above
order eight, not an enumeration of all minimizer classes at those orders.
The table records exact rational results; `d_max` is minimum allowable
maximum support slack and `d_mean` is minimum mean slack.

| Stored case | Cap | Positive / negative projective grounds | d_max | d_mean |
|---|---:|---:|---:|---:|
| n3_class0 | 3 | 1 / 0 | 2 | 3/2 |
| n4_class0 | 4 | 1 / 1 | 2 | 2 |
| n5_class0 | 4 | 5 / 5 | 4 | 2/3 |
| n6_class0 | 5 | 6 / 6 | 0 | 0 |
| n7_class0 | 9 | 0 / 4 | 6 | 12/5 |
| n7_class1 | 9 | 7 / 0 | 6 | 19/6 |
| n7_class2 | 9 | 0 / 3 | 4 | 2 |
| n8_class0 | 10 | 4 / 4 | 2 | 6/5 |
| n8_class1 | 10 | 4 / 4 | 2 | 2/3 |
| n9_stored | 12 | 9 / 12 | 4 | 4/3 |
| n10_stored | 13 | 20 / 20 | 2 | 1 |
| n11_stored | 17 | 0 / 5 | 4 | 2 |
| n12_stored | 18 | 10 / 10 | 0 | 0 |
| n13_stored | 20 | 39 / 39 | 4 | 10/7 |
| n14_stored | 21 | 78 / 78 | 0 | 0 |

All cap values are replayed by exact integer enumeration of the
`2^(n-1)` spins with first coordinate one. Feasible laws have nonnegative
rational weights, sum exactly one, and every pair moment exactly zero.
Each infeasible support threshold has a rational quadratic separator
whose value is strictly positive on EVERY allowed spin and whose
expectation under isotropy is zero. Minimum mean slack has both an
exact rational primal law and a pointwise rational dual inequality,
with exactly equal objectives. Floating-point LPs only discover these
certificates; they do not certify the conclusions. Numerical support
rank, also retained in the JSON, is not used as a proof.

### Data provenance and imported optimality status

For orders three through eight all stored orbit representatives come
from `computations/results/m{n}_minimizer_orbits.json`, with JSON pointers
`/classes/{class}/representative_matrix`. Above that range the inputs are:

| n | Input in computations/results/ | Matrix key |
|---:|---|---|
| 9 | exact_m9.json | matrix |
| 10 | exact_m10.json | matrix |
| 11 | nested_10_in_11_cap17.json | matrix |
| 12 | extension_nested_m11_to_12.json | parent_matrix |
| 13 | bridge_6_7_sign1_cap20.json | parent_matrix |
| 14 | conference_completion_m13.json | conference_matrix |

Every case in the output stores its input path, SHA-256, JSON pointer,
imported minimum cap, and imported certificate status. The lower
certificates for orders eleven/twelve and thirteen/fourteen are in
`certified_m11_m12.json` and `certified_m13_m14.json`, respectively.
Their solver-dependent global lower bounds were NOT rerun here and are
not converted by this screen into standalone proof objects. The
order-four theorem in Section 1 has no such dependency. The screen
depends on NumPy, SciPy/HiGHS and SymPy for discovery and rational
reconstruction; the identities checked afterward are exact.

### Failed preliminary runs, retained as research provenance

The first exploratory run stopped at the order-ten slack-two feasible
law because independent bounded-denominator rounding did not reconstruct
the exact pair constraints. The correction added exact SymPy affine
row reduction on the discovered support; no tolerance was loosened to
accept an approximate law. A second run stopped at order eleven because
the input key was mistakenly `parent_matrix`; the stored file actually
uses `matrix`. That input mapping was corrected. Neither stopped run
wrote a final result file, because output is written only after every
case passes. The final complete rerun checked all 15 cases and the
elementary order-four identities.

## 4. Boundary of the conclusion

The finite table is compatible both with success and with failure of an
eventual exact-minimizer small-slack isotropy principle. In particular,
fixed positive finite slacks are not an asymptotic obstruction to
`o(sqrt(n))` slack. No pattern in these orders is promoted to a theorem
about large-order exact minimizers, a nested sequence, or convergence
of `M_n/n^(3/2)`. The separate planted-clique theorem rules out the
corresponding assertion for ALL asymptotically minimizing sequences,
not for all exactly minimizing sequences.
