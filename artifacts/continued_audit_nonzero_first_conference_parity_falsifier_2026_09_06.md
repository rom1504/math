# An actual nonzero-first conference parity and covariance falsifier

Date: 2026-09-06. This is an actual-signing obstruction to transferring
the zero-first sine residual covariance unchanged to a nonzero-first
source. It does not demonstrate cancellation of the entire innovation;
in this example the residual instead has positive even-noise mass.

Let B be any normalized symmetric conference signing, so B is hollow,
its off-diagonal entries are signs divided by sqrt(n-1), and B^2=I.
There are arbitrarily large actual examples: the Paley conference
construction over F_q for q=5^k uses index set {infinity} union F_q,
unit infinity row/column, and quadratic characters chi(a-b) on distinct
finite indices. Symmetry follows from chi(-1)=1, and the elementary
quadratic-character sum gives the exact unnormalized identity C^2=qI.
For completeness, for distinct a,b the finite-index row inner product
is `1+sum_c chi((c-a)(c-b))=0`. Translate c by (a+b)/2; the remaining
character sum is `sum_x chi(x^2-d)=-1`, d nonzero. Counting the solutions
of `x^2-y^2=d` by the nonzero factor `x-y` gives q-1, proving that sum.
The diagonal row inner products are q, and the infinity/finite inner
products vanish because the nontrivial character sums to zero.

Fix `f(g)=kappa sin(g)` with `0<kappa<=1/4`, and a constant mask
`0<h<=1/4`. Its Gaussian first coefficient is
`b=kappa exp(-1/2)>0`. Let

```math
r(g)=f(g)-bg,\qquad G=BS,\quad Z=Br(G),\quad
t=\frac{\pi}{2b}.
```

The source f is bounded, but r has a linear part; it belongs to the
proved nonzero-first coherent-return scope, not the bounded zero-first
source hypothesis by itself. The literal coherent return is

```math
Bf(G)=bS+Z.
```

Consequently the actual first response satisfies the pointwise identity

```math
C=h\sin(t(bS+Z))=hS\cos(tZ).                         (1)
```

This is not a Gaussian substitution for QS: QS=S exactly. In particular
EVERY odd conditional-noise Hermite coefficient of C is identically
zero. The cubic residual from the pure zero-first sine response is absent;
the first nonconstant conditional-noise term is quadratic times the
own spin S. Thus its disappearance is caused by the retained coherent
shift, not by an abstract covariance model or an operator-limit guess.

The ideal old-noise covariance is especially explicit. Since Q=I,

```math
\tau^2=\kappa^2e^{-1}(\sinh(1)-1)>0,\qquad T=\tau^2I,
\quad v=t^2\tau^2=\frac{\pi^2}{4}(\sinh(1)-1).
```

The exact coherent conditional mean and mean derivative used by the
first-history theorem are

```math
c_0=h e^{-v/2}S,\qquad a=0.
```

Let `L_0=Bc_0` and `eta=B(C-c_0)`. Orthogonality of B gives exact finite
trace identities. The already proved LOCAL Gaussian comparison of Z
at each root then yields

```math
\frac1n E\|\eta\|^2
\longrightarrow
h^2\operatorname{Var}(\cos(\sqrt v N))
=\frac{h^2}{2}(1-e^{-v})^2>0,                         (2)
```

and `n^-1 E eta^T L_0 ->0`. No full new cross-root covariance theorem is
needed for either conclusion, since B^2=I reduces their traces to
one-root expectations. In fact S disappears from both trace expressions
because S_i^2=1, so only the marginal noise comparison is necessary.
It is part of the bounded nonzero-first retained-QS theorem; no
Gaussian replacement of the actual spin S is made.

If one incorrectly reused the pure zero-first sine residual matrix,
its normalized returned variance here would instead be

```math
h^2e^{-v}(\sinh(v)-v).
```

The difference between (2) and this value is

```math
h^2e^{-v}\,[v-1+e^{-v}]>0,                            (3)
```

because `v>0` and `e^{-v}>1-v`. Thus this is a genuine nonvanishing
trace discrepancy on actual signings, which also rules out normalized
nuclear convergence to that incorrectly transferred residual matrix.

Equation (3) falsifies the covariance FORMULA, not every inequality that
one might derive from a cubic lower matrix; an unrelated lower inequality
can hold accidentally even when its proposed Hermite mechanism is absent.
The example does not show `E eta_i^2+E eta_i(L_0)_i=0` on all roots.
Indeed its positive averaged variance and vanishing signed averaged
cross show that complete cancellation is impossible here. What remains
open is a general positive weighted innovation theorem for nonzero-first
or marked sources with their full literal coherent history retained.

Reproducible finite checks are in
`computations/continued_audit_conference_parity_falsifier_2026_09_06.py`
and its same-stem JSON under `computations/results/`. The conference
matrix identities are checked in integer arithmetic. All Boolean seeds
are enumerated at orders 6,14,18; order 30 uses 16,000 fixed-seed samples.
The trigonometric identity (1) is checked to floating error below 10^-12.
At kappa=h=1/4, the proved limiting variance is approximately .00384959,
whereas the incorrectly transferred prediction is .000551283, a
positive limiting discrepancy of approximately .00329831. The order-18
exhaustive average is .00344416. These finite transcendental evaluations
are diagnostics, not substitutes for the analytic positive difference
in (3), and the order-30 sampling uncertainty is recorded explicitly.
