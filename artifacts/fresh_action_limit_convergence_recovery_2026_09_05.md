# A precise action-limit recovery target for original convergence

Date: 2026-09-05. This is a structural reduction, not a proof of recovery
or of convergence. It was checked after the new bounded-operator and
weighted Gaussian modules were completed.

## 1. The primary theorem actually available

Primary source: A. Backhausz and B. Szegedy, *Action convergence of
operators and graphs*, Canadian Journal of Mathematics 74 (2022),
72--121; [primary preprint](https://arxiv.org/pdf/1811.00626).

In the inspected preprint, Theorem 2.16 gives compactness of action
equivalence classes with uniformly bounded `L^p -> L^q` norm for
`p<infinity`, `q>1`. Theorem 2.14 represents limits by operators, and
Proposition 3.3 preserves self-adjointness under these norm bounds.
These apply with `p=q=2`. The action profile records distributions of
bounded input functions together with their images under the operator.
Proposition 11.1 supplies subsequential action convergence for normalized
iid sign matrices; the introductory random-matrix discussion explicitly
leaves full-order convergence open. Thus no all-order sign-matrix
sampling theorem is imported here.

## 2. The original quadratic optimum is continuous in this topology

View `B=A/sqrt(n-1)` as an operator on the probability space consisting
of `n` equally weighted points. Define for any bounded self-adjoint
operator `T`

`E_quad(T)= (1/2) sup_(|f|<=1) |E[f(Tf)]|`.

For an actual hollow matrix,

`E_quad(B)=Q(A)/(n sqrt(n-1))`.

Indeed the quadratic form is affine in each coordinate separately, so
its maximum absolute value over the continuous cube equals that over
its vertices. This equality would need a diagonal correction for
arbitrary nonhollow matrices.

The functional `E_quad` is continuous on an action-compact class with
`||T||2->2<=L`. Its proof uses only the one-profile. Every law of
`(f,Tf)` has `|f|<=1` and `E|Tf|^2<=L^2`. Therefore the unbounded test
function `xy` is uniformly integrable: truncating its second coordinate
at magnitude `R` loses at most `L^2/R`. The truncated test is bounded
and Lipschitz, so its integral varies continuously under the
Levy--Prokhorov metric. Hausdorff convergence of the one-profiles then
passes both the upper and lower bounds of their suprema. Equivalently,
a coupling argument gives a modulus of order `C_L sqrt(delta)` after
optimizing the truncation radius against a profile error `delta`.

Consequently every bounded-operator sequence of original signings has
a subsequence on which its normalized optimum converges, with the limit
represented by `E_quad(T)` for an actual action-limit operator. This is
more informative than scalar subsequential compactness but still does
not compare different orders.

## 3. A sufficient recovery statement, independent of the scalar optimum

The following would suffice:

> For every action limit `T` of normalized hollow symmetric sign matrices
> with some fixed operator cap, there are normalized hollow symmetric
> sign matrices `C_N` at every sufficiently large order, with a fixed
> operator cap (allowed to depend on `T`), whose closed one-profiles
> converge to that of `T` in Hausdorff distance.

Full action convergence is stronger than needed. Even one-sided
containment of the approximating one-profiles in vanishing neighborhoods
of the target one-profile is enough for the required upper bound on
`E_quad(C_N)`. This is an all-order flat-sign recovery theorem, not merely
an arbitrary finite-matrix approximation theorem for bounded operators.

To check its implication, use the already proved same-order spectral
regularization of a subsequence approaching the original liminf. For
each fixed deletion/refill parameter `epsilon`, it yields signings at
those same subsequence orders, with a fixed operator cap and with
normalized optimum at most the original liminf plus `O(sqrt(epsilon))`.
The refill restores each selected order; it does not itself supply all
large orders. Extract an action-convergent subsequence
of these signings, with limit `T_epsilon`. The recovery statement and
Section 2 would produce signings at every large order with optimum at
most `E_quad(T_epsilon)+o(1)`. Thus the original limsup is at most its
liminf plus `O(sqrt(epsilon))`; then let epsilon decrease to zero.

## 4. What the new Gaussian modules do and do not add

The marked-tree and weighted unmarked identities constrain observables
that every such limit must realize. The variance-normalization inequality
is uniform and gives a new legitimate nonlinear test of these limits.
The matrix-chaos theorem controls the operator norms of fixed polynomial
derivative matrices. None of these statements constructs the all-order
flat-sign matrices required in Section 3.

In particular, ordinary finite-rank approximation of a bounded operator
does not preserve entries of modulus `1/sqrt(N)`. Disjoint copying fails
the same constraint and changes the critical normalization. Regular
Hadamard tensoring respects entry magnitudes but may enlarge the action
profiles and increase the quadratic optimum; it is not a proved recovery
operation for the original limit. These are specific failures of the
available constructions, not a general impossibility theorem for action
limits or for original convergence.

No inspected primary source currently supplies the stated recovery
mechanism. The exact remaining constructive target is the one-profile
upper recovery in Section 3, with the flat symmetric sign constraint and
the `sqrt(N)` normalization retained.
