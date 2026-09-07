# Six-hour favorable-flatification campaign: final mathematical synthesis

2026-09-07, campaign begun 05:27:43 UTC, closing audit through 11:28 UTC.
Status: **convergence and nonconvergence remain open**. Preservation is not
verification; the distinctions below are part of the claims.

## 1. Verified original-problem progress

The independently reconstructed interval is

```math
0.4333221116640807\le\liminf_n M_n/n^{3/2}
\le\limsup_n M_n/n^{3/2}<0.493608094.
```

This campaign improved the all-order upper bound from the reported
0.494515125. Its exact upper expression is

```math
U=\frac{97/20+(24/25)\log2-5151/6250}
{(97/10)\sqrt{24/25}}
=0.4936080935887486527\ldots.
```

See [upper proof](flatify_independent_2026_09_07_ternary_upper_proof.md)
and [independent reconstruction](flatify_adversary_2026_09_07_ternary_upper_reconstruction.md).
Dependencies: precision-Schur conditional-variance supersolution; direct
Gaussian-boundary stopping; exact ternary interval certificate; full-spin
and both-polarity counting; terminal-uniform Fock estimate; H2/H12 all-order
realization. The director replayed the exact certificate, and a separate
implementation checked its endpoint. Fix the desired margin, then finite
depth, then send order to infinity, then remove the margin. This is not
an assertion that the recursive ensemble is optimal.

The lower-bound dependencies were reconstructed earlier in the campaign;
see the standalone marked-response and width audits in the preceding
checkpoint map. No new asymptotic result depends on unrerun small-order
solver lower certificates. Decimal improvement alone does not prove a limit.

## 2. Strongest new construction and structural input

### Actual all-order subspace bridges — proved

For arbitrary subspaces with orthogonal projectors P,Q and total rank r,
there is a full n-by-n sign bridge C such that, simultaneously for all
Boolean x,y,

```math
|x^TCy|\le\sqrt n\|(I-P)x\|_2\|(I-Q)y\|_2
+O\left(n^{3/2}(r\log n/n)^{1/4}+n^{41/40}\right).
```

Thus any r=O(n^(1-delta)) has a power-saving error, up to a logarithm.
[Proof](flatify_adversary_2026_09_07_random_gauge_projector_localization.md),
especially its final diagonal-contraction refinement;
[independent audit](flatify_construct_2026_09_07_nearlinear_projector_audit.md).
Random Hadamard gauges bound projector entries, diagonal contraction makes
the projected matrix entry-feasible, its squared-mass deficit pays uniform
sign rounding, and nearby Hadamard orders pay the displayed all-order error.
The earlier nuclear-deficit proof gives an alternative
O(n^(5/4)sqrt(r)) error. Neither construction assumes spectral flatness.

### Actual-input fourth-moment rigidity — proved

For a real rectangular matrix A with maximum entry magnitude M,

```math
\|A\|_{S_4}^2\le G(A)M\le K_G\beta(A)M.
```

Here G is the vector Grothendieck relaxation and beta is the bilinear
Boolean norm. In particular Q(A)<=C n^(3/2) implies
tr(A^4)<=16 K_G^2 C^2 n^3. The proof uses SDP factorization and a
gamma_2 factorization of AA^T A, not a presumed spectral law.
[Proof and sharpness scope](flatify_construct_2026_09_07_spectral_fourth_moment.md).
A planted clique in an exact minimizer gives additive near-minimizers
whose higher normalized Schatten moments diverge for every exponent>4.
This is not a theorem about every exact minimizer's largest eigenvalue.

### Why a scalar residual envelope is insufficient — proved, scoped

For every bounded-cap actual signing and EVERY rank-o(n) projection,
near-orthogonal Boolean states retain energy half-width at least

```math
\left[\frac23\sqrt{\frac{2}{3\pi}}-o(1)\right]n^{3/2}.
```

[Proof](flatify_independent_2026_09_07_rademacher_residual_width.md).
An independent/random split, greedy signs, scalar Berry--Esseen and the
new fourth-moment bound give this uniformly over the projection. Hence
separately maximizing child energy plus the residual-norm envelope costs
at least (1+2c_R-o(1))n^(3/2), above the desired equal-child target.
This does NOT lower-bound the actual constructed parent: its bridge may
be much smaller than that envelope on the offending child states.

## 3. New full-maximum transfer: a paid Gaussian comparison

For any Gaussian correlation matrix epsilon I<=R<=KI, let C=sign(G_R)
and let Y be Gaussian with the EXACT covariance Cov(C)=(2/pi)arcsin[R].
For arbitrary bounded features |u_i(omega)|<=1, arbitrary deterministic
offsets b(omega), and ANY finite configuration set,

```math
\left|\mathbb E\log\sum_\omega e^{b_\omega+\lambda C\cdot u_\omega}
-\mathbb E\log\sum_\omega e^{b_\omega+\lambda Y\cdot u_\omega}\right|
\le C K^2\epsilon^{-2}|\lambda|^3d\log^{3/2}(d+1).
```

[Full proof](flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md),
[constructive audit](flatify_construct_2026_09_07_quenched_gaussian_universality_audit.md),
[adversarial reconstruction](flatify_adversary_2026_09_07_quenched_gaussian_sign_audit.md).
The director checked the conditional replacement, exact covariance,
block-local Stein kernel, two-star contraction, arbitrary Gibbs dependence,
gap dependence, and all soft-max normalizations. The proof is uniform in
the number of configurations, not an annealed pairwise union estimate.

For n-by-n bridges, fixed epsilon,K, this gives expected absolute PARENT
cap error O(n^(4/3)sqrt(log n)), with the same actual child energies on
both sides. Thus sign realization costs a power-saving error once a
suitable Gaussian-process parent bound is proved. It does not prove that
bound. The earlier one-source MGF theorem remains valid but is weaker
for this application; its central-shell first-moment obstruction is not
an obstruction to the new quenched transfer.

For the explicit opposite-spectral law
R=I-rho A tensor D/(||A||op||D||op), no flatness assumption is needed.
Letting 1-rho=n^(-1/7) also pays approach to the endpoint GAUSSIAN process
with O(n^(10/7)sqrt(log n)) expected-cap error. This is not a theorem
about the endpoint singular sign law itself.

## 4. Sharp falsifier on actual near-minimizers

If p=||A||op||D||op and p/n tends to infinity, the preceding particular
opposite-spectral law becomes iid-like in its expected parent cap:
uniform cross-covariance error O(n^3/p) gives cap error O(n^2/sqrt(p)).
The sign/Gaussian transfer and a conditional-variance concentration proof
give, with probability tending to one, normalized parent cap at least
P_SK/sqrt(2)>0.53. The independently certified heat-control lower bound
P_SK>3/4 suffices. This is worse than the current original upper bound.

Such actual additive near-minimizing child sequences exist: plant a +1
clique of size floor(n^(3/4)/log n) in an exact minimizer. Its edit cost
is o(n^(3/2)), while p/n diverges for two copies.
[Proof](flatify_independent_2026_09_07_opposite_spectral_nearmin_law_obstruction.md).
This rejects typical sampling/expectation-based existence for THIS law
on THESE near-minimizers. It excludes neither rare favorable outputs,
selectable exact minimizers, another covariance, nor global old-edge changes.

Finite checks are kept separate: optimal order-12 children with a fixed
Paley-12 bridge row/column ordering have exact minimum cap62 over all
row/column sign phases and both relative child polarities. Full 24-spin
replay verifies witnesses. The random permutation search is heuristic;
the exact claim does not cover all bridge permutations or all bridges.
Order-8 opposite-spectral sampling is exploratory, not an optimality proof.

## 5. Exact remaining convergence obligation

For every comparable sufficiently large split N=m+l, find selectable
actual optimal children and a hollow full signing C_N with

```math
Q(C_N)\le Q\!\left(\operatorname{diag}\left(
\sqrt{\frac{N-1}{m-1}}A_m,
\sqrt{\frac{N-1}{l-1}}A_l\right)\right)
+O(N^{3/2-\delta})
```

for a fixed delta>0. Old internal edges may all change. Complementing
children aligns their positive maxima. With u_n=M_n/sqrt(n-1), the
inequality yields u_N<=u_m+u_l+O(N^(1-delta)); the audited balanced-tree
argument then proves convergence. Bare o(N^(3/2)) without a summable
rate is not silently substituted. Global child reversal always gives
max_epsilon |a+b+epsilon c|=|a+b|+|c|.

This obligation is **not proved**. The Gaussian transfer removes a real
realization/entropy payment for a class of correlated bridges, but leaves
their entire joint Gaussian parent bound unpaid. The subspace construction
removes another actual realization step; its scalar residual certificate
cannot finish the argument. Neither is presented as a completed strict
reduction to an evidently easier optimization problem.

## 6. Director judgment and next direction

Strongest original-problem theorem: the strict all-order upper U above.
Strongest new mechanism: full quenched covariance transfer, usable with
arbitrary child offsets, combined with actual low-dimensional annihilation.
The work has gone beyond formal reformulation but has not established
near-optimality preservation or a convergence recurrence.

The single next target should be a **quantitative favorable Gaussian
parent comparison for actual selectable optimal children**, with an
explicit admissible covariance and a power-saving total cap payment.
The new transfer theorem supplies a genuine sign realization if such a
comparison is achieved. A p-normalized covariance must not be justified
by bounded cap alone: the actual-near-minimizer falsifier forbids that.
Global dense operations remain available if that expectation mechanism
cannot beat its iid floor. No assertion of convergence to a chosen
constant, and no nonconvergence conclusion, is supported.

All intermediate proofs, unsuccessful mechanisms, scripts, seeds,
outputs and independent reports are retained in canonical artifacts and
dated research archives. They need not be polished to be preserved.
