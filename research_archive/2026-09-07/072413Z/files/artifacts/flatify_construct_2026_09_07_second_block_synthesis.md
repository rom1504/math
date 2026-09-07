# Constructive track: second sustained research block

Date: 2026-09-07; continuation after the 06:28 first-block handoff.
Scope: leave MUB lifts and seek an actual global old-edge surgery operation.
No cross-order flatification or convergence theorem was obtained.

## Main positive operation

The constructive track first developed uniform signed finite-rank surgery
and optimizer-derived covariance stationarity. Independent adversarial
auditing simplified the Schur majorant to rank one, and then supplied the
stronger target-contraction operation. The latter is the definitive proof:

    artifacts/flatify_adversary_2026_09_07_target_contraction_surgery.md

For any actual hollow full signing A and symmetric Delta with nuclear norm
L, there is an actual full signing A' with

    Q(A')<=Q(A-offdiag(Delta))+2sqrt(nL(n+2)log2)
                                      +4(n+2)log2/3.

The coefficient mean is S(A-offdiag(Delta))S, where
S_ii=(1+|Delta|_ii)^(-1/2). This mean is inside the edge coefficient cube,
and contraction of the TARGET gives its cap bound exactly. Independent
rounding has variance at most2nL. Thus L=o(n) gives a lower-order error;
fixed-rank, order-sqrt(n) perturbations have O(n^(5/4)) error. No cap bound,
operator norm, incoherence, or favorable-ramp assumption is needed.

More generally it suffices to supply nonnegative ell_i with
|Delta_ij|<=sqrt(ell_i ell_j), and use L=sum ell_i. Nuclear leverage and a
uniform entrywise bound are two choices. The weighted same-amplitude
extension, including forbidden zeros, also passed our independent check.

## Optimizer consequence

Actual optimality is applied only AFTER sign recovery. Minimax then gives
signed near-ground covariance balance on each selected finite feature
space. The independent track extended the feature budget to a convex
nuclear ball and a single full-space covariance law. None of these laws
forces equal polarity masses, isotropy, or midpoint balance: a broad
one-polarity near-ground ensemble can escape every small feature space.

Original derivation and failures, preserved with subsequent improvements:

    artifacts/flatify_construct_2026_09_07_finite_rank_surgery_stationarity.md

Independent audit of original proof:

    artifacts/flatify_adversary_2026_09_07_finite_rank_stationarity_audit.md

Nuclear-ball extension is owned by the independent track:

    artifacts/flatify_independent_2026_09_07_nuclear_budget_stationarity.md

## Entrywise operation

Independently of the rank route, for |G_ij|<=1 and t=theta/sqrt(n), use
mean (A-tG)/(1+t). It has cap Q(A-tG)/(1+t), and rounding error O(n^(5/4))
for fixed theta. Minimax gives a full signed covariance law with

    E[q-sigma h]+(theta/n²)sum_(i<j)|E[sigma x_i x_j]|=O(n^(-1/4)).

The precise inequality includes near-minimality excess and explicit
constants. This is also a special case of the supplied-envelope operation.

    artifacts/flatify_construct_2026_09_07_entrywise_surgery_stationarity.md

    artifacts/flatify_adversary_2026_09_07_entrywise_surgery_audit.md

Constant-density versions do not close composition: the rounding entropy
constant is too large, and regimes giving small covariance cost cease to
force concentration near the ground states. No operative cap recurrence
was inferred from the necessary condition.

## Higher-rank favorable-direction investigation

An explicit bounded-cap actual signing sequence shows that a favorable
rank-o(n), operator-bounded direction need not itself have small entries:

    artifacts/flatify_construct_2026_09_07_favorable_direction_coherence.md

The construction mixes a weak rank-one positive bias into a symmetric
Hadamard signing and appends a coherent block-projector penalty. The given
direction lowers cap by a fixed leading amount while some correction
entries diverge. It is NOT a near-minimizer example and has an obvious
alternative rank-one descent. Thus a replacement-descent theorem remains
open; merely inferring regularity of every favorable direction is false.

For op-bounded rank-o(n) directions, large correction entries occupy only
o(n²) edges, while the small-entry remainder is recoverable. The unresolved
step is cap-safe removal or replacement of the exceptional correction.
Sparse edge count alone is insufficient. In coherent block models, a large
coarse diagonal can genuinely center an energy interval on a block-reversal
orbit, so discarding it runs directly into the midpoint problem.

## Explicit non-closures

1. Direct comparison with the weighted two-child target has Frobenius
   discrepancy Omega(n), hence cannot have the o(n) envelope/nuclear budget
   required by this operation.
2. Adding one vertex and deleting its bridge has nuclear cost Theta(sqrt(n));
   generic recovery costs O(n^(5/4)), worse than the trivial O(n) insertion
   bound. Even directly rounding the feasible target diag(A,0) only gives
   O(n). No summable insertion estimate follows.
3. A rank-one completion-of-squares penalty suppresses one energy polarity
   while enlarging the other. Signed covariance balance does not reverse
   that sign on the two extremal landscapes.
4. None of the MUB obstruction proofs from the first block is used to claim
   that this new operation cannot eventually contribute to flatification.

All scripts and drafts remain in the repository; no Git commit was made by
this track. This checkpoint records a genuine stronger actual sign operation
and its audited optimizer consequence, not a claimed resolution of the
requested cross-order construction.
