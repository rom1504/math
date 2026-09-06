# Rich feedback extension: verified analytic bridge and exact remaining test

2026-09-06, final extension checkpoint. This note records the outcome of
the bounded-trigonometric branch. It does NOT claim a new rich-frame cap
theorem, a larger universal lower constant, or convergence of M_n/n^(3/2).

## 1. What is already an actual cap theorem

`transfer_seed_high_degree_actual_feedback_probe_2026_09_06.md`, audited
in `transfer_adversary_high_degree_actual_feedback_audit_2026_09_06.md`,
proves a strict gain over the first marked response for each fixed
operator cap and fixed nontrivial ternary response of the first Gaussian
old variables. In the notation of that proof,

```
liminf [Q(A)/(n sqrt(n-1)) - average E H|BF|] >= Delta(f,L)>0.
```

Its exact literal returned-query comparison, not merely a positive
correlation with a probe, is essential. The explicit certified Delta
does not currently pay the unrestricted spectral-core loss; see
`transfer_fresh_uniform_first_marked_family_and_core_barrier_2026_09_06.md`.

## 2. Verified new analytic ingredients

The following are finite, independently checked statements.

* `transfer_seed_bounded_trigonometric_global_cuts_2026_09_06.md`:
  the entrywise exponential has a quadratic operator bound, and exact
  finite Boolean cover identities preserve all fixed global cuts.
* `transfer_seed_trigonometric_high_degree_transport_2026_09_06.md`:
  if a finite old polynomial frame has degree at most M and polylogarithmic
  global cuts, every fixed q>M component of B exp(itW) has proper cuts
  O(max|B_ij| polylog(n)). The independent full audit is
  `transfer_adversary_trigonometric_transport_audit_2026_09_06.md`.
* `transfer_seed_bounded_local_full_contractions_2026_09_06.md`:
  bounded smooth composition preserves LOCAL primitive mixed Gamma
  smallness; fixed original degree isolates each full contraction.
  For fixed trig coefficients, every original Sobolev degree moment
  is bounded, and a polynomially weighted multiplication bound permits
  ordered product L2 approximation. No polynomial density in the
  numerical values of W is used.
* `transfer_adversary_bounded_response_nuclear_dual_split_2026_09_06.md`:
  the sparse/dense dual split proves its stated mixed nuclear comparison
  for a finite frame with uniformly bounded global cuts and its explicit
  no-alias hypothesis. My independent scope audit is
  `transfer_seed_bounded_response_dual_split_audit_2026_09_06.md`.

These statements address different norms. They must not be substituted
for one another solely because all are called contraction estimates.

## 3. The old rich-frame primitive Gamma hypothesis is available

The exact relevant dependency is
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`,
Sections 10--11. Its fully injective forest argument proves deterministic
proper/full kernel-contraction smallness against exceptional first-chaos
transports BX_T=QV_T+error. In particular its root-hit argument proves
E||J||F^2=o(n) before bounded transports and diagonal-root projection.

Those deterministic contractions also control the Boolean Gamma product:
the exact same diagonal-free coefficient tensors occur, followed by
symmetrization and free-label distinctness projections. The finite-degree
Boolean product formula therefore yields averaged L2 primitive Gamma
smallness. One must use this explicit contraction conclusion, NOT infer
Gamma-square control from the old Section 12 linearly tested input
replacement. The adversarial agent independently checked this mapping.

For averaged maximum influences delta_i and primitive Gamma norms
gamma_i, the new bounded chain lemma gives the rate-free estimate

```
average ||Gamma(X_i,A(W_i))||_2^2
 <= C[average gamma_i^2 + sqrt(average delta_i)].
```

This matches the old proof's averaged collision/influence errors.
No uniform logarithmic rate is silently added. Scope remains a fixed
finite marked-tree/fully injective forest frame and a fixed operator cap.

## 4. The concrete remaining returned-field check

For a selected source degree P and odd k>=3, the probe is an exact
degree-p noise with p=kP. To pass from its positive full-return
correlation to actual spin ascent, one needs the literal joint query
comparison for

```
(old coherent W, fixed old nonlinear channels, BC-beta_i E_i).
```

At every fixed larger original degree d>p, this requires the full
probe-into-return contraction to vanish in averaged squared Hilbert norm.
Equal-degree covariance is retained by beta, not discarded. A sufficient
way to establish the difficult larger-degree case is to keep both output
root axes open until one has an o(sqrt(n)) Frobenius tensor estimate;
bounded root transports and the diagonal-root projection can then be
applied safely. Entrywise O(n^(-1/2)) is not enough: n^2 such entries
have Frobenius norm of order sqrt(n).

The natural decomposition, at fixed response approximation stages, is

```
BC=B c0(W)+B[(A(W)-a)Z+R_(>=2)]+B D_a B r,
a_i=E A(W_i).
```

The new high-degree trig theorem handles Bc0. The local Gamma/Sobolev
lemmas remove the local coefficient/noise product-approximation obstacle
for the middle terms. Their finite positive-degree factors can use the
global-cut/two-factor transport bounds. In the final term a is a
DETERMINISTIC diagonal coefficient; an open k-branch Hall argument
against the literal old source, keeping both roots open, is the intended
completion. The full chain through every fixed approximation, collision
error and root cap has not been reconstructed end-to-end here. Hence
this checkpoint leaves the rich returned-field theorem open.

Late bounded addition: the finite-stage deterministic mean-linear branch
has now been checked separately in
`transfer_seed_rich_mean_linear_open_hall_2026_09_06.md`. The old open-Hall
proof extends to P>M: k left branches admit k distinct proper/influence
merges, giving a two-root open Frobenius bound O(n epsilon_n^k).
It retains an explicit uniform main-kernel rate and separate averaged
error passage. This removes that branch's fixed-stage algebraic gap,
not the complete ordered bounded-response assembly.

Treating BF itself as a finite polynomial query of degree D and selecting
p>D is not a shortcut: the cutoff D needed for a sufficiently accurate
bounded-response approximation can depend on and exceed the already
selected probe degree p. A valid argument must fix its cutoffs in order,
not assume this circular degree choice can be solved.

## 5. Independent source-tail and alias limitations

Even completing Section 4 would not justify a universal claim for every
rich ternary policy. The actual twin-Hadamard example in
`transfer_director_rich_holes_without_high_degree_2026_09_06.md` has
positive holes but an exactly degree-one response. Thus a positive
original source block P above the fixed old primitive degrees must be
proved for the particular policy or made an explicit hypothesis.

Nor can literal same-old-source structure be replaced by generic good
cuts. The exact later-channel alias in
`transfer_adversary_nuclear_alias_counterexample_2026_09_06.md` has
controlled cuts and small entrywise local errors but a positive nuclear
gap. It violates the required no-alias condition and does not contradict
the actual fixed-old-frame theorem.

## 6. Verification and preservation

The trigonometric checker passes 192 matrix bounds, 3200 scalar cover
identities and 96 derivative-matrix bounds. The new local Sobolev/product
checker passes 90 direct identities, 90 multiplication tests and 8726
factorial-coefficient tests. Both scripts are self-contained in
`computations/`; neither reads tmp inputs or writes output artifacts.
The mathematical scopes above do not rely on these finite tests as
proofs of their asymptotic conclusions. No commits or archival moves
were performed by this branch.
