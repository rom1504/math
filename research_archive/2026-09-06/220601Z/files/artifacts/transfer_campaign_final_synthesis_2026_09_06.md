# Renewed six-hour campaign: final mathematical synthesis

2026-09-06. Director synthesis of independently reconstructed proofs and
reproducible finite checks. The campaign started at 15:15 UTC; the roughly
16:05--17:02 service interruption is excluded from active time. Its closing
checkpoint is around 22:15 UTC. Preservation runs occurred during research,
not in place of it. This synthesis is not itself a proof dependency.

## 1. Original question: still open

```math
0.4333221116640807\le\liminf M_n/n^{3/2}
\le\limsup M_n/n^{3/2}<0.499432211<1/2.
```

The upper construction and its exact certificate were reconstructed. A
small additional entropy term improves its expression to approximately
0.4994322102336964. This excludes convergence to 1/2, not convergence to
another value. Neither convergence nor nonconvergence was established.
The exact order-15/16 witnesses independently replayed here give only
M15<=27 and M16<=30; external exhaustive lower catalogues were not replayed.

## 2. Strongest original-signing structural results

**Same-order balanced selection.** At every sufficiently large order there
is an actual hollow signing B with

```math
M_n\le Q(B)\le M_n+O(n^{5/4}),\qquad
|\max H_B-\max(-H_B)|=O(n).
```

The [construction](transfer_director_same_order_orientation_repair_2026_09_06.md)
uses a simultaneous diagonal majorant, low-field vertices of an actual
ground state, and an exact same-order block replacement. Its iterative
constants and finite inequalities were independently checked. It removes
orientation imbalance for a selectable near-minimizing sequence, but does
not supply isotropy, bounded operator norm, or transfer between scales.
Balanced near-minimizers can still have isotropic-response slack of order
n^(4/3); that is subleading, not a macroscopic impossibility theorem.

**Full-parent source approximation at energy scale.** For
B=A/sqrt(n-1), Lambda(B)<=C and the actual first-marked pair (G,Y), a fixed
Gaussian-L2 polynomial approximation with error delta satisfies

```math
\limsup_n {1\over n}\mathbb E\|B[f(G,Y)-p(G,Y)]\|_1
\le4K_G C\delta.
```

The [proof](transfer_fresh_untruncated_local_rate_and_energy_transport_2026_09_06.md)
uses uniform-in-root source approximation, not average row approximation.
An explicit local smooth-test bound is
`400||D3 phi||/sqrt(n-1)+32||D2 phi|| ||B||op/sqrt(n-1)`.
The archived elementary spectral bound gives ||B||op=O_C(n^(1/4)), yielding
uniform rectangle error O_C(n^(-1/12)). The qualitative local Gaussian law
was already proved in the archive; the quantitative calculation and cap-only
energy passage are the new conclusions. This is NOT transported L2 or nuclear
covariance control, nor a Gaussian law for the returned field BF.

## 3. Actual feedback, with the scope kept intact

For every fixed odd ternary first-marked f(G,Y) with both positive occupied
and hole masses, and every fixed operator cap L, the actual signed-feedback
argument gives a strict positive gain gamma(f,L) over its first response.
Its [proof](transfer_seed_high_degree_actual_feedback_probe_2026_09_06.md)
retains FULL BC and the joint literal regressed query, rather than inferring
energy from a covariance sign or independent marginal limits.

The dependencies are: exact finite source cuts and collision surgery;
mixed nuclear comparison; the signed finite-rank Gaussian return identity;
open-mark Hall contractions; joint Boolean Stein comparison; and exact
feasible spin flips in both energy orientations. Fixed degrees and cutoffs
precede the matrix-order limit. Two independent readers reconstructed the
mixed and endpoint steps, and the director checked their normalization.

Compact fixed-complexity rectangle families have a uniform finite probe
catalog and gain. However the PARTICULAR conservative gain now certified
never repays the stated principal-deletion loss, for any deletion fraction.
This is a limitation of that guarantee, not an upper bound on true feedback.

The deterministic endpoint itself can be made operator-free: instability
d gives gain d^2/(4C). The [cap-only moment theorem](transfer_fresh_untruncated_endpoint_and_row_moment_tests_2026_09_06.md)
also supplies simultaneous moment-capped roots for the entire old query
list and actual BF,BC without deleting the signing. What is not proved is
that useful return covariance survives on those roots together with the
necessary joint regressed-query comparison.

## 4. What survived the richer-history audit

Several finite lemmas are proved, but they do not yet make a rich-history
or untruncated feedback theorem:

- [Sparse/dense nuclear dual splitting](transfer_adversary_bounded_response_nuclear_dual_split_2026_09_06.md)
  upgrades two DIFFERENT local comparisons under explicit uniform cut and
  alias hypotheses. Its limit order is n, dual threshold, root cutoff.
- [Bounded trigonometric flat transport](transfer_seed_trigonometric_high_degree_transport_2026_09_06.md)
  supplies small high-original-degree cuts without polynomial density in
  the coherent law. Exact cover identities and all-slot derivative bounds
  retain overlapping labels. The primary matrix-concentration and graph-norm
  ingredients were read and mapped, not imported by analogy.
- [Local bounded-coefficient contractions](transfer_seed_bounded_local_full_contractions_2026_09_06.md)
  and weighted Sobolev tails permit fixed-degree local product surgery.
  Local smallness does not automatically control two open output roots.
- [Cubic Boolean chain cancellation](transfer_director_boolean_cubic_chain_nuclear_remainder_2026_09_06.md)
  gives nuclear O(n sqrt(kappa)) CHAIN remainder under vanishing aggregate
  coordinate influence. Actual random Stein brackets remain in its formula.

Three falsifiers prevent overextension. A bounded continuous response of
X^3 is orthogonal to every polynomial moment. An actual bounded-operator
Sylvester-twin signing has exactly half marks and holes but an EXACT degree-
one rich response, so holes alone cannot provide a high-degree probe.
A colored Boolean alias example has valid two-root local limits but a
leading nuclear discrepancy; the same-source alias condition is essential.
Their [moment](transfer_director_continuous_moment_approximation_counterexample_2026_09_06.md),
[actual-signing](transfer_director_rich_holes_without_high_degree_2026_09_06.md),
and [alias](transfer_adversary_nuclear_alias_counterexample_2026_09_06.md)
proofs preserve those distinct scopes.

## 5. The actual convergence gap and next direction

Let ell=liminf M_n/n^(3/2). No proved construction takes selectable seeds
A_j with normalized cap tending to ell and produces all sufficiently large
orders with limiting cap at most Q(A_j)/|A_j|^(3/2)+epsilon_j, epsilon_j->0.
Relatively dense target orders would suffice by principal restriction.
Target order tends to infinity FIRST, then j. The current upper ensemble
has not been shown asymptotically optimal or seed-faithful.

Full weaving, even tensor powers and several random selectors have genuine
loss theorems. They do not rule out all structured transfer. Pointwise
powered child-cap recurrences have counterexamples; the corresponding
original-M comparable-split recurrence remains unproved and unfalsified.

The next defensible bounded task is to test whether cap-only energy
approximation and moment-capped roots can retain quantitative SIGNED feedback
on actual optimizing signings without principal deletion. State the needed
covariance retention and joint query comparison before any computation.
It would improve our access to actual minimizers, but would not by itself
prove convergence. A convergence campaign must additionally deliver the
liminf-to-all-order comparison above or another explicit limit mechanism.

## 6. Preservation and verification boundaries

Canonical proofs and certificates remain in their original directories.
The dated research archive preserves unfinished work, failed solvers,
counterexamples, source programs, seeds, raw reports and informative outputs.
Environments, reviewed reproducible products and pinned downloaded software
dependencies are excluded with manifests. No research was deleted to clean
Git status. Two unavailable historical modules are documented as missing;
neither is a dependency of the present rigorous interval.

The final checkpoint replay records exact integer/rational finite tests
separately from floating trigonometric or Gaussian diagnostics. Finite
regressions do not certify an asymptotic theorem. Current proof and replay
dependencies were checked for undocumented ignored inputs; the automatic
literal-path scan is not represented as a proof of arbitrary dynamic code.
