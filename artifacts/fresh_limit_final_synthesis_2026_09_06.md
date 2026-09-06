# Six-hour campaign: strongest result and unresolved original question

Closing audit for the campaign authorized on 2026-09-05 at 18:29 UTC.
Substantive work and independent verification continued through 00:29 UTC
on 2026-09-06. The proofs, certificates, and final state are committed.
The original convergence/nonconvergence question is **not solved**.

## Strongest verified original-problem result

For symmetric hollow sign matrices, with the original absolute same-spin
normalization, the independently audited proof and exact interval replay give

```math
0.4306581794055286\le\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}\le\frac12.
```

The lower endpoint is a universal analytic result, not a finite-order fit.
The exact 21-anchor Gaussian certificate lies in

```
[.430658179405528602724053804634711026327173238336325190455581,
 .430658179405528602724053804634711026327173238336325190455840].
```

The final replay of both this certificate and its nonzero cubic coefficient
passed. The proof chain can be kept short:

1. Classical diagonal Grothendieck regularization retains a principal
   `(1-epsilon)n`-vertex signing with a fixed normalized operator bound.
2. Exact injective-tree moments and the polynomial paired-energy identity
   realize a Gaussian creation isometry on each fixed finite family.
3. Bounded-operator L2 approximation realizes the finite-anchor certificate.
   No unbounded-operator Sobolev extension is needed for this numerical bound.
4. Principal monotonicity and then `epsilon -> 0` give the displayed
   universal lower endpoint.

The two feasible means are evaluated on the **same** hollow signing.
Their exact half energy difference, not presumed child-bridge cancellation,
gives the bound. The original bridge-reversal identity remains respected.

[Certificate theorem](fresh_finite_anchor_fixed_point_2026_09_05.md),
[minimal dependency audit](fresh_final_lower_dependency_algebra_audit_2026_09_05.md).

## Stronger result, with a deliberately unevaluated increment

Weighted odd Hermite transport and a Schur mean-standard-deviation inequality
prove `liminf M_n/n^(3/2)>C_scalar`, where `C_scalar` is the supremum over
the entire old scalar central-mask hierarchy. In particular the liminf is
strictly above the exact finite-anchor value. The extra positive increment
has not been evaluated. Every finite construction is fixed before matrix
order grows; no generic AMP universality is assumed.

[Uniform escape](fresh_uniform_scalar_hierarchy_escape_2026_09_05.md) and
[independent weighted audit](fresh_weighted_projection_final_adversarial_audit_2026_09_05.md).

A final elementary inequality gives a degree-417 nonlinear direction for
every arbitrary marked mask with certificate at least .43. It does not
by itself provide the local feasible slack needed to improve every such
mask. A proposed uniform intrinsic-innovation shortcut is explicitly false,
even at high value and with minimal finite Gaussian input spaces.
[Theorem](fresh_general_mask_nonlinear_direction_2026_09_06.md),
[falsifier and exact remaining slack question](fresh_arbitrary_mask_zero_strip_obstruction_2026_09_05.md).

## Strongest structural falsifier and corrected inherited inference

Exact order-12 certificates give two bounded-operator hollow signing
families with asymptotically equal squared operators and separated normalized
Boolean caps. Their present odd-channel one-root covariances also agree.
The certified gap is at least `437/(2000*12^(3/2))`. These are not
near-minimizers and do not prove nonconvergence. Cross-root signed action
is not asserted equal. [Exact construction](fresh_cosquare12_scalable_gap_2026_09_05.md).

Schmidt's primary Walsh theorem settles the H2 regularized seed at sqrt(2).
An inherited inference from a stagnant greedy trajectory was too strong:
that trajectory did not describe the global optimum.
[Checked primary mapping](fresh_schmidt_odd_walsh_regularization_2026_09_05.md).

## Precise gap to the original convergence question

No matching all-order upper mechanism was proved. One precise sufficient
missing statement, using already verified constructions, is the following.
Put `c_*=liminf M_n/n^(3/2)` and let R(B) be the same-spin regularized
cap of a finite full symmetric sign seed B under the specified regular
order-4/order-144 Hadamard semigroup. Prove

```
For every epsilon>0, there is a finite full sign seed B of order k
with R(B)/k^(3/2) <= c_*+epsilon.
```

The proved all-order full-seed lifting theorem would then give
`limsup M_n/n^(3/2)<=c_*`, establishing convergence. This seed-recovery
statement is **not proved**, not asserted necessary for every conceivable
convergence proof, and not demonstrated to be easier than the original
problem. The exact majorant floor and solved H2 seed do not establish it.

Alternatively, a bounded-op liminf action object needs all-order flat-sign
upper recovery. The square-state falsifier blocks a general replacement of
that signed-action obligation by squared-operator data alone. No theorem
here rules out special recovery for actual near-minimizers.

Thus the campaign improves the rigorous interval and supplies new checked
mathematics, but provides neither a convergence proof nor two separated
infinite subsequences of the minimizing values. Solver timeouts, finite
plateaus, and surrogate counterexamples have not been promoted to either.
