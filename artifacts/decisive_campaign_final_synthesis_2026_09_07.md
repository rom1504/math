# Six-hour original-limit campaign: closing synthesis

Campaign ran from 2026-09-06 22:44:53 UTC through 2026-09-07 04:45 UTC.
Closing verification and research preservation are complete. This document records results,
not a proof of convergence. The chronological checkpoints are in
[the campaign index](decisive_limit_campaign_2026_09_06.md) and ledger 10.151.

## 1. Strongest result about the original problem

For hollow symmetric full sign matrices, put
`H_A(x)=sum_(i<j) A_ij x_i x_j`, `Q(A)=max_x|H_A(x)|`, and `M_n=min_A Q(A)`.
The current independently reconstructed interval is

```math
0.4333221116640807
\le \liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
<0.494515125.
```

The new upper constant is

```math
u_*={4+(31/32)\log2+g_4(1)\over8\sqrt{31/32}},
\qquad
u_*\le {7787631971809\over15748015748016}<0.494515125.
```

Its value is approximately 0.494515124725174. This improves the previous
0.499432211 upper endpoint. **Neither convergence nor nonconvergence is
proved.** Convergence to 1/2 is excluded. No numerical value for a putative
limit is inferred from this particular construction.

### Minimal upper proof dependency map

1. **Precision-side Schur inequality.** Arithmetic/harmonic precision pivots,
   the information chain rule, and log-majorization prove `BE<=E` for the
   mean conditional-variance envelope. This removes a genuine analytic
   alignment obligation; it is not an assertion about the optimal signing.
2. **Direct Gaussian-boundary stopping.** Mutual-information drift, strict
   Gaussian equality, and weak-moment tail control bound deep finite Bellman
   iterates by E. The stopping budget is `Phi-G`. The stronger identity H=E
   is proved separately but is NOT a necessary upper-bound premise.
3. **Exact ternary Gaussian phase.** For p=31/32 and t=4, an analytic
   reproduction-support reduction plus 212505 outward-rounded rectangles
   proves `E_4(nu_p)=g_4(1)`. The full certificate was independently replayed
   several times; a separate rational calculation encloses the final number.
4. **Actual sign realization.** Full-spin counting includes both objective
   signs and all `2^N` spins. Independent fibre averages give `(E Z)^m`.
   The Gaussian orbital estimate is uniform in terminal orthogonal matrices;
   diagonal removal costs only O(N). H2/H12 terminal orders and principal
   restriction supply all sufficiently large orders, not just a subsequence.

The limit order is: fix the desired margin, choose finite depth, send matrix
order to infinity through the all-order construction, then remove the margin.
No depth-uniform CLT, negligible rare-spike assumption, prime-gap theorem,
or lossless arbitrary-seed transfer is used.

Read the [standalone proof](decisive_audit_standalone_direct_E_upper_2026_09_07.md),
[phase and final conversion](decisive_director_gaussian_phase_all_order_upper_2026_09_07.md),
and [fresh independent reconstruction](decisive_independent_upper_realization_audit_2026_09_07.md).

### Lower proof reconstruction

The existing lower endpoint was rebuilt from the actual marked two-Gaussian
response, tree-frame/covariance identities, nonlinear replacement estimates,
and spectral deletion. Its 256-cell policy certificate was replayed
byte-for-byte. Approximation parameters are fixed before matrix order and the
deletion cutoff is removed last. Prior positive verdicts were not substituted
for these dependencies. See the
[fresh lower audit](decisive_audit_fresh_full_lower_chain_2026_09_07.md).

## 2. Additional proved mathematics, with scope separated

- **Half-range and square bipartite lower bounds.** For
  `W(A)=(max H_A+max(-H_A))/2`, the same c*=0.4333221116640807 lower bound
  holds asymptotically for its minimum. On square m-by-m full sign matrices,
  `liminf min_C ||C||_(infinity->1)/m^(3/2)>=2c*>.8666442233281615`.
  The complete-bipartite support is explicitly retained in the proof; this
  is not a theorem about all weighted row-normalized matrices.
  [Width](decisive_audit_certified_minimum_width_lower_2026_09_07.md),
  [bipartite theorem](decisive_audit_bipartite_marked_response_lower_2026_09_07.md).
- **Actual cap controls quadratic fluctuations.** A cap bound, without
  spectral flatness, gives quantitative Gibbs variance for delocalized
  quadratic observables. Integration gives an extensive bridge gain for
  actual minimizing children at contracted temperatures. Temperature
  restoration is NOT paid by this result.
  [Proof](decisive_director_delocalized_quadratic_fluctuations_2026_09_07.md).
- **Exact actual-optimizer interpolation.** A log-cosh path splits the
  width-pressure derivative into a favorable PSD term and a signed imbalance
  of nonnegative edge-flip costs. Higher-order Taylor errors are sublinear;
  the leading integrated imbalance remains uncontrolled.
  [Identity](decisive_bridge_width_logcosh_interpolation_2026_09_07.md).
- **Conditional actual-sign midpoint repair.** A proportional quadratic
  ramp permits sparse sign surgery. General trace `r=o(sqrt(n))` costs
  `O(n^(4/3)r^(1/3))+e`; sign-compatible positive contractions permit
  `r=o(n)` with normalized error `e/n^(3/2)+O((r/n)^(1/4))+O(n^(-1/2))`
  and explicit probability feasibility. This is a real sign construction,
  not a weighted intermediate silently declared rounded. Such certificates
  have NOT been obtained for selected width minimizers.
  [Proof and precise hypotheses](decisive_audit_low_rank_midpoint_sign_surgery_2026_09_07.md).
- **Continuous approximation with quantified error.** Soft-flatness replaces
  the discrete value within `log2/(4tau)` uniformly in dimension; spectral
  regularization preserves leading outer entropy. The corresponding
  all-order variational limit is still open, not supplied by an ultrafilter
  random-matrix LDP.
  [Approximation](decisive_independent_soft_flatness_variational_bridge_2026_09_07.md).
- **A new consequence using BOTH endpoints.** On any actual original
  near-minimizing sequence, every admissible scalar-affine Gaussian law
  `Z~N(0,Id+sA)` has `limsup E|H_A(sign Z)|/n^(3/2)<.427688`, strictly
  below the universal width lower bound. Cap-only spectral control and a
  four-sign Taylor calculation give variance `O(n^(5/2))`; no spectral
  flatness is assumed. This excludes that particular mean-absolute rounding
  witness, not anisotropic/nonlinear Gaussian constructions or the signings.
  The previous upper endpoint did not suffice for this strict gap.
  [Theorem](decisive_director_affine_gaussian_near_minimizer_gap_2026_09_07.md).

## 3. Strong falsifiers, not universal impossibility claims

Prescribed dense Hadamard stabilization has a half-floor even for optimal
seeds, hence cannot be lossless at the now strict-subhalf optimum. This does
not rule out other lifts or nonlinear retention rules.

Generic bounded-cap temperature payment is false on an explicit full-sign
sequence with bounded normalized operator norm. Actual GLOBAL minimizer
payment is not falsified. Contracted-temperature scalar inequalities alone
admit convex abstract profiles with oscillating zero-temperature slopes;
these profiles are not asserted to be signing realizations.

Two-sided adaptive variance universality is false: sparse variance-one
amplitudes can force optimized normalized cap at least 2/pi-o(1), despite
the usual fixed-sign third-moment scale. This does not kill favorable
one-sided flatification with fixed bounded amplitudes.

Requiring almost complete projection of every high-energy state forces
almost full rank or an already vanishing midpoint. The weaker proportional
ramp avoids that circularity, but Steiner-frame full-sign examples show
that even it can require trace Omega(n) for bounded-cap parents. They are
NOT width minimizers. Exact finite skew-lift tests also fail the desired
factor on actual small width minimizers; no scalable orientation obstruction
was proved. Cut reversal turns that lift into a sum of absolute responses,
so it supplies no favorable signed cancellation.

Proofs and exact replay evidence are linked from the checkpoint ledger and
[director closing derivations](decisive_director_closing_derivations_2026_09_07.md).
Finite solver results, exhaustive certificates, and asymptotic proofs remain
separately labeled. The historical n=11,13 optimum lower solvers were not
rerun; none of the new asymptotic bounds depends on those finite values.

## 4. Exact remaining gap and next direction

The new recursion constructs excellent signings but does not propagate an
arbitrary near-liminf seed without a fixed leading loss. No theorem currently
transfers a liminf object to all large orders. The width approach additionally
needs BOTH integrated optimizer payment AND width-to-absolute recovery.

One precise sufficient original-value target is global sign replacement of
the two-block weighted matrix

```math
B=\operatorname{diag}\left(
\sqrt{\frac{N-1}{m-1}}A_m,
\sqrt{\frac{N-1}{n-1}}A_n\right),\quad N=m+n,
```

where A_m,A_n are selectable actual minimizing children at comparable large
orders. Every row of B has squared norm N-1. A full hollow sign replacement
with `Q(A_N)<=Q(B)+C N^(3/2-delta)`, uniformly for some delta>0, would give

```math
u_{m+n}\le u_m+u_n+O((m+n)^{1-\delta}),
\qquad u_k=M_k/\sqrt{k-1},
```

and the verified balanced-tree almost-subadditivity argument would force
convergence. See the
[flatification scope and exact implication](decisive_independent_favorable_flatification_scope_2026_09_07.md).
This is an OPEN constructive lemma, not a newly completed strict reduction.
Global changes to old coefficients are allowed; independent entry rounding
and a free bridge-cancellation claim are not its proof.

Recommended next work: one focused attempt to prove or scalably falsify this
fixed-oversaturation, selected-child replacement. Do not substitute more
conditional ramp statistics for it. The precision-Schur construction should
be preserved as a successful all-order upper-bound mechanism, not presumed
to be the asymptotically optimal ensemble.

## 5. Preservation

Canonical proofs and scripts remain in their existing locations. Unfinished
derivations, finite failures, parameters, informative outputs, and subagent
research are retained in dated `research_archive/` snapshots with hashes.
Reproducible build products, environments, caches and pinned external
dependencies are excluded under the reviewed manifests, not deleted.
The [final publication audit](../research_archive/publication_audit_2026_09_07.md)
records the exact snapshot and remaining files. The post-publication audit
matched all 1,824 remaining research originals to saved hashes. The ordinary
untracked scan executable is a reviewed rebuild product, not lost research.
