# Active research state

Evidence cutoff: ledger Section 10.148.4, resumed campaign in progress;
full nonlinear response theorem independently audited. Use the ledger for archive
comparison; use the linked proof files to reconstruct new statements.

## Exact problem and verified frontier

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,\quad
M_n=\min_A Q(A),
```

where A is hollow symmetric with off-diagonal signs. The objective is
convergence or nonconvergence of M_n/n^(3/2), not specifically the value 1/2.

```math
0.4314603928237005\le\liminf_n M_n/n^{3/2}
\le\limsup_n M_n/n^{3/2}\le1/2.
```

Recorded M_3,...,M_14 are (3,4,4,5,9,10,12,13,17,18,20,21).
The n=11,13 lower bounds are solver-certified infeasibility results.
The fresh verification script does not rerun those solvers.

The improved lower bound is analytic, with an exact rational interval
certificate for its numerical evaluation, not a finite-order extrapolation.
See the [full nonlinear theorem](artifacts/resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md),
[conditional-projection certificate](artifacts/resumed_response_conditional_v_full_center_2026_09_06.md),
and [independent numerical reconstruction](artifacts/resumed_director_third_checkpoint_audit_2026_09_06.md).
The upper is unchanged. The new theorem uses the whole nonlinear local Gaussian
response norm, not only selected edge coefficients. Nuclear covariance, actual
endpoint parity and full mixed contractions justify a fresh one-step channel;
neither a full BF conditional law nor an iterative refresh is assumed.

For N=binom(n,2), let C_n^+={(c+b_i+b_j)_(i<j)} in binary coordinates.
Then Q(a)=N-2d(a,C_n^+) and M_n=N-2rho(C_n^+). For n>=3 its dimension
is n. It is exactly RM(1,n) punctured to the weight-two slice.

## Strongest new quantitative module

For mu_r=E|epsilon_1+...+epsilon_r|,

```math
M_{n+r}\le M_n+M_r+
\min\left\{
n\mu_r+\sqrt{2nr(r-1)\log2},
r\mu_n+\sqrt{2nr(n-1)\log2}
\right\}.
```

All orders and exact minimizing children are covered. The associated bridge
minimum is n mu_r(1+O(sqrt(r/n))) whenever r=o(n). If r also diverges,

```math
M_{n+r}\le M_n+(\sqrt{2/\pi}+o(1))n\sqrt r.
```

This improves the archived random-completion coefficient, not its exponent
or the qualitative o(n) order-transfer principle. See
[the proof](artifacts/mesoscopic_completion_2026_09_05.md).

## Exact bridge identity and missing recurrence

```math
Q\!\begin{pmatrix}A&B\\B^{\mathsf T}&D\end{pmatrix}
=\max_{x,y}\bigl(|H_A(x)+H_D(y)|+|x^{\mathsf T}By|\bigr).
```

Large values of its two nonnegative terms would have to occur in different
places. A scalar bound on each term cannot capture that. The archived
separated certificate has an equal-split b-scale defect of at least
(0.218646...+o(1))n, for b_n=M_n^(2/3).

A sufficient open target is a bridge for exact minimizing children giving
b_parent <= b_m+b_n+K(m+n)^(1-delta), uniformly at comparable large orders.
Balanced merge trees would make this defect summable. No such theorem, or
strictly simpler sufficient minimizer property, was proved.

## New scoped geometry results

- A translation-invariant pseudometric preserving the normalized deficit
  with a uniform modulus on the full coset quotient has covering entropy
  (log2/2+o(1))n^2. It cannot be a uniformly compact carrier.
- Every leader support contains an isometric quotient-Hamming cube.
  Deepest leaders give exponentially many near minimizers, but shrinking
  halos have shrinking separation; this alone does not obstruct their
  fixed-resolution compactness.
- Two explicit low-rate code families have identical primal/dual
  enumerators and separated N^(3/4) deficits. Their dimension is
  Theta(N^(3/4)), so they do not falsify a special cut-code theorem.
- All translated queries are stronger than physical disjoint-continuation
  queries. The archived physical compiler has quadratic calibration.
  Its information lower bounds are not low-cap parent constructions.

Proofs and precise assumptions:
[critical-scale code audit](artifacts/critical_scale_code_audit_2026_09_05.md).

## Strongest older constraints

1. Bounded moment, local-profile, restriction, and fixed-level SOS statements
   miss leading extrema on their proved classes. Do not generalize beyond
   each counterexample's hypotheses.
2. Action compactness lacks all-order lossless sign recovery; projective
   exchangeability and uniform sampling have specific iid obstructions.
3. Sign-near weighted recovery has a rounding theorem but its existence
   is not a strict reduction and retains almost all edge bits.
4. Finite-temperature, posterior, transport, and sparse-repair routes have
   no current summable recurrence. New quantities must remove an obligation.
5. Linear-degree radial moments approximate the maximum but retain
   high-order signed cancellation. Raw Bonferroni needs exponential rank.
6. Bounded local stationarity permits suboptimal signings.
7. Arithmetic/conference examples do not provide separated universal order
   classes, and hence do not prove nonconvergence.

## Current independent campaign: verified additions and gaps

The latest user authorized six substantive hours starting 18:29 UTC on
2026-09-05. The six-hour campaign and its final independent audits are complete.
The user has now authorized a second six-hour campaign, from 01:13:55 to
07:13:55 UTC on 2026-09-06; it is active. Its first checkpoint is below.

- Proved polar-Gram inequality:
  Q(A)>=n(n-1) asin(n/||A||_*)/pi. This constrains the nuclear mass of
  low-cap signings but does not improve the universal lower constant.
- The same-spin regularization R(B)=sup_s Q(H_s tensor B)/s^(3/2), for
  explicit regular Hadamard orders s=4^a144^b, has a convergent minimum.
  Its limit c_R lies in [1/sqrt(2pi),1/2] and bounds the original limsup
  from above. Equality with the original liminf is unproved. Full sign
  seeds remove the leading diagonal-completion cost. PSD-majorant
  certificates cannot give a full-seed regularized ratio below 1/2.
- Corrected entropy interpretation: convergence of ordinary
  n^-2 log(1+good-signing count) for dense fixed cap thresholds implies
  convergence of M_n/n^(3/2). Independent edge noise makes every seed
  exponentially numerous after any fixed cap relaxation. A genuine
  lower-tail LDP also suffices; adjacent published Gaussian/spectral
  theorems do not establish this Bernoulli Boolean assertion.

- A universal rooted Gaussian theorem now gives the improved original
  lower bound above. If B=A/sqrt(n-1), S is uniform, and h is fixed smooth
  even, then (B[S h(BS)])_i converges uniformly to N(0,E h(Z)^2) for
  every bounded-normalized-cap sequence. Explicit chaos contractions and
  endpoint-spin transport avoid assuming generic AMP universality.

An elementary injective-tree moment theorem and own-spin-free transport now
extend the energy identity to every fixed finite odd-degree tree hierarchy.
Its Gaussian isometry has stable scalar and anchored fixed points. An exact
two-field mask improvement gives the displayed lower bound, with a rational
interval certificate enclosing all one-dimensional integration errors. Every
infinite-space operation is realized by finite approximations before n grows.
Finite ancestor anchors now strengthen the fixed point to the displayed
0.4306581794055286 bound; the 21-anchor certificate is exact and audited.
The nonnegative-mask certificate has a proved ceiling below 0.4495; this does not
bound all actual algorithm outputs. A new joint unmarked channel escapes that
certificate class. For Q=B^2 and R=Q^(circ3), Schur Jensen proves
sum_i sqrt((BRB)_ii)>=n. A weighted cubic projection theorem and exact paired
cube means give a gain c/log L over the banked Gaussian certificate for
every sufficiently large fixed operator cap L. Principal deletion transfers
this to a strict universal improvement. The increment is not evaluated.
[Complete theorem](artifacts/fresh_normalized_unmarked_gain_and_limit_order_2026_09_05.md).
Classical diagonal factorization plus random refill also shows that restricting
||A||op<=L sqrt(n) changes the normalized minimum by at most O(L^(-1/2));
fixed-L convergence remains unproved.
Unverified numerics are not frontier.
The strict gain now exceeds the supremum of the entire scalar central-mask
hierarchy, uniformly, not only the selected finite-anchor certificate.
All fixed odd channels can be combined before clipping; a finite-degree
nonlinear response bound is proved. The independently audited headline
numerical lower bound can bypass unbounded-operator smooth transport by
working after spectral deletion and using ordinary Gaussian L2 approximation.

Two exact full-sign seeds of order 12 have the same square, but rational
PSD-majorant certification and regular Hadamard lifting give an asymptotic
normalized cap gap at least 437/(2000*12^(3/2)). Hollowing preserves the
gap and makes their normalized squares differ by o(1) in operator norm.
Every fixed odd-channel one-root covariance also agrees asymptotically.
This is not a near-minimizer example, full action-limit equality, or
nonconvergence. [Exact theorem](artifacts/fresh_cosquare12_scalable_gap_2026_09_05.md).
Schmidt's primary theorem settles the formerly unresolved H2 regularized
seed at sqrt(2); greedy stagnation did not obstruct its global optimum.
[Mapping](artifacts/fresh_schmidt_odd_walsh_regularization_2026_09_05.md).
A final monomial argument uniformly selects a nonlinear odd direction
through degree 417 for every marked mask with certificate at least .43.
It needs neither a regularity hypothesis nor a covariance compactness
argument. Uniform local slack beyond the scalar class remains open;
high value does not force even intrinsic output innovation to stay positive.
[Closing synthesis and exact original-problem gap](artifacts/fresh_limit_final_synthesis_2026_09_06.md).
The resumed campaign sharpens the degree to 51, with an exact rational replay.
It also falsifies general local slack: a genuine mask has J>.4301875 and
H=1 throughout |UH|<=.001. A newly proved restricted nonlinear-channel theorem
bypasses this obstruction. Exact feasible means +/-F+H sign(BF), nonlinear
joint tests, and Schur mean-standard-deviation Jensen give
liminf M_n/n^(3/2)>=E H E|U*F+sN| for fixed feasible old pairs and finite
nonlinear edge mass s. Two independent reconstructions check the proof.
Consequently the original liminf strictly exceeds the ENTIRE nonnegative
marked-mask supremum, by a uniform symbolic positive increment. No new
decimal was evaluated at that checkpoint. The subsequent full nonlinear
variance and usable innovation theorem now give the headline decimal above.
The full variational supremum admits ternary F and binary H=1-|F| by
purification after finite approximation. Iteration is not licensed: exact
coordinate-ascent plateaus and an actual twin-spin non-product limit refute
simple refreshing, with explicitly limited scope.
[Checkpoint and proof links](artifacts/resumed_limit_campaign_2026_09_06.md).
Convergence/nonconvergence and a Level-6 recurrence remain open.
The third checkpoint proves a uniform increment of 1/40000 above the
ENTIRE old marked-mask supremum, and a one-coordinate strict ascent for
finite canonical responses (also for the particular cyclic core). The
full gradient admits one-new-Gaussian ascent with an explicit positive
line-search bound; there is no uniform convergence rate. Finite strictly
causal stationary frames have value at most 1/sqrt(2pi), with genuine
zero-gradient ties retained. Cyclic frames are outside that theorem.
Crucially, fixed-rule GFOM universality on flat involutions bounds the
ENTIRE full Gaussian certificate supremum by sqrt(15)/8<1/2. This does
not bound actual signings or growing-complexity/nonlocal algorithms.
The independent nonlocal track realizes arbitrary finite-probability
averaging reflections by exact H144 three-cell and simplex quotients.
Products on separate interfaces are valid; overlapping composition and
general polar recovery are not proved. These are tests, not a matching
upper or an all-order realization theorem.
The fourth checkpoint removes the finite Gaussian feasibility assumption:
normalized even vector features admit their unique canonical creation frame
exactly when the matrix-height kernel iterates to I. A strict subfixed
barrier also gives a quantitative approximate-realization obstruction.
Continuous cyclic cores plus finite causal extensions have a uniform
two-stage escape above 1/sqrt(2pi) on each FIXED frame; no rate uniform over
growing frames follows. Arithmetic phase carriers and vectorial bent maps
give new nonlocal regularized-norm tests, not upper-preserving recovery.
An exhaustive q=16 balanced-profile test found exactly 7/2 for its weighted
seed; that finite result is not an upper bound on the full regularization.
[Proof audit and exact scopes](artifacts/resumed_director_fourth_checkpoint_audit_2026_09_06.md).
Historical route judgments are evidence to audit, not binding directives.
See [the continuing campaign](artifacts/fresh_limit_campaign_2026_09_05.md).

The verifier uses exact finite enumeration for the small code examples,
leader cubes, slice mapping, and reversal identity, and reports analytic
bound evaluations separately from those finite exact facts.
