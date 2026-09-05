# Active research state

Evidence cutoff: ledger Section 10.147.3, active 2026-09-05 campaign;
second-response theorem independently audited. Use the ledger for archive
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
0.3857858769087784\le\liminf_n M_n/n^{3/2}
\le\limsup_n M_n/n^{3/2}\le1/2.
```

Recorded M_3,...,M_14 are (3,4,4,5,9,10,12,13,17,18,20,21).
The n=11,13 lower bounds are solver-certified infeasibility results.
The fresh verification script does not rerun those solvers.

The improved lower bound is analytic, with an exact rational interval
certificate for its numerical evaluation, not a finite-order extrapolation.
See [optimized second-response proof](artifacts/fresh_limit_response_variational_2026_09_05.md)
and the two linked independent audits of its underlying energy identity.
The upper bound is unchanged.

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
2026-09-05. The campaign is active. Checkpoints do not terminate it.

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

The second-response energy formula is now independently verified. Collective
indirect-path control uses the actual bilinear cap and maximal spin influence,
not termwise cancellation. Optimizing its smooth first mask proves the
displayed bound. A six-vertex counterexample rules out exact finite-rule
positivity but does not contradict the smoothed asymptotic theorem.
The active target is richer joint response rounding and its quantitative
ceiling; no unverified numerical improvement is part of the frontier.
Convergence/nonconvergence and a Level-6 recurrence remain open.
Historical route judgments are evidence to audit, not binding directives.
See [the continuing campaign](artifacts/fresh_limit_campaign_2026_09_05.md).

The verifier uses exact finite enumeration for the small code examples,
leader cubes, slice mapping, and reversal identity, and reports analytic
bound evaluations separately from those finite exact facts.
