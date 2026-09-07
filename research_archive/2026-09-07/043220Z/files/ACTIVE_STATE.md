# Active research state

Evidence cutoff: closing work of the six-hour decisive-limit campaign,
2026-09-06 22:44:53 through 2026-09-07 04:45 UTC (still active).
Read the [closing synthesis](artifacts/decisive_campaign_final_synthesis_2026_09_07.md)
and STEERING; use ledger 10.151 and proof dependencies as needed.

## Exact problem

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad Q(A)=\max_x|H_A(x)|,\quad
M_n=\min_{A\text{ hollow symmetric full signing}}Q(A).
```

Determine convergence or nonconvergence of M_n/n^(3/2), without fixing a value.

```math
0.4333221116640807
\le\liminf M_n/n^{3/2}
\le\limsup M_n/n^{3/2}<0.494515125.
```

Neither convergence nor nonconvergence is proved. Convergence to 1/2 is excluded.
The upper is approximately .494515124725174; the old upper was .499432211.

For d=binom(n,2), the AUGMENTED cut code is
C_n^+={(c+b_i+b_j)_(i<j)}. Then Q(a)=d-2dist(a,C_n^+) and
M_n=d-2rho(C_n^+). The constant coordinate is essential for absolute energy.

Recorded exact M_3,...,M_14=(3,4,4,5,9,10,12,13,17,18,20,21).
The n=11,13 lower solver proofs were not rerun here. External witnesses give
M_15<=27,M_16<=30; their global lower-catalogue completeness was not replayed.
None of the new asymptotic theorems depends on these finite optima.

## Upper theorem: necessary dependencies only

[Standalone proof](artifacts/decisive_audit_standalone_direct_E_upper_2026_09_07.md)
and [independent reconstruction](artifacts/decisive_independent_upper_realization_audit_2026_09_07.md).

1. Arithmetic/harmonic precision Schur pivots and information chain rules:
   the mean-conditional-variance envelope E satisfies BE<=E.
2. Gaussian-boundary direct stopping controls finite-depth Bellman iterates.
   Its tail/stopping budget is Phi-G, not an assumed global Phi-E gap.
3. The exact phase E_4(nu_(31/32))=g_4(1) is certified by an analytic
   reproduction-support reduction and 212505 directed-interval rectangles.
4. Full-spin and both-polarity counting, independent fibre averages (E Z)^m,
   a terminal-orthogonal-uniform Fock bound, and O(N) diagonal correction
   produce actual hollow full sign matrices.
5. H2/H12 terminal orders and restriction give all orders. Fix desired margin,
   choose finite depth, send order to infinity, then remove margin.

The final upper is [4+(31/32)log2+g_4(1)]/[8sqrt(31/32)] and is at most
7787631971809/15748015748016<.494515125. The phase and final number have
multiple independent full replays. H=E is a separately proved certificate
identity, not a premise equating this ensemble to M_n.

## Lower theorem and new consequences

[Fresh lower reconstruction](artifacts/decisive_audit_fresh_full_lower_chain_2026_09_07.md).
The actual marked/two-Gaussian response, tree-frame and covariance estimates,
nonlinear replacements and spectral deletion yield c*=.4333221116640807.
The 256-cell exact policy certificate was replayed byte-for-byte.
Finite approximations precede order; deletion cutoff is removed last.

Writing P=max H_A,R=max(-H_A),W=(P+R)/2, the same asymptotic lower holds
for min_A W(A). The square bipartite full-sign extension yields
liminf min_C ||C||_(infinity->1)/m^(3/2)>=2c*>.8666442233281615.
It retains the complete-bipartite support, not an arbitrary weighted graph.

## Useful new mechanisms and their actual gaps

- Bounded actual cap gives delocalized quadratic Gibbs fluctuations throughout
  a fixed bridge interpolation. Contracted-temperature gains do not supply
  same-temperature almost-subadditivity.
- The actual weighted WIDTH-minimizer log-cosh path has a favorable PSD
  derivative term plus an uncontrolled signed edge-flip-cost imbalance.
  Taylor remainders are sublinear; leading integrated payment remains open.
- Sparse actual sign surgery repairs the midpoint under a proportional
  semidefinite ramp. It works for trace o(sqrt(n)), or sign-compatible trace
  o(n) with explicit entrywise feasibility. No such certificate has been
  obtained for selected width minimizers.
- Soft-flatness approximates the original optimum within log2/(4tau)
  uniformly in dimension. Spectral regularization preserves leading outer
  entropy. The required all-order variational limit is not proved.
- Scalar-affine Gaussian-sign expected ABSOLUTE energy falls below width by
  a fixed gap for actual low-cap sequences, using BOTH current endpoints and
  a cap-only four-sign variance theorem. This is a method restriction,
  not a limitation of actual Boolean signings.

## Strongest scoped falsifiers

1. Child reversal forbids unpaid bridge cancellation.
2. Prescribed dense Hadamard stabilization has a half-floor for any full
   seed; this is not a universal construction obstruction.
3. Generic bounded-cap temperature payment is scalably false, not a theorem
   against actual global-minimizer payment.
4. Sparse variance-one amplitudes disprove two-sided adaptive universality;
   fixed bounded favorable flatification is not refuted.
5. Strong high-layer capture is circular. Steiner-frame full-sign matrices
   force linear trace even for proportional PSD ramps; they are NOT
   near-minimizers.
6. Finite skew-lift failures do not imply an asymptotic Boolean floor.
7. Contracted scalar recurrences allow abstract oscillating slopes; no
   realizability of those counterprofiles is asserted.
8. Existing compactness/ultrafilter limits do not supply all-order recovery.

## Exact next construction target

For comparable N=m+n and selectable actual minimizing children, globally
round the row-square-normalized weighted block diagonal
diag(sqrt((N-1)/(m-1)) A_m,sqrt((N-1)/(n-1)) A_n)
to full signs with extra cap O(N^(3/2-delta)), delta>0.
Then u_N<=u_m+u_n+O(N^(1-delta)), u_n=M_n/sqrt(n-1), forces convergence.

No such rounding or equivalent original-value recurrence is currently proved.
The original liminf-to-all-order gap remains. A nonconvergence proof instead
needs genuinely separated infinite order subsequences, not route falsifiers.

## Research preservation

Follow README: never /tmp; preserve research-bearing ignored/untracked files
in dated tracked archives even when unfinished. Keep canonical proofs in place.
See research_archive/ for reviewed exclusions, dependency recovery and hashes.
