# Active research state

Evidence cutoff: ledger Section 10.148.6, resumed six-hour campaign active.
Use the linked proofs for reconstruction and the ledger only for archive detail.

## Exact problem and verified frontier

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,\quad
M_n=\min_A Q(A),
```

where A is hollow symmetric with off-diagonal signs. The objective is
convergence or nonconvergence of M_n/n^(3/2), not specifically the value 1/2.

```math
0.4333221116640807\le\liminf_n M_n/n^{3/2}
\le\limsup_n M_n/n^{3/2}\le1/2.
```

Recorded M_3,...,M_14 are (3,4,4,5,9,10,12,13,17,18,20,21).
The n=11,13 lower bounds are solver-certified infeasibility results;
this campaign did not rerun those solvers.

For N=binom(n,2), the augmented cut code is
C_n^+={(c+b_i+b_j)_(i<j)} in binary coordinates. Exactly
Q(a)=N-2d(a,C_n^+) and M_n=N-2rho(C_n^+). For n>=3 its dimension is n.
It is RM(1,n) punctured to the weight-two slice, including the affine constant.

## Strongest current quantitative result

The new original lower endpoint is analytic, with exact rational interval
evaluation and two independent source reconstructions/replays. The finite
Gaussian approximation is chosen before the matrix order tends to infinity;
the spectral-deletion cutoff is removed last.

The full nonlinear response theorem gives, for a fixed feasible odd F/even H,

```math
\liminf_n M_n/n^{3/2}\ge
\mathbb E\,H\,\mathbb E_N|U^*P_1F+\|(I-P_1)F\|_2N|.
```

U is the actual canonical creation isometry, not an arbitrary distributional
fixed point. Nuclear covariance, endpoint parity, and mixed contractions are
proof dependencies. A 21-anchor finite-resolvent cyclic core is realizable.
A second genuinely correlated Gaussian gives an actual two-coordinate policy.
Its old inverse g has norm one; neither its scalar conditional mean nor an
independence approximation replaces it in the norm or masked expectation.

The exact 256-rectangle pure policy gives the lower endpoint
.433322111664080753415812928897579346558634648033693413106996.
The entire omitted conditional Hermite tail is charged. The numerical
optimizer supplies only rational policy endpoints, not trusted moments.

- [Full nonlinear proof audit](artifacts/resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md)
- [Correlated actual birth](artifacts/resumed_response_rich_core_response_birth_2026_09_06.md)
- [Independent optimized-policy audit](artifacts/resumed_bound_audit_rich_core_cell_ascent_2026_09_06.md)
- [Exact certificate](computations/results/resumed_response_rich_core_optimized_policy_certificate_2026_09_06.json)

## Response theory: what is now known and what is not

1. Every old Gaussian inverse has a conditional new-coordinate tail bounded
   by theta^(K+1)||g||_2 when the old/new first-chaos principal angle is at
   most theta<1. Weighted Jensen transfers it without inverse-mask-mass loss.
2. Fixed-frame response optima exist and are uniformly approximated by finite
   cell optima. Numerical KKT stationarity is not global optimality.
3. Analytic cyclic cores with finite causal appendages have strict unrestricted
   high-value escape; the increment is not uniform over expanding frames.
4. The ENTIRE terminal one-response functional is below .45, proved by the
   Gaussian moment body, weighted-square concavity, and an exact 109-interval
   envelope. This is NOT a ceiling for actual signings or all algorithm output
   energies. Their common quadratic endpoint energy was discarded.
5. The older sqrt(15)/8 fixed-GFOM ceiling has broader algorithmic but different
   scope. Neither ceiling covers all nonlocal or growing-complexity methods.

[Current proof map and scopes](artifacts/resumed_director_sixth_checkpoint_audit_2026_09_06.md).

## New nonlocal and actual-minimizer theorems

**Ramsey Fourier saturation.** Every centered even Fourier involution over a
growing finite abelian group has an asymptotically mean-zero Boolean vector
with absolute Rayleigh quotient tending to one, uniformly over multipliers.
Deuber sets and projective vector-space Ramsey supply a finite monochromatic
squarefree frequency pattern; uniform cosine/Hermite approximation supplies
the Boolean vector. Exact multiplicative Fourier structure is essential.

The robust finite-field version proves: every additive-Cayley hollow signing
with ||A||op<=(1+o(1))sqrt(n) has Q(A)/n^(3/2)->1/2. Native Paley cores and
bordered conferences saturate at ALL admissible prime-power orders, with both
Rayleigh signs. No theorem for arbitrary eigenbases, non-Cayley matrices, or
Cayley families with a fixed spectral excess is claimed.

**Gap-two geometry.** A single-coefficient locally minimizing signing with
cap Q has a near-state outside projective radius r of every oriented ground
whenever n(2r-1)<Q and Q>r(n-r)+1. Thus exact minimizers have gap-two states
separated by at least (c/2-o(1))sqrt(n) from every ground, for any universal
lower coefficient c. No operator bound is needed.

If all gap-two states have the same orientation, their number is at least
1+log_2(n^2/(2Q+n)), and their collective changed support is linear. Mixed
orientation instead gives a twisted code-overlap bound; it does not imply
the same collision count. A separate abundant-pattern determinant argument
gives an unconditional (1-o(1))log(n)/loglog(n) count in every C sqrt(n) window.
These are exact-minimizer facts, not claims for arbitrary asymptotic near-minima.

## Retained cross-order theorem and precise gap

For mu_r=E|sum_(j=1)^r epsilon_j|, every n,r satisfies

```math
M_{n+r}\le M_n+M_r+
\min\{n\mu_r+\sqrt{2nr(r-1)\log2},
r\mu_n+\sqrt{2nr(n-1)\log2}\}.
```

For r->infinity, r=o(n), its remainder is
(sqrt(2/pi)+o(1))n sqrt(r). It improves a coefficient, not the exponent or
the previously known qualitative near-order continuity.

The exact bridge identity is

```math
Q\!\begin{pmatrix}A&B\\B^{\mathsf T}&D\end{pmatrix}
=\max_{x,y}(|H_A(x)+H_D(y)|+|x^{\mathsf T}By|).
```

A sufficient open target is a sign bridge for actual minimizing children with
b_parent<=b_m+b_n+K(m+n)^(1-delta), b_n=M_n^(2/3), uniformly at comparable
large orders. No such recurrence or demonstrably simpler sufficient child
property has been proved. Near-state counting/spread alone does not bound
joint row discrepancy. Convergence and nonconvergence remain open.

The Hadamard-regularized seed minimum has a limit, but equality with the
original liminf is missing. Arbitrary controlled Fourier-index translations
are realizable; a single such gate is exactly blockwise phase-only and adds
no seed-norm power. Arbitrary overlapping composition is not established.

## Strongest scoped older constraints

- Translation-invariant coset metrics uniformly preserving normalized deficit
  require exp((log2/2+o(1))n^2) covering states. This is not a near-minimizer-only
  or physical-continuation impossibility.
- Deep leader cubes give many near-minima at vanishing separation, not failure
  of fixed-resolution compactness.
- Spectral squares and fixed one-root covariance do not determine Boolean
  cap: exact cosquare seeds have a scalable cap gap. This is not full action
  equality and not a near-minimizer example.
- Action compactness lacks all-order lossless sign recovery. Sign-near weighted
  recovery retains nearly all edge bits.
- Local-profile, fixed-gadget, Gibbs, scalar-channel and algebraic-family
  obstructions retain their exact hypotheses; do not broaden them by analogy.

## Current campaign and next discriminating work

Authorized six-hour interval: 2026-09-06 01:13:55–07:13:55 UTC. Work continues.
Audit the common endpoint energy omitted by the <.45 terminal certificate,
and the reach of the new nonlocal Ramsey mechanism. Preserve improvements
without equating a certificate ceiling with a signing ceiling. Every
decisive claim needs independent reconstruction and explicit limit order.
