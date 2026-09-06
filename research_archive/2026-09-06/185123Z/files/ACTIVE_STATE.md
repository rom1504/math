# Active research state

Evidence cutoff: ledger Section 10.150.3, renewed campaign checkpoint 3.
The six-active-hour campaign remains ACTIVE, approximately through 22:15 UTC,
excluding the service interruption around 16:05–17:02.
Read this compact state and current STEERING; use the ledger only as needed.

## Original problem and rigorous interval

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad Q(A)=\max_x|H_A(x)|,\quad
M_n=\min_{A\text{ hollow symmetric signing}}Q(A).
```

Determine convergence or nonconvergence of M_n/n^(3/2), without fixing its value.

```math
0.4333221116640807\le\liminf M_n/n^{3/2}
\le\limsup M_n/n^{3/2}<0.499432211<1/2.
```

The new upper is C_delta=(1-2delta)^(-2)[1/2-(a+p h(delta))/(8sqrt(p))],
p=31/32, delta=2^(-24), a=91470529542342299/20460000000000000000.
Its value is approximately .4994322102337004. Convergence to 1/2 is impossible;
convergence itself is OPEN. Recorded M_3,...,M_14 are
(3,4,4,5,9,10,12,13,17,18,20,21); n=11,13 lower solvers were not rerun.

For N=binom(n,2), augmented cut code
C_n^+={(c+b_i+b_j)_(i<j)} has Q(a)=N-2d(a,C_n^+) and M_n=N-2rho(C_n^+).
The affine constant is essential for the absolute objective.

## Current strict upper proof

[Standalone reconstruction](artifacts/transfer_reconstruction_standalone_2026_09_06.md)
independently rebuilds and exactly replays the entire chain:

1. Exact restricted rank-one weave, both objective signs and all 2^N spins.
2. PSD folded Gaussian kernels and graph contraction give the one-row
   permanent criterion; independent fibre averages give (E Z)^m.
3. Gaussian Fock orbital bound has exp(O(sqrt(m))) terminal loss uniformly.
4. Finite-depth exact type counting gives the Bellman operator.
5. Unbounded Gaussian self-transport is strictly source-concave; zero drift
   is exactly Gaussian. No controlled-policy uniform CLT is assumed.
6. Conditional-variance envelope T=sup_L[E g(Var(X|L))-I(X;L)] supersolves
   Bellman. Sequential labels, Schur pivots, log-majorization and convexity
   prove this without equal-temperature alignment.
7. Direct stopping at Phi-T bounds deep iterates by T. The older smaller
   latent envelope and its unproved supersolution are not dependencies.
8. Rational 2501-by-2501 posterior grid plus an outward continuous modulus
   gives offset <=-a at p=31/32,t=4.
9. Gibbs-cluster extraction: Z_++Z_- >= exp[N h(delta)+b(1-2delta)^2 Q].
   Select one signing with partition no greater than its ensemble mean.
   This gives C_delta with the same fixed-depth certificate and O(m) diagonal cost.
10. H2/H12 Kronecker terminal orders are multiplicatively asymptotically dense.
   Principal restriction fills all orders. No prime-gap theorem is needed.

Limits: fixed margin, then fixed depth, then all large orders, then margin.
The current campaign's independent reconstruction and exact checks passed.
This construction is not asserted optimal.
[New extraction proof](artifacts/transfer_reconstruction_gibbs_cluster_extraction_2026_09_06.md).

## Original lower proof retained

The .4333221116640807 bound uses the actual nonlinear marked/two-Gaussian
response and a 256-rectangle rational policy with interval evaluation.
The old Gaussian inverse remains norm one; coherent masked dependence and
omitted Hermite tails are retained. Nuclear covariance, parity and mixed
contractions are dependencies. Fixed approximation precedes n, and the
spectral cutoff is removed last.

- [Nonlinear proof audit](artifacts/resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md)
- [Actual correlated birth](artifacts/resumed_response_rich_core_response_birth_2026_09_06.md)
- [Policy audit](artifacts/resumed_bound_audit_rich_core_cell_ascent_2026_09_06.md)
- [Exact lower certificate](computations/results/resumed_response_rich_core_optimized_policy_certificate_2026_09_06.json)

## New verified structural results and strict scopes

- Conditional-variance factorization works for every orthogonal gate and
  heterogeneous marginals, charging total correlation. Every fixed real
  Hadamard gate has Gaussian-only zero drift. This does not accept an
  arbitrary signing as an orthogonal mixer.
- Full rank-one weaving has exact cap (m^3+m|tr S|)/2 after hollowing.
  Balanced retention p has cap >=sqrt(p)/2 in normalized units. Outer signs
  are gauge-erasable when arbitrary Hadamard column signs are available.
- The UNCORRECTED full-spin row-permanent criterion has floor sqrt(15)/8.
  Typical retained-Hadamard spectra are uniformly Gaussian in W2, and exact
  homogeneous endpoint counting proves the same actual annealed pressure floor.
  Actual expected high-spin count has rate >=p log2+(1/4)log(1-4c²p).
  These are NOT typical-signing lower bounds; entropy-corrected extraction escapes
  the one-extremizer counting implication, not its valid partition lower bound.
- Complete even self-tensor classification: fixed Hadamard seeds tend to cap
  coefficient 1/2; every other full symmetric sign seed tends to infinity.
  A fixed-positive random restriction does not repair the latter divergence.
  Every diagonal completion of actual strict-upper weave seeds loses a fixed
  cap amount already at the single power 76.
- Extremely thin chosen restrictions of any non-rank-one seed tensor encode
  ALL signings; sufficiently thin uniform random restrictions approach iid
  signs and have greedy cap >=.531923-o(1). Intermediate scales remain open.
- Near-optimal twin surgery: o(n) copied pairs cost o(n^(3/2)) in Q but can
  force divergent normalized precision-determinant penalties under EVERY
  sign diagonal/precision/order. The displayed defects are o(n)-removable;
  existence of a good seed sequence remains possible.
- Every bounded-cap parent, including exact minimizers, has uniform sparse
  restrictions of cap >=2/pi-o_P(1) for ALL n->infinity,n/D->0. The proof
  uses local graph moments, fixed-degree Gaussian rounding and spectral cores.
  At n0<=n<=rho D, good c<2/pi selectors have fraction at most
  min(K/n,3 exp(-kappa n^2/D)), with constants depending on parent cap and c.
- Rarity is not impossibility: any prescribed bounded-cap child of size o(D)
  can be planted exactly in an asymptotically minimizing parent at o(D^1.5)
  cost. At high retention 31/32, a particular Hadamard parent family has
  typical child cap <.4995; the SAME family has sparse child cap >=2/pi.
- Fixed-marginal clipped correlations give an extensive pressure credit,
  including at the current tiny flip density. Raw Gaussian logdet and
  unattenuated fixed-variance comparisons have explicit counterexamples.
  A threshold-channel resummation is valid with O(sqrt N) error but its
  attenuation makes it inferior to clipping at the current parameters.
- Hadamard-stabilized seed norm R is <= the absolute-PSD majorant T.
  T(A)>=n sqrt(n-1)/2 for every hollow signing. R=T is OPEN. All-win
  parallel repetition, ordinary quantum bias and parity catalysts differ.

[Current campaign index](artifacts/transfer_campaign_2026_09_06.md).

## Exact original-problem gap and active work

A sufficient seed-transfer theorem needs ONE near-liminf seed sequence and,
for each fixed seed, all-order/dense-order constructions preserving its
normalized cap within epsilon_seed->0. Take target order first, seed order
second. Current realization preserves flatness/signs, not arbitrary Q.

Active discriminating tests:

1. Retain the positive correlated-cluster theorem; test sharper discrete
   entropy control only with an explicit pressure/cap consequence.
2. Seek a genuinely cap-preserving original-problem mechanism. Random
   selector loss is now quantitative, but exceptional selectors remain real.
   Neither a necessary local moment defect nor planting into an already-good
   large parent supplies the missing upward transfer.
3. Secondary actual-feedback transport: retain coherent returns QS,QD and
   actual energy. The general high-value marked innovation is still open.

The exact bridge identity remains
Q(parent)=max_(x,y)(|H_A(x)+H_D(y)|+|x^TBy|).
Mesoscopic completion costs (sqrt(2/pi)+o(1))n sqrt(r) for r=o(n), r->infinity;
it does not provide a comparable-order sublinear recurrence.

## Constraints against overinterpretation

The terminal feedback certificate below .45 and fixed-GFOM sqrt(15)/8 ceiling
are class/algorithm limitations, not universal upper bounds on signings.
Finite-field near-flat Cayley and Paley saturation at 1/2 are scoped class
theorems; arbitrary near-minimizers need not be flat or Cayley.
Local profiles, spectral cosquares, action recovery, Gibbs and code no-gos
retain their exact proved hypotheses. Full history is archived in Git.
Research-bearing scratch files are durable in research_archive/; preserve new
work at every substantive checkpoint under the README's updated rule.

Continue after checkpoints. Reassess failed implementations rather than
ending the campaign. Genuine nonconvergence requires fixed positive separation
of two infinite subsequences. No such separation is established.
