# Active research state

Evidence cutoff: ledger Section 10.151.3, third decisive-limit checkpoint.
The additional campaign began 2026-09-06 22:44:53 UTC and remains active,
targeting 2026-09-07 04:45 UTC. The precision-side Schur argument below is new.
The six-active-hour campaign completed its closing audit around 22:15 UTC,
excluding the service interruption around 16:05–17:02. Read the
[final synthesis](artifacts/transfer_campaign_final_synthesis_2026_09_06.md).
Read this compact state and current STEERING; use the ledger only as needed.

## Original problem and rigorous interval

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad Q(A)=\max_x|H_A(x)|,\quad
M_n=\min_{A\text{ hollow symmetric signing}}Q(A).
```

Determine convergence or nonconvergence of M_n/n^(3/2), without fixing its value.

```math
0.4333221116640807\le\liminf M_n/n^{3/2}
\le\limsup M_n/n^{3/2}<0.494515125<1/2.
```

The new upper is `[4+p log2+g_4(1)]/(8sqrt(p))`, with p=31/32;
approximately .494515124725174. The exact Gaussian phase E_4(nu_p)=g_4(1)
is proved by an analytic three-atom reproduction-support theorem and a
directed-interval exclusion of all high-precision channels.
[Phase and all-order consequence](artifacts/decisive_director_gaussian_phase_all_order_upper_2026_09_07.md).
The first checkpoint's looser E bound gave .4968760940775481.
[Independent all-order reconstruction](artifacts/decisive_bridge_improved_all_order_cap_audit_2026_09_06.md)
uses the [precision Schur supersolution](artifacts/decisive_director_precision_schur_supersolution_2026_09_06.md)
and the replayed rational E certificate. This does not identify the optimal
signing ensemble or prove convergence. The older upper was
(1-2delta)^(-2)[1/2-(a+p[h(delta)+g])/(8sqrt(p))],
p=31/32, delta=2^(-24), a=91470529542342299/20460000000000000000.
where g=(18725/32768)[4delta(1-delta)]^2. Its value is about
.4994322102336964. Convergence to 1/2 is impossible;
convergence itself is OPEN. Recorded M_3,...,M_14 are
(3,4,4,5,9,10,12,13,17,18,20,21); n=11,13 lower solvers were not rerun.
Independent external-witness replay adds M_15<=27 and M_16<=30 only;
the author's M_15 lower catalogue completeness was not replayed here.

For N=binom(n,2), augmented cut code
C_n^+={(c+b_i+b_j)_(i<j)} has Q(a)=N-2d(a,C_n^+) and M_n=N-2rho(C_n^+).
The affine constant is essential for the absolute objective.

## Current strict upper proof

**New simplification:** the mean-conditional-variance envelope
`E_t(X)=sup_L[g_t(E Var(X|L))-I(X;L)]` itself supersolves Bellman.
Child precisions lambda1,lambda2 are replaced by arithmetic/harmonic Schur
pivots; log-majorization pays the constant and the information chain rule
pays the dependence. Thus the already proved Gaussian-boundary limit H is
EXACTLY E. This removes the previous temperature-alignment obligation and
strictly improves the old T-based certificate. It is equality of a Bellman
certificate, NOT of actual ensemble pressure or the original optimum.
The finite-dimensional matrix extension has exact independent tensor
additivity and dependence error at most mutual information; no arbitrary
seed-to-all-orders construction follows from that statement alone.

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
9. Entropy-aware extraction: log(Z_++Z_-) >= N[h(delta)+g]+b(1-2delta)^2 Q-o(N).
   Select one signing with partition no greater than its ensemble mean.
   Threshold/projection/resampling proves g; the same certificate and O(m)
   diagonal cost give the upper coefficient above.
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
- OU smoothing plus threshold projection and coordinate resampling gives
  D(discrete law||product)<=22 eta^2 s^2 N/4+log(N+1) at certified fixed
  parameters. This improves clipped pressure credit with O(sqrt N) error.
  A weighted finite-quantizer extension holds without entrywise flatness.
  Raw Gaussian KL/logdet replacements remain invalid; this is not universality.
- Every fixed retention has SOME liminf-realizing near-minimizers with a
  positive-probability leading random-child loss. Exact minimizers and a
  specially chosen seed sequence are not covered by this strengthened claim.
- Every Paley conference has an isotropic law on absolute ground states;
  every appended signed row costs at least sqrt(n), with parity rounding.
  This is not an asymptotic exact-minimizer obstruction.
- Some actual near-minimizers exclude isotropy on every o(n^(4/3))-slack
  shell, via a coherent clique edit costing O(n^(4/3)). But selectable balanced
  near-minimizers EXIST at the SAME order: Q<=M_n+O(n^(5/4)), |P-R|=O(n).
  Iterated switched-clique repair proves this. Some balanced near-minimizers
  still have Q-I>=Omega(n^(4/3)), where I=max_isotropic E|H|; no exact-minimizer
  or macroscopic-gap claim. Scalar balance is not isotropy or an op-norm bound.
- Pointwise powered aggregation fails even with strict-subhalf parent/children
  and one EXACT minimizing child at asymmetric splits. A balanced capped variant
  also fails. Neither defeats an inequality solely for the actual M_n.
- Exact-minimizer ground isotropy already fails at n=4. Strong edge-local
  optimality cannot replace global minimality: hollow tensor Hadamards have
  orientation gap n and stability against every f<sqrt(n)-1 edge change.
- Hadamard-stabilized seed norm R is <= the absolute-PSD majorant T.
  T(A)>=n sqrt(n-1)/2 for every hollow signing. R=T is OPEN. All-win
  parallel repetition, ordinary quantum bias and parity catalysts differ.
- Actual first marked feedback now has a strict fixed-f,fixed-L gain over
  its first response j_n: liminf[Q(A)/(n sqrt(n-1))-j_n]>=gamma(f,L)>0.
  A high-original-degree probe, signed tensor-rank inequality, mixed nuclear
  comparison and JOINT literal-query Stein law retain the full BC return.
  This does not yet cover the rich 200-degree frame or pay O(1/L) deletion.

[Current campaign index](artifacts/transfer_campaign_2026_09_06.md).

## Exact original-problem gap and active work

Latest scoped conclusion: prescribed dense Hadamard stabilization R(B)
is at least sum|Bij|/(2sqrt(n)), without assuming R=T. It therefore has a
half-floor on full sign seeds, including actual exact minimizing seeds.
See `decisive_independent_prescribed_hadamard_half_floor_2026_09_07.md`.
Actual-optimizer finite-temperature interpolation now has exact edge/cut
cavity bounds and O(sqrt(n)) Taylor errors. Its leading block defect is
uncontrolled, and exact finite monotonicity is false even at global minima.
No cross-order recurrence is proved. The current attack uses stronger global
optimality and integrated comparisons, not another Hadamard surrogate.

A sufficient seed-transfer theorem needs ONE near-liminf seed sequence and,
for each fixed seed, all-order/dense-order constructions preserving its
normalized cap within epsilon_seed->0. Take target order first, seed order
second. Current realization preserves flatness/signs, not arbitrary Q.

Active discriminating tests:

1. The original-value lower criterion b_(m+n)>=b_m+b_n-C sqrt(m+n),
   b_n=M_n^(2/3), would imply convergence by a balanced-tree argument for all
   comparable splits. No such recurrence is proved. Exact-parent partition
   selection passes finite tests, but exchange blockers do not control large edits.
2. Test the precise prescribed-Hadamard catalyst equality R=T rather than
   assuming ordinary quantum or all-win product theorems establish it. The
   strict upper makes this a useful bounded seed-family falsifier if proved.
3. Fixed-L first-marked gain is uniform on compact rectangle families, but its
   current guaranteed constant cannot pay principal deletion. Cap-only full-
   parent ENERGY source approximation is now proved, using uniform row L2 and
   an explicit local Gaussian rate. This is not transported L2 or covariance.
   Rich finite-frame dual splitting and trig small-cut lemmas pass; actual
   high-degree signal and the joint open-root return remain separate gaps.

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
