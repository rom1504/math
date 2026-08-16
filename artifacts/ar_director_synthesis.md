# Research-director synthesis: all-order action recovery

Date: 2026-08-16.

Status: second focused AR checkpoint.  The rigorous interval remains

```math
0.336493364431\ldots
\le\liminf\frac{M_n}{n^{3/2}}
\le\limsup\frac{M_n}{n^{3/2}}
\le\frac12.
```

This synthesis ranks mathematical architectures, not the amount of work put
into them.  An `A module` is a proved reusable theorem; it is not necessarily
an `A architecture` for convergence.

## 1. Ranked architecture table

| Rank | Architecture | Class | Exact missing lemma | Director judgment |
|---:|---|:---:|---|---|
| 1 | Extremal-envelope recovery with energy-tail control (`EER_UI`) | B | realize the union of all purified liminf one-profile phases on order sets satisfying the tolerance-dependent covering condition, with uniform integrability of `xy` | strict coherence reduction from selected-phase AR; no constructor known |
| 2 | Near-order transfer | A module / B support | produce good orders with `gamma_eta^(3/2)(alpha+eta)->alpha` | deletion and `o(n)` insertion are solved; density of good orders is not |
| 3 | Outer-profile design plus sparse repair | A repair module / C route | a uniform outer-design theorem controlling every fixed-alphabet coloring before an `O(n)`-edge leave | repair is solved; enforcement is still EER's universal quantifier |
| 4 | Sign-near weighted recovery | A rounding module / C route | construct the weighted recovery input by an independent low-complexity rule | existentially equivalent to exact recovery; almost all edge signs are already exposed |
| 5 | Nonprojective microcanonical/profile pressure | C | an all-order pressure or shell-entropy no-gap theorem at fluctuation scale | stronger entropy regularity than recovery; no interpolation or LDP supplies it |
| 6 | Exposed-face one-vertex absorption | C | sharp discrepancy of the entire near-cap slack landscape at every extension step | exact collision with ledger Section 10.44, equation (10.254) |
| 7 | Uniform mesoscopic induced sampling | D | proposed sampling invariance | false: every bounded-op parent has iid fixed-size induced limits, with a `0.063846...` gap in `Phi` |
| 8 | Projectively exchangeable recovery | D | one nested exchangeable all-order law | false: tight operator scale forces iid signs, whose greedy cap exceeds `1/2` |
| 9 | Ordinary blow-up/conference-fibre residual | D/C | a nonlinear joint absorber not reducible to paid residual channels | ordinary versions are killed by exact Frobenius/ANOVA residual mass; an unspecified absorber renames AR |

There is no class-A convergence architecture and no published theorem that
provides all-order recovery.  The strongest live mathematical statement is
the class-B envelope quotient, not a construction mechanism.

## 2. Candidate card: extremal-envelope recovery

- **Domain:** action convergence, Gamma recovery, compact variational limits.
- **Imported theorems:** Backhausz--Szegedy action compactness; the project's
  verified quantitative continuity of `Phi` under a common `2 -> 2` bound.
- **Problem translation:** for a tolerance `eta`, take the compact cluster set
  `K_eta` of one purified liminf sequence and the closed envelope
  `E_eta=cl union_(T in K_eta) S_1(T)`.
- **Proposed mechanism:** recover only the envelope, not one action phase.
  Each profile of a target-order matrix may match a different cluster phase.
- **Exact missing lemma:** on order sets with upward covering ratios
  `gamma_eta`, construct exact signings with directed one-profile distance to
  `E_eta` tending to zero, uniformly integrable energy products `xy`, and
  `gamma_eta^(3/2)(alpha+eta)->alpha` along a null tolerance sequence.
- **Why it proves convergence:** weak profile matching plus uniform
  integrability passes every quadratic energy integral to the envelope;
  principal deletion multiplies the bound only by `gamma_eta^(3/2)`.
- **Why it is weaker than archived AR:** it forgets cluster identity,
  compatibility of different profiles inside one operator, reverse inclusion,
  all joint profiles, and ratio-density at each fixed tolerance.
- **Why it is not yet executable:** it still excludes a bad profile for every
  one of exponentially many bounded colorings.  No finite-model theorem
  supplies that outer assertion.
- **Ledger collisions:** none for the profilewise phase envelope itself; the
  remaining outer-tail obstruction is the archived Boolean-spike problem.
- **Falsification test:** a real obstruction must occupy multiplicatively
  nonnegligible order intervals or violate every admissible covering sequence;
  isolated parity failures do not suffice.
- **Specialist confidence:** medium in the implication, low in realizability.
- **Verifier confidence:** high in the implication and error budgets.
- **Director judgment:** retain as the correct optimizer-free statement of the
  bottleneck; do not launch proof agents without a new constructor or invariant.

The absolute formal minimum is smaller.  It is enough to match only one
target-order profile attaining `Phi`, with energy-tail control.  Identifying
that profile, however, solves the full target Boolean maximum and violates the
project's information criterion.  Dropping even that law leaves scalar
recovery, which is equivalent to convergence.  Thus `EER_UI` is the weakest
optimizer-free profile theorem found, not the weakest existential sentence.

## 3. Candidate card: near-order transfer

- **Domain:** probabilistic combinatorics, restriction/extension.
- **Imported theorem:** none needed beyond Hoeffding/Bernstein and principal
  restriction.
- **Proved mechanism:** `Q(A[S])<=Q(A)` for every principal submatrix.  Adding
  `h` random signed vertices to order `n` has a deterministic supported outcome
  with error `O(sqrt((nh+h^2)(n+h)))`; hence `h=o(n)` costs `o(n^(3/2))`.
- **Exact remaining lemma:** obtain recovery orders whose upward covering
  ratios obey the tolerance-dependent condition in the first card.
- **Why it is weaker than all-order recovery:** a recovery family may omit
  orders, and its multiplicative gaps need only approach one as the objective
  tolerance vanishes.
- **Falsification test:** exhibit a selected extremal envelope for which every
  favorable order set has a covering ratio bounded above one by a fixed amount.
- **Director judgment:** the order-mismatch subproblem is solved; no theorem
  makes an arbitrary liminf sequence sufficiently dense.

## 4. Candidate card: design/absorption

- **Domain:** dense designs, approximate decomposition, absorption.
- **Imported theorems:** Keevash/Glock--Kuehn--Lo--Osthus-type fixed-template
  decompositions and typical-host packings with `O(n)`-edge leaves.
- **Proved module:** any `O(n)`-edge completion changes `Q` by `O(n)` and has
  normalized `L^infinity -> L^1` error `O(n^(-1/2))`; bounded leave degree is
  unnecessary.
- **Exact missing lemma:** a packing whose signed core satisfies all directed
  outer-profile inequalities simultaneously for every fixed-alphabet vertex
  coloring.
- **Why the lemma would prove convergence:** sparse repair preserves the
  already-enforced profile, then EER applies.
- **Why it is not weaker yet:** fixed gadget counts and local densities do not
  control the `q^n` colorings.  The missing uniform outer-design assertion is
  the whole envelope-recovery obligation.
- **Falsification test:** color by gadget roles, absorber membership, or fibre
  modes and exhibit a separating profile law.
- **Director judgment:** retain sparse repair as a theorem; reject design
  enforcement until a nonlocal concentration or potential-drop lemma is
  stated independently.

## 5. Rejected-idea map

| Rejected proposal | Precise obstruction or collision |
|---|---|
| Uniform induced sampling | `mesoscopic_induced_sampling_no_go.md`: op boundedness forces dense cut limit zero, fixed samples become iid, and greedy energy gives `Phi>=1.063846...` |
| Projective exchangeability | `exchangeable_recovery_obstruction.md`: Aldous--Hoover mean kernel must vanish; the law is iid and nonextremal |
| Objective sign-near weighted recovery | `sign_near_weighted_recovery.md`, (WR.B1)--(WR.B6): exact recovery gives `V=0`, while rounding gives the converse at the same orders |
| Low-cost tilted barycentre | entropy sharpening: `V=o(n^2)` forces `D(mu||U)>=N log2-o(N)` and `H(mu)=o(N)` |
| Canonical profile pressure | alternating zero-state multiplicities make fixed-temperature pressure oscillate even when a zero state exists at every order; archived `good_signing_entropy_threshold.md` has the same endpoint program |
| Hard profile-shell entropy | existence of the all-order entropy limit is stronger than shell nonemptiness and adds an independent no-gap obligation |
| Exposed-face insertion | exact identity and derivative-scale target are ledger (10.254)--(10.256); thick-cap states defeat ground-face balancing |
| Fixed-template design counts | bounded local statistics miss planted fluctuation-scale resonances; universal coloring enforcement remains open |
| Conference/regular fibres | `bounded_op_signed_realization.md` and `regular_microblock_absorption_audit.md`: forced residual Frobenius mass and nonnegative ANOVA channels survive at leading scale |
| Ordinary graphon/Gamma sampling | the kernel for `A/sqrt(n)` has amplitude `sqrt(n)`; published recovery theorems assume the kernels already converge and recover vertex states, not sign kernels |
| Parity, Witt, spectral multiplicity | all located effects have bounded/vanishing repair or ratio-dense carriers; none yields a fixed-positive action gap |
| Fourth-moment or bounded-cycle rigidity | a one-vertex conference deletion has `Delta_4 -> 0` while retaining a complete Boolean kernel spike; fixed statistics miss planted leading maxima |
| Spectral or row-field delocalization | adjoining one universal vertex preserves normalized near-optimality and operator scale but destroys square-field uniform integrability; flat eigenspaces still admit Boolean resonance |
| Ground-state row-law realization | two order-six signings have the same positive labelled row sums but caps `11` and `7`; certifying globality restores the full Boolean maximum |

## 6. Research decision

The campaign executed the two plausible terminal modules and the canonical
sampling proposal.  The modules were proved; the sampler was rigorously
falsified.  Independent specialists supplied no constructor for `EER_UI` and
no intrinsic fixed-positive AR counterexample.  A final rigidity audit tested
fourth-moment/cycle flatness, spectral and field delocalization, and switched
ground-state row laws.  Weak forms lose rare Boolean extremizers; useful
strong forms are false or reconstruct the directed one-profile.  No property
passed forcing, limit closure, and objective-safe all-order realization.

Further autonomous execution inside this architecture is therefore not
justified after checkpoint two.  The appropriate next input is one of:

1. **External action/local-global review:** ask whether profilewise envelope
   recovery with uniform integrability has an inverse finite-model theorem or
   a non-sofic/non-approximability obstruction not present in the current
   literature packet.
2. **A concrete nonlocal constructor:** resume only when a proposed design,
   microcanonical, or algebraic rule states how it controls every coloring
   without carrying an almost complete sign skeleton.  Then assign one proof
   agent and one falsifier to that exact lemma.

Absent one of these inputs, more pressure limits, weighted restatements,
sampling variants, or finite-template absorbers should not be generated.
