# Strategic steering

Evidence cutoff: first focused AR checkpoint, ledger Section 10.143
(2026-08-16), based on commit `c4bc2e0`.
Status: **focused all-order action-recovery campaign, second checkpoint active**.

## User-stated objective and workflow directives

Determine whether `M_n/n^(3/2)` converges.  Convergence to any constant is
success; `1/2` is conjectural but is not the objective.  A rigorous proof of
nonconvergence is also success.

The user authorized a retrieval-grounded panel with ledger-blind specialists,
archive verification, cross-domain criticism, contrarians, and severe
selection.  That panel is complete.  On 2026-08-16 the user explicitly
authorized this focused restart on all-order action realization/recovery
(`AR`) as a standalone problem.  The requested campaign includes the weakest
sufficient theorem, literature-first specialists, adversarial verification,
five realization architectures, Gamma recovery, a falsification track, and
near-minimizer rigidity.  It must consolidate after two substantive
checkpoints without a theorem, strict reduction, or scalable obstruction.

The README verification, Git, stopping, and blank-slate-audit rules remain in
force.  These paragraphs record user objectives and workflow directives only;
all route choices and mathematical judgments below are agent-authored.

## Agent-authored rigorous frontier

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.
```

The interval is unchanged.  The first AR checkpoint did prove a strict
relaxation of exact-sign recovery and two scalable obstructions.

## Leading route: sign-near weighted recovery

For a symmetric hollow weighted matrix `W in [-1,1]^(m x m)`, define

```math
V(W)=\sum_{i<j}(1-w_{ij}^2).
```

The selected target is the following tolerance-by-tolerance statement.

> **WAR.**  For each selected bounded-operator liminf action cluster `T`,
> there is an upward ratio-dense set of orders and weighted matrices `W_m`
> with `V(W_m)=o(m^2)` whose directed one-profiles approach that of `T`.

A weaker scalar version may replace profile convergence by
`Phi(T_Wm)<=Phi(T)+o(1)`.

Why it proves convergence: biased rounding and scalar Bernstein give an exact
hollow signing `A_m` with

```math
Q(A_m)\le Q(W_m)+C\bigl(\sqrt{mV(W_m)}+m\bigr).
```

For profile recovery, deleting `o(m)` high-variance rows makes every row
variance `o(m)`; Bandeira--van Handel rounding then gives
`||A_m-W_m||_op=o(sqrt(m))`.  Directed action continuity transfers `Phi`, and
principal deletion fills every order below the ratio-dense sequence.  The
entire implication has been independently verified.

Why this is a strict reduction: it removes exact target-order integrality.
Fractional entries are allowed on as many edges as desired provided their
total variance is subquadratic.  It does not assume `M_m` or an exact signing
at the target order.

Why tractability is still unproved: the scalar version still quantifies over
all Boolean spins, and the profile version may still encode the dangerous
microscopic spikes.  Naive blow-ups have `V=Theta(m^2)`.  Therefore WAR is a
rigorous reduction of the realization layer, but is not yet demonstrably
simpler than the complete tail problem.

Exact falsifiers:

1. prove every weighted realization of some selected cluster with the needed
   upper profile has `V(W)>=c m^2` on infinitely many target orders; or
2. prove that specifying its directed one-profile to the required accuracy
   reconstructs the complete Boolean energy landscape.

## Established obstructions at this checkpoint

1. The literally weakest scalar recovery condition is equivalent to
   convergence after purification and deletion, so it is circular.
2. One projectively consistent exchangeable sign array cannot be extremal.
   Tight normalized operator norm forces iid Rademacher edges, whose greedy
   cap is at least `(2/3)sqrt(2/pi)>1/2` almost surely.
3. If `W=E_mu A` for an edge-sign law `mu`, then
   `D(mu||U)>=(binom(m,2)-V(W))/2`.  A sign-near barycentre costs
   `Theta(m^2)` entropy; an `O(m)`-cost tilt cannot construct WAR.
4. Fixed-kernel graphon sampling, published dense-graph Gamma recovery, and
   fixed-template absorption have the wrong quantifiers or scale.  They do
   not supply inverse exact-sign recovery for `A/sqrt(m)`.
5. Paley compression gives all-order operator-flat signs, but only for the
   conference comparison family, not an arbitrary selected liminf profile.

## Ranked alternatives

1. **Order-dependent microcanonical recovery.**  It evades the exchangeable
   obstruction, but must pay quadratic entropy or use a non-barycentric
   mechanism while preserving the extreme profile.
2. **Dense design/absorption.**  Promising only if the necessary directed
   action information admits finitely many constraints at fixed accuracy;
   fixed local graph statistics alone are too coarse.
3. **Direct `AR_min^->`.**  Exact signs on upward ratio-dense orders with only
   one-sided one-profile control suffice, but no construction is known.
4. **Near-order insertion.**  Adding `h=o(n)` random signed vertices costs
   `o(n^(3/2))`; useful only after obtaining a ratio-dense good sequence.
5. **Genuine nonconvergence.**  Logically open.  It requires a fixed positive
   recovery/objective gap on infinitely many orders, not vanishing parity or
   design residue effects.

## Checkpoint decision

Checkpoint one counts as theorem-level progress: it established WAR and two
scalable no-go results, so the focused campaign proceeds to checkpoint two.
The live question is whether nonprojective microcanonical, design, or Gamma
machinery constructs WAR using less than the full Boolean landscape.

At checkpoint two, archive-check all independent proposals and select at most
one theorem for proof/disproof.  If none survives as a strict, testable
reduction, consolidate and stop this campaign.  Do not restart broad search.
If ordinary waves ever resume, Wave 61 remains the next scheduled blank-slate
boundary.
