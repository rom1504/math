# Strategic steering

Evidence cutoff: the consolidation after Wave 56, later Wave 57--58
refinements, and the completed bounded diagnostic campaign (2026-07-30).
Ordinary wave generation is paused. The committed recommendation is to seek
external mathematical review before deciding whether any research resumes.
If ordinary waves later resume, Wave 61 remains the next scheduled refresh
and must include a blank-slate abstraction audit.

## User-stated objective and workflow directives

The research objective is to determine whether `M_n/n^(3/2)` converges. The
conjectural value `1/2` is not an additional user objective.

The pause policy, consolidation rule, requirement to investigate convergence
and genuine nonconvergence, and bounded three-agent campaign are user workflow
directives. Route rankings and mathematical conjectures below are the main
agent's judgments. Suggestions from subagents, external reviewers, or other
model instances are feedback to evaluate, not directives.

## Current decision

Continuous autonomous waves remain stopped. The campaign produced exact
diagnostic constraints but no primary progress: no improved bound, completed
recurrence, uniformly proved exact-minimizer bridge, strict reduction,
scalable asymptotic counterexample, or convergence/nonconvergence mechanism.

The evidence supports **seeking external mathematical review**. It does not
support resuming one named internal target, because no target was shown
strictly weaker than the original tail. It does not support changing
architecture, because no new mechanism survived audit. A permanent stop
would overstate the evidence because no impossibility theorem was proved.
Further autonomous research requires an explicit later decision; it must not
restart automatically.

## Rigorous frontier

The interval remains

~~~math
0.336493364431\ldots\le\liminf\frac{M_n}{n^{3/2}}
\le\limsup\frac{M_n}{n^{3/2}}\le\frac12.
~~~

No fixed-density restriction recurrence is complete. The verified final
chain is still

~~~text
uniform recurrence with O(n^(3/2-c)) error
  -> summable geometric landing
  -> convergence of q_n/n^(3/2)
  -> convergence of M_n/n^(3/2).
~~~

## Selected-prior and common-active-face verdict

The selected-prior package is the cleanest exact characterization, but the
diagnostic audit proves it is not a reduction. If a prior has favourable mass
`Z` and favourable-captured row `Rbar`, then for every `a>1` one supported cut
has mass at least `(1-1/a)Z` and row at most `a Rbar`. Conversely a good cut
gives the point-mass prior. At the project scales, this is exactly the bare
favourable low-row tail up to fixed constants.

Captured parent deficit is optional for the bare-tail theorem because deficit
already appears inside `widehat ell`. Requiring the event inside the
particular common-active-face law and canonical escape adds compatibility; it
does not weaken the tail.

The common-active-face theorem is correct and noncircular, but it controls
polynomial-scale first moments. Its certificate is insensitive to arbitrary
`e^{-Theta(H)}` perturbations and controls neither saved bare favourability
nor row conditioned on that rare event. Scalar averaging is circular because
`E widehat ell` is exactly the unknown restriction excess. This route should
not resume without a new entropy-scale theorem stated in a quantity not
already equivalent to favourable mass plus row.

## Blank-slate diagnostic audit

The diagnostic generated three formulations before comparing them with the
ledger.

1. **Block `2/3`-power composition.** A power-saving, Hammersley-summable
   defect would prove convergence. Bare `o(N)` error permits slow log-log
   oscillation, and known bridges have a leading-order defect.
2. **Far-down principal restriction.** Competitive signings are
   cut-quasirandom and contain every fixed signed pattern. For fixed child
   order and parent order tending to infinity, the best restriction is
   eventually exactly the order-child optimum, so the proposal becomes the
   scalar conclusion itself. Coupled scales return to the known
   hidden-versus-revealed adaptive-selector gap.
3. **Fixed-temperature pressure.** A same-temperature power-saving
   near-additivity theorem would suffice. The exact random bridge has
   extensive curvature, changes child temperatures, and switching makes
   ordinary covariance interpolation independent of the signing.

These are not new architectures with leverage. Fixed-pattern compactness,
bare little-oh composition, and covariance-only interpolation should not be
reintroduced as independent routes.

## Genuine nonconvergence alternative

Genuine nonconvergence remains a standing alternative. It requires a fixed
`epsilon>0` and two infinite subsequences separated by at least `epsilon`, or
an equivalent proof that liminf is strictly smaller than limsup. Failure of
an intermediate lemma does not count.

The campaign found no such construction. It proved an explicit
`O(sqrt(beta))` relative-order modulus for the actual normalized optima, where
`beta` is relative order separation, and hence

~~~math
\operatorname{Clust}\left(M_n/n^{3/2}\right)
=[\liminf M_n/n^{3/2},\limsup M_n/n^{3/2}].
~~~

If nonconvergence is real, low and high values must occupy alternating
macroscopic intervals of positive width on the `log n` axis. Ratio-dense
arithmetic classes cannot form distinct phases. A surviving proposal must
start with both a strict low construction and an order-sensitive larger lower
bound for every signing on another sparse hierarchy. No such universal lower
mechanism is known.

## External-review questions

The recommended review should ask:

1. Is there a known proportional principal-submatrix or discrepancy theorem
   controlling the adaptive selector gap with a power saving for exact
   minimizers?
2. Can exact-minimizer structure yield entropy-scale favourable low-row mass
   without assuming an equivalent bare-tail statement?
3. Does covering-radius, Boolean quadratic-form, spin-glass, or design theory
   provide a thermodynamic-limit mechanism for the minimum over signings, or
   an all-signings lower theorem capable of genuine log-scale oscillation?

The complete implication chains, falsifiers, proofs, and campaign synthesis
are in `ledger.md`; Git history preserves earlier rankings. Any external
proposal must be checked against those obstructions before it can justify a
new target or architecture.
