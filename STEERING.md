# Strategic steering

Evidence cutoff: the completed post-consolidation diagnostic campaign and the
start of the computational--composition campaign (2026-07-30). Status:
**computational--composition campaign active**. This is a new sustained
campaign architecture, not a return to the prior short-wave sequence.

## User-stated objective and workflow directives

The research objective is to determine whether `M_n/n^(3/2)` converges. The
conjectural value `1/2` is not an additional user objective.

The user has explicitly authorized research to resume through a sustained
computational--composition campaign. Its three connected tracks are exact
computation and certificates, structural analysis of exact and heuristic
signings, and cross-order composition as the principal mathematical target.
The campaign should run roughly four to six hours before a global assessment,
integrating the tracks regularly rather than reverting to short proof waves.
The existing consolidation and stopping discipline remains in force.

The success criteria and requirement to preserve reproducible programs,
solver evidence, classifications, ledger updates, and Git checkpoints are
user workflow directives. Mathematical hypotheses, route ranking, and the
choice of a precise composition target are agent-authored judgments.

## Agent-authored campaign assessment

The preferred target is to derive the exact summable-defect criterion first,
then use certified finite data to discover what a bridge block must accomplish
in a construction from orders `m` and `n` to order `m+n`. Candidate
`2/3`-power subadditivity is a test case, not an assumed theorem or exponent.
Exact optimization and structural statistics are support tracks whose purpose
is to constrain that bridge lemma or reveal genuine subsequence structure.

Update this agent-authored assessment when the preferred target changes, a
major obstruction is proved, or at each substantial campaign checkpoint; do
not rewrite it after minor experiments.

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

## Inactive routes: selected prior and common active face

These are not current targets. The selected-prior package is a clean exact
characterization, but the diagnostic audit proves it is not a reduction. If a prior has favourable mass
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
`E widehat ell` is exactly the unknown restriction excess. Neither route
should reactivate unless a genuinely new ingredient avoids these equivalences
and the entropy-scale loss.

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

## Connected campaign tracks and checkpoint

1. **Exact computation.** Reuse the existing programs before writing new
   ones. Seek exact values or rigorous intervals beyond the present range via
   symmetry-aware MIP, SAT/pseudo-Boolean, branch-and-cut, or constraint
   generation, preserving certificates and logs.
2. **Structural analysis.** Compare spectra, rows, ground states, principal
   restrictions, automorphisms, residues, and nearby-order relations for
   exact minimizers and clearly labelled heuristic constructions.
3. **Cross-order composition.** Test bridge-block objectives exhaustively on
   certified finite data, derive the deterministic lemma each objective would
   need, and reject laws whose accumulated defect is not summable.

Coding-theoretic or external results count only after the mapping to the
absolute quadratic maximum, all hypotheses, and every normalization are
verified. One-sided frustration or maximum-cut statements are not silently
interchangeable with this problem.

At the four-to-six-hour checkpoint, integrate all three tracks in `ledger.md`,
update this assessment, and name the strongest defensible next target. If no
track yields one, commit a concise negative report specifying the missing
resource or idea. The next scheduled ordinary-wave blank-slate boundary would
still be Wave 61, but this campaign is not being counted as ordinary waves.
