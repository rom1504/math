# Benchmark validation campaign: director scorecard

**Date:** 2026-08-16
**Scope:** commits `660ad18` through `e11e281`, plus the final audited
clarifications in this checkpoint.  All benchmark experiments were rerun
together; every advertised finite check passed.  Proofs, not enumeration,
carry the theorem claims.

## Scorecard

| problem | state predicted by contextual response | classical state | independent discovery? | exact / approximate minimality | verdict |
|---|---|---|---|---|---|
| bounded-separator Max-Cut / binary CSP | conditional optimum profile on boundary assignments; for pure Max-Cut, a projective table on `{±1}^w/{±1}` | treewidth separator table, modulo global flip | qualified yes: derived solution-hidden, although the table is nearly forced; lookup and resource theorems came later | every projective coordinate is exposed and every table is realizable by an unrestricted private compiler; under unit load the response class is the full projective Lipschitz ball, with exponential macroscopic response bits; an `m`-edge grammar has only `O_epsilon(m^2+m log(w+m))` bits | strong A/B success; compiler size and the entropy-exponent gap remain |
| 1D / fixed-width Ising | width one: baseline and cavity gap; two-ended fragment: projective `2 by 2` max-plus kernel; fixed width: boundary profile | zero-temperature transfer matrix / cavity field | yes, but low-surprise | endpoint pins prove exact coarseness and max-plus multiplication proves reuse; width-one covers are polynomial in precision; full profile-cube lower bounds apply to gadget-complete treewidth fragments, not yet strict strips | clean sanity-test A success; strict-strip approximate lower remains open |
| heterogeneous / fixed-rank mean field | discrete-concave top-occupancy profile, equivalently sorted slopes; finite-grid histogram; under fixed pair score, the concave response roof | empirical field distribution, magnetization counts, Legendre roof | not cleanly blind; fixed-rank closure was known, but the heterogeneous roof, collapse threshold, and rate law were derived here | exact biconjugacy and metric identity; histogram addition; `eta N/2` depth-independent error; finite atoms have `Theta(n^r_Z)` exact states; equally spaced atoms have exact binomial cover below half-grid error; `J>=4B/n` sharply collapses mass `n` to total field | strong A/C success; arbitrary-continuous-field macroscopic rate remains gapped |
| weighted automata / weighted languages | restricted future residual; under tropical lumpability, nonlinear block maxima | weighted Nerode/Hankel residual and max-plus lumping | qualified yes: solution-hidden derivation, then exact classical collision | residual is coarsest exact state; quotient updates at every depth; robust reachable pins give `Theta(r log(1+B/epsilon))` bits on an exposed aggregate box | successful harder A/B validation; sharp rate needs reachability |
| directed response under interaction | metric-shell strength plus isometry holonomy; anisotropic transported loads; on selector cells, twisted gauges and tangent-reset data | tropical distance transforms, cocycles, finite semigroup actions | agent-authored live target, not solution-hidden | shells close exactly under min-plus interaction; gauges telescope; bounded normal forms absorb coherent defects; on a fixed factorial tie-free selector language, syndetic tangent resets are necessary and sufficient against fresh residuals; a finite suffix-product lift decides the regular case | strong specialized D/E success; generic switching weak coupling remains open |

## New general theorems

### Arithmetic feature-algebra growth

For a fixed finite set of additive response atoms, the exact mass-`n`
contextual quotient has

```math
Theta(n^(r_Z))
```

states, where `r_Z` is the integer rank of the atom-response differences.
Real conditioning controls coarse response covers; lattice margin controls
the resolution below which all exact states remain distinguishable.  This
predicts the polynomial mean-field branch and explains why exact arithmetic
rank can differ from robust exposed dimension.

### Response image and update congruence are separate resources

Landmark covers and balanced exposure price a realizable response image.
Shared-parameter presentation can be much smaller than semantic carrier
rank.  None of those static bounds makes a state reusable: the quotient must
also descend under every continuation.  Histogram addition, max-plus
lumpability, metric shells, gauges, resets, and finite semigroup relations
are distinct verified congruence mechanisms.

The finite-atom theorem also has a continuous upper extension: one common
root-scale net of atom responses gives an additive type histogram with
`n eta` error and binomial state count.  This yields simultaneous sublinear
bits and error for both heterogeneous fields and fixed-rank signed-sum
zonotope support queries.  Adaptive child nets are not allowed: a parent
cannot recover distinctions already collapsed below.
Coordinate atoms attain the full binomial cover, so `D(eta_n)=o(n)` is the
sharp universal atom-net threshold for vanishing summary bits; better rates
must use algebraic dependence.

### Finite robust selector decision

For a regular tie-free selector language, store all products of nonempty
suffixes in a finite lifted graph.  Acyclic height gives a uniform error
bound.  A reachable cycle pumps reset-free words and, under arbitrary
residuals, produces linear drift.  This turns one important depth-stability
hypothesis into a finite verifier/falsifier.

## Revisions forced by counterexamples

1. Unit interface sensitivity is not compression: unit-load Max-Cut still
   realizes the full Lipschitz ball.
2. Low query-parameter dimension is not exposure dimension: one affine
   suffix line can expose every automaton coordinate.
3. A small static cover is not a reusable automaton state.
4. One-step approximate idempotence can accumulate a transition toll.
5. “Kernel gauge or small full-image reset” is false for coherent maps:
   nearby idempotent clamps stay close by a shared algebraic relation.
6. The clamp does not refute every broad zero-holonomy recurrence; coherent
   and fresh-adversarial perturbations require different converses.
7. Exact arithmetic rank is not robust information when real conditioning
   or lattice margin vanishes.
8. A response roof depends on the declared future algebra: changing
   curvature normalization or using non-biaffine interactions can resurrect
   discarded fibres.
9. Literal and projective response norms differ by a factor of up to two;
   the mean-field constants are recorded separately.
10. Static nets at successive scales do not form an update algebra: a coarse
    child atom quantizer can collapse type counts that a finer parent needs.
    The continuous-atom theorem therefore fixes one root-scale net.

## Director judgment

The framework is **Level 2 globally and Level 3 on several finite-feature or
bounded-interface classes**.  It is coherent around one operational law:

```math
usable future compression
=small realizable response image
+future-semigroup congruence.
```

The campaign did not merely rename transfer matrices: it predicted the
Max-Cut Lipschitz obstruction, a strict mean-field roof quotient and sharp
collapse threshold, the arithmetic growth exponent, and the coherent versus
adversarial stability split.  It is not yet a universal representation
theory.  The correct decision is to deepen it, not branch again and not
reconnect to the signing problem.

## Single strongest next theorem

**Paired-selector skew-product theorem with switching cells.**  Given two
locally close finite families of piecewise-affine projective max-plus maps,
build the reachable automaton of paired active cells, including a finite tie
refinement.  After factoring exact common semigroup relations, return one of:

1. twisted coboundaries on every pumpable recurrent selector cycle plus a
   bounded gap to tangent resets, implying uniform `C epsilon` error; or
2. a pumpable paired cycle with unequal transported means, giving an explicit
   `Omega(T epsilon)` drift word.

The fixed-itinerary suffix-product theorem proves the non-switching case.
The missing theorem must allow perturbations to change the active cell
without pretending that coherent fixed kernels are independent noise.

## No reconnection yet

None of the benchmark states is presently both smaller than the full Boolean
response landscape and closed under the signing problem's growing cross-order
interface.  Separator entropy is cautionary, mean-field relies on a finite
anonymous atom algebra, automaton compression assumes exact lumpability, and
metric shells are highly structured.  Development should therefore remain
inside `extremal_information/`.
