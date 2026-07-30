# Strategic steering

Evidence cutoff: the first substantial computational--composition checkpoint
(2026-07-30), including certified values through order 12. Status:
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

The leading route is now a scalable algebraic/entropy-aware bridge theorem,
prompted by a certified exact `6+6` composition at order 12.  The exact
order-12 minimizer splits into order-6 conference minimizers `S,T` and a sign
bridge `C` satisfying

~~~math
S^2=T^2=5I,\quad SC+CT=0,\quad
CC^{\mathsf T}=6I+2S,\quad C^{\mathsf T}C=6I-2T.
~~~

This is genuine cross-order structure, but it is not yet an asymptotic
mechanism.  The exact sufficient result remains: for
`b_n=M_n^(2/3)`, construct comparable-order parents with
`b_(m+n)<=b_m+b_n+e(m+n)` where
`sum_j e(2^j k)/(2^j k)` tends to zero as `k` tends to infinity.  Balanced
merging then forces `b_n/n`, and hence `M_n/n^(3/2)`, to converge.

The concrete target is to determine whether the conference/bridge identities
have a scalable closure whose **Boolean cap**, not merely spectral norm,
meets that defect condition.  A complementary target is an energy-level-count
lemma strong enough to choose a state-dependent bridge.  The route is
falsified as an asymptotic strategy if the identities are sporadic, if every
extension has linear `b`-defect, or if they control only spectrum while the
Boolean cap retains a leading-order loss.  Fixed `5+6` children already show
that the ideal zero-defect law is false: their exact bridge cap is 17 and the
`2/3`-power defect is `1.167629...`.  A rigorous random-bridge union bound is
also much too weak at these orders.

Ranked alternatives are: broader exact/nested computation to identify the
algebraic family; refined energy-histogram bridge selection; a correctly
mapped covering-code theorem; and genuine nonconvergence.  Update this
assessment only when the target changes, a major obstruction is proved, or at
a substantial checkpoint.

## Rigorous frontier

The interval remains

~~~math
0.336493364431\ldots\le\liminf\frac{M_n}{n^{3/2}}
\le\limsup\frac{M_n}{n^{3/2}}\le\frac12.
~~~

The newly solver-certified finite values are

~~~math
M_{11}=17,\qquad M_{12}=18.
~~~

They do not change the asymptotic interval. No fixed-density restriction
recurrence is complete. The verified final chain is still

~~~text
uniform recurrence with O(n^(3/2-c)) error
  -> summable geometric landing
  -> convergence of q_n/n^(3/2)
  -> convergence of M_n/n^(3/2).
~~~

## Inactive routes: selected prior and common active face

These are not current targets. The selected-prior package is exactly
equivalent, up to constants, to the bare favourable low-row cut: Markov
extraction converts the prior to one cut, and a point mass gives the converse.
Captured parent deficit is already present in `widehat ell`, so it does not
weaken that tail obligation.

The common-active-face theorem is correct but controls polynomial-scale first
moments, not exponentially rare bare favourability or its conditioned row.
Scalar averaging returns exactly the unknown restriction excess. Neither
route should reactivate without an ingredient that avoids both equivalence
and entropy loss.

## Blank-slate diagnostic audit

The blank-slate audit found three formulations. Power-saving `2/3`-power
composition suffices, but bare `o(N)` defect permits slow oscillation.
Far-down fixed-order restriction becomes the scalar conclusion, while coupled
scales return to the adaptive-selector gap. Fixed-temperature pressure needs
power-saving near-additivity, but a random bridge changes temperature and has
extensive curvature. Thus fixed-pattern compactness, bare little-oh
composition, and covariance-only interpolation are not independent routes.

## Genuine nonconvergence alternative

Genuine nonconvergence remains a standing alternative. It requires a fixed
`epsilon>0` and two infinite subsequences separated by at least `epsilon`, or
an equivalent proof that liminf is strictly smaller than limsup. Failure of
an intermediate lemma does not count.

No such construction is known. The diagnostic proved an `O(sqrt(beta))`
relative-order modulus for normalized optima and hence

~~~math
\operatorname{Clust}\left(M_n/n^{3/2}\right)
=[\liminf M_n/n^{3/2},\limsup M_n/n^{3/2}].
~~~

Nonconvergence would therefore need alternating macroscopic log-scale epochs,
not ratio-dense arithmetic classes. It still needs both a strict low
construction and a larger order-sensitive lower bound for every signing on a
second hierarchy; no such universal mechanism is known.

## Connected campaign tracks and checkpoint

Continue exact certification, structural comparison, and bridge composition
as connected tracks. Preserve programs and logs; label exact, certified,
heuristic, and open claims. External coding results count only after checking
the absolute two-sided maximum and all normalizations; one-sided frustration
is not interchangeable with this problem.

At the four-to-six-hour checkpoint, integrate all three tracks in `ledger.md`,
update this assessment, and name the strongest defensible next target. If no
track yields one, commit a concise negative report specifying the missing
resource or idea. The next scheduled ordinary-wave blank-slate boundary would
still be Wave 61, but this campaign is not being counted as ordinary waves.
