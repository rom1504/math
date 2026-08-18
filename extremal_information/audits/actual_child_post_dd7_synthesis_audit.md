# Post-`dd7e43d` actual-child synthesis audit

Status: **adversarial scope audit**.  This note classifies only the results in
the effective-row-support, deletion/extension, convex-superconcentration,
radial-curvature, entropy/orientation, and extension-escort notes.  It does
not modify their proofs.

## 1. Strongest optimizer-specific theorem

The strongest new quantitative statement which actually uses contracted-
temperature child optimality is the **neutral one-vertex extension
envelope** in EE.2.  If `D` is an exact order-`n` minimizer of the augmented
pressure and `t=beta/sqrt(N)`, then every new row `b` satisfies

```math
z_{D,t,t}^{0}(b)\ge e^{-\delta_n(t)},
\qquad
0\le\delta_n(t)\le n\log\cosh t=O_\beta(1).
```

Consequently the neutral inverse extension escort has

```math
D_\infty(r_{D,t,t}^{0}\Vert U_n)
\le\lambda\delta_n(t)=O_{\beta,\lambda}(1).
```

This is a genuine dimension-free optimizer theorem.  For a fixed orientation
the actual canonical row is instead biased by the opposite child's sector
parameter `gamma`; the proved bound becomes

```math
D_\infty(r_{D,t,t}^{\gamma}\Vert U_n)
\le\lambda\min\{\delta_n(t)+2|\gamma|, 2tn\}.
```

Thus the dimension-free conclusion is presently neutral (or bounded-bias),
not uniform over all conditioned orientation sectors.  The exact orientation
KL records large simultaneous sector bias, but sublinear orientation KL by
itself does not force `gamma=O(1)`.

DER.2 supplies the other genuinely optimizer-specific fact: every row of an
exact minimizing parent is an optimal same-temperature augmented reinsertion
into its own deletion.  This fact is exact, but its state is the complete
function `b -> R_C^aug(b;t,t)` and therefore is not yet a compression.

## 2. Scope classification

| Result | Exact scope | Audit classification |
|---|---|---|
| ES.0--ES.3 | Every finite child pair; hence also optimized children | Rigorous actual-channel regularity and localization, but **not caused by child optimality**.  Linear `J` cannot hide in escaping conditional `D_2` or a vanishing row set.  The antecedent still contains the unresolved interaction/product gap. |
| ES.21--ES.29 | Every exponential tilt of the canonical product | Exact KL chain rules.  They identify TC versus marginal retuning but are an identity, not a reduction. |
| DER.1 / EE.8 | Every signing | The same extension identity in two normalizations. |
| DER.2 / EE.3 | Exact minimizing parent | The same Bellman/reinsertion mechanism.  It retains the full extension-response table. |
| EE.1 | Every actual forward channel | Exact sector-bias algebra.  Optimizer-specific content begins only with the lower envelope EE.2. |
| CSA.1 | Every finite actual channel | Exact one-row curvature sufficient condition.  It is IC.4 plus exponential-family calculus, not yet an optimizer theorem. |
| CSA rank-one/block-parity examples | Generic child priors | Scalable method falsifiers only; they are explicitly not optimizing children. |
| EO.1 | Exact pressure minimizer | Correct entropy ceiling, but only `O(N)` at physical scaling and obtained by summing the existing edge-flip inequality. |
| EO.2 / RP.1 | Fixed finite child orders, `u -> 0` | Essentially the same tangent theorem: the oriented four-spin tensor gives leading row TC and marginal drift starts later.  It is not uniform in order. |
| EO.3 | Order-three exact minimizers | Exact finite Taylor expansion on the physical diagonal near zero temperature parameter.  It has no large-order implication. |
| EO.4 / RP.10--RP.11 | Two certified order-eight minimizer classes | Exact finite actual-minimizer no-go for pressure/entropy/radial data at the response tangent. |
| RP.13 | The same finite classes, `t -> infinity` | Exact tropical leading slopes.  The displayed `t=3` TC/retuning split is numerical, not an interval certificate. |
| EE.4 | Exact optimizing children | Exact orientation-or-density statement, but mere `H_infty=n log 2-o(N)` was already generic from bit oscillation.  The optimizer gain is the sharper neutral `O(1)` density bound. |

The main duplications which should not be counted as separate frontier
advances are DER.1 versus EE.8, DER.2/DER.12 versus EE.3, and EO.2 versus
RP.1.  CSA.1 and the ES decomposition are useful localizations, but neither
removes the quantity it is meant to bound.

## 3. Sharpest actual-child falsifiers

The sharpest statistic-level falsifier is RP.10--RP.13.  Two certified
order-eight pressure-minimizer classes have the same complete signed energy
histogram, hence the same entire radial pressure/entropy data, but different
leading row-dependence curvature (`7 lambda^2` versus `5 lambda^2` in the
stated normalization) and different exact tropical physical-ray slopes
(`25/12` versus `9/5`).  Thus no scalar, homogeneous-flip, or complete radial
consequence of minimization decides the actual row resource.

The sharpest route-specific falsifier is DER.26.  Using only exact thermal
minimizers, an actual deleted row which is optimal for neutral reinsertion is
strictly disfavored by the sector-biased canonical erased-row inverse escort.
Therefore actual-row cavity optimality cannot simply be transferred into the
two-child bridge law.  A sector synchronization, orientation compensation,
or adjacent-order stability theorem is genuinely necessary.

Neither falsifier is asymptotic.  They eliminate exact implications from the
declared coarse data; they do not prove that a large optimizing sequence has
linear row dependence.

## 4. Audit of the proposed smallest missing lemma

The proposed **external-field cavity stability** statement
`Xi_N=o(N)` is logically sufficient and does not mention a target-order
optimizer, so it is noncircular as a bare implication.  It is also more
specific than “understand the actual child law.”  But it is not yet a
strict information reduction: `Xi_N` takes a supremum over the complete
row external-field environment, and evaluating that object can retain the
full bridge response landscape.  Likewise, the whole table
`b -> z_{D,t,t}^gamma(b)` is an exponentially large noisy transform of the
child law.  Calling either table the new state would fail the project's
information criterion.

The **high-transport positive-density resource** is a rigorous conclusion
conditional on linear `J` or on a linear canonical-to-optimal product gap.
It becomes a noncircular SML only after naming a child-only observable which
certifies that resource without computing `J`, the optimal product `p*`, or
the full hybrid laws.  In its current ES.42 form it is a sharp diagnostic
dichotomy, not a deciding theorem.

Accordingly the honest narrowed SML is:

> Find a sector-oriented statistic of exact optimizing children, provably
> lower-information than the complete external-field response, which either
> gives a power-saving bound on the probability-weighted cavity curvature or
> directly exposes a positive density of coherent row information/retuning.

This formulation is noncircular.  Its strictness is a requirement on the
still-missing statistic, not something already established by the current
cavity or extension tables.  EO.2 identifies the first tangent candidate,
the oriented four-spin tensor, but the nonuniform cumulant barrier prevents
its present finite-order form from meeting the requirement at physical
scale.

## 5. Level judgment

No audited claim creates a Level-5-to-6 route.  The results make Level-5
contact substantially sharper: conditional component escape, sparse-row
concentration, neutral extension support, radial summaries, and naive
deletion-to-canonical transfer are now classified.  But no optimizer-specific
statistic selects the `o(N)` versus linear row-resource branch, no pressure
recurrence is completed, and the exact Bellman recurrence retains the full
extension-response state.  A Level-6 claim would therefore be an overclaim.
