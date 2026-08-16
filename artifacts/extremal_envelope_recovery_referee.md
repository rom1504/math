# Referee report on `extremal_envelope_recovery.md`

Date: 2026-08-16.

## Verdict

The central implication

```math
\mathrm{EER}\Longrightarrow \lim_n M_n/n^{3/2}\ \text{exists}
```

is correct.  The compactness argument behind (EE.3)--(EE.4), the directed
continuity estimate (EE.6), the weakened covering condition
`gamma_(eta_l) -> 1`, and the deletion proof (EE.10)--(EE.12) all work with
the stated order of quantifiers.

The note should nevertheless not be accepted under its present title or with
its final paragraph unchanged.  EER is not the weakest retained profile
theorem: it is enough to recover one extremizing source profile at each target
order.  Also, the `2 -> 2` rate in (EE.8) can be replaced by uniform
integrability of the source energy products.  Finally, the claimed "exact
falsifier" is not the logical negation of EER and can be avoided by the
allowed order sets.

## 1. Compactness and (EE.4)

This part is valid, subject to spelling out which quantity is uniformly
integrable.

Uniform `2 -> 2` boundedness gives action precompactness.  Its bound is
closed under action convergence: if laws of `(f_j,T_j f_j)` converge to the
law of `(f,Tf)`, then lower semicontinuity of the second-coordinate square,
together with convergence of the bounded first-coordinate square, gives the
corresponding `L^2` inequality in the limit.  Quantitative action continuity
then also passes the cap `Phi <= alpha+eta` to every member of `K_eta`.
Thus `K_eta` is a nonempty compact cluster set with (EE.2).

All laws in the union in (EE.3) satisfy

```math
|x|\le1,
\qquad \int y^2\,d\nu\le C_\eta^2.
```

Consequently they form a uniformly tight family, since

```math
\nu\{|y|>R\}\le C_\eta^2/R^2.
```

Prokhorov compactness makes the closure in (EE.3) compact.  The support
condition is closed, and the second-moment upper bound passes to a weak limit
by lower semicontinuity.  The energy integrand is handled separately: the
relevant family, `xy`, is uniformly integrable because

```math
\int |xy|1_{\{|xy|>R\}}\,d\nu
 \le \int |y|1_{\{|y|>R\}}\,d\nu
 \le C_\eta^2/R.                                         \tag{R.1}
```

Hence `int xy` is continuous along the closure, proving (EE.4).  The sentence
at lines 53--54 is correct only if "uniform integrability" refers to `xy`
(or to `y` in `L^1`).  A bounded second moment does **not** by itself make
the squares `y^2` uniformly integrable; their upper bound passes by lower
semicontinuity instead.  The note should include (R.1) to remove this
ambiguity.

## 2. Directed envelope continuity

Equation (EE.6) has the correct direction, constant, and quantifiers.  Put

```math
K:=\max\{D,C_\eta\}.
```

For a source law `mu` and a matching envelope law `nu` at LP distance
`epsilon`, the truncation calculation gives, for every `R>0`,

```math
\left|\int xy\,d\mu-\int xy\,d\nu\right|
 \le {D^2+C_\eta^2\over R}+(3R+1)\epsilon
 \le {2K^2\over R}+(3R+1)\epsilon.                       \tag{R.2}
```

Taking `R=K/sqrt(epsilon)` yields

```math
\left|\int xy\,d\mu-\int xy\,d\nu\right|
 \le 5K\sqrt\epsilon+\epsilon.
```

The cases `K=0` and `delta=0` follow directly or by an auxiliary positive
error.  The LP distance is at most one, so no missing large-`delta` case is
needed.  Since (EE.5) gives a match for each source law at distance
`delta+o(1)`, (EE.4), the triangle inequality for the absolute energy, and a
supremum over source laws give precisely (EE.6).

It is legitimate for the target phase to depend on the source profile.  The
proof is pointwise in the source law and uses from its target only a moment
bound and the scalar cap (EE.4); it never combines two target laws or asks
that they arise from one operator.  Because (EE.3) is a closure, a completely
literal formulation should say that a target law is a limit of laws from
possibly varying cluster objects, unless closedness of the profile graph is
invoked.  This does not affect the proof.

## 3. The covering condition and the convergence proof

Fix `l`.  Equations (EE.6) and (EE.8), including `delta_m -> 0`, give
(EE.10).  For `m=s_(eta_l)(N)`, unboundedness of the order set implies
`m -> infinity` as `N -> infinity`; repetition of some values of `m` causes
no problem.  Positivity then gives

```math
\limsup_N
 \left({s_{\eta_l}(N)\over N}\right)^{3/2}
 \Phi(T_{A_{s_{\eta_l}(N)}})
 \le \gamma_{\eta_l}^{3/2}(\alpha+\eta_l),
```

which proves (EE.12).  Since `gamma_(eta_l) -> 1` and `eta_l -> 0`, its right
side tends to `alpha`.  This establishes the claimed convergence.  The
sequential quantifiers in lines 116--117 are essential and are stated
correctly: first take `m -> infinity` at fixed `l`, then take `l -> infinity`.

The displayed covering hypothesis is sufficient and is a genuine weakening
of upward ratio-density at every fixed tolerance.  If the heading "weakest
retained order requirement" is intended literally, the exact condition
delivered by this proof is only

```math
\liminf_{l\to\infty}
 \gamma_{\eta_l}^{3/2}(\alpha+\eta_l)=\alpha,             \tag{R.3}
```

after passing to the favorable subsequence of tolerances.  When `alpha>0`,
(R.3) and `eta_l -> 0` force a subsequence with `gamma_(eta_l) -> 1`, so the
form in the note is essentially exact.  Without stating or citing
`alpha>0`, the zero case permits the weaker requirement
`gamma_(eta_l)^(3/2) eta_l -> 0`, even with `gamma_(eta_l)` not tending to
one.  Thus either cite the standard positive lower bound for `alpha`, or
replace the literal minimality wording by "a sufficient weakened order
requirement."

## 4. A weaker continuity hypothesis: `EER_UI`

The global source operator norm and the rate
`max(D_m,C_eta)*sqrt(delta_m) -> 0` are not intrinsically needed.  Weak
profile matching only has to preserve the unbounded observable `(x,y) -> xy`.
A clean replacement is asymptotic uniform integrability of the source
profile products.

For fixed `l`, define

```math
\mathfrak U_l(R):=
\limsup_{\substack{m\to\infty\\m\in\mathcal N_{\eta_l}}}
\ \sup_{\mu\in\mathcal S_1(T_{A_m})}
 \int |xy|1_{\{|xy|>R\}}\,d\mu.
```

Consider the following variant.

> **Uniform-integrability envelope recovery (`EER_UI`).**  There are
> `eta_l -> 0`, unbounded order sets `N_(eta_l)`, and exact hollow signings
> `A_m` on those sets such that, at each fixed `l`,
>
> ```math
> \partial_1(T_{A_m},\mathcal E_{\eta_l})\longrightarrow0,
> \qquad
> \lim_{R\to\infty}\mathfrak U_l(R)=0,                   \tag{R.4}
> ```
>
> and `gamma_(eta_l) -> 1`.

No source `2 -> 2` bound occurs in (R.4).  To prove the implication, suppose
uniform energy closeness failed.  Choose offending source laws `mu_m` and
matching `nu_m in E_(eta_l)` with `d_LP(mu_m,nu_m) -> 0`.  Compactness of the
envelope gives a subsequence `nu_m -> nu`, and then also `mu_m -> nu`.
Products under the target laws are uniformly integrable by (R.1), while
(R.4) gives uniform integrability under the source laws.  The usual weak
convergence plus uniform-integrability lemma therefore gives

```math
\int xy\,d\mu_m-\int xy\,d\nu_m\longrightarrow0
```

uniformly over possible source choices.  Using (EE.4) recovers (EE.10), and
the deletion proof is unchanged.

Output uniform integrability,

```math
\lim_{R\to\infty}\limsup_m\sup_{\mu}
 \int |y|1_{\{|y|>R\}}\,d\mu=0,
```

is an easier sufficient condition, but is stronger than product uniform
integrability because `|x|` can suppress a large output.  Product uniform
integrability is the natural Vitali condition for the actual energy
observable.

This is materially weaker operator information than (EE.8), not merely a
rephrased spectral bound.  Indeed, (EE.8) implies (R.4): under a Strassen
coupling, on the good event

```math
|XY-X'Y'|\le\delta_m(1+|Y'|),
```

while on the exceptional event Cauchy--Schwarz bounds the source and target
contributions by `D_m sqrt(delta_m)` and
`C_(eta_l) sqrt(delta_m)`.  Thus the source products are uniformly `L^1`
close to the uniformly integrable envelope products.  Conversely, product
uniform integrability need not control the source operator norm or even the
second moments of its outputs.

There is one qualification: (R.4) is precisely scalar **tail** control for
the observable `xy`.  It retains substantially less operator structure than
a spectral estimate.  It is not the circular scalar conclusion
`Phi(T_A_m) <= alpha+o(1)`, however: convergence of signed integrals does not
imply uniform integrability of their absolute tails, nor does it imply the
profile matching in (R.4).

## 5. A still weaker noncircular profile theorem

EER controls every law in every source one-profile, but the convergence proof
uses only a law attaining (or nearly attaining) `Phi`.  Since each `A_m` is a
finite hollow matrix, coordinatewise affinity supplies a Boolean
`f_m in {+1,-1}^m` with

```math
\left|\langle f_m,T_{A_m}f_m\rangle\right|
=\Phi(T_{A_m}).
```

Write

```math
\mu_m:=\operatorname{Law}(f_m,T_{A_m}f_m),
\qquad
b_m^2:=\int y^2\,d\mu_m,
\qquad
\delta_m^*:=\inf_{\nu\in\mathcal E_{\eta_l}}d_{LP}(\mu_m,\nu).
```

The following **extremizing-profile envelope recovery** is sufficient:

```math
\delta_m^*\to0,
\qquad
\max\{b_m,C_{\eta_l}\}\sqrt{\delta_m^*}\to0             \tag{R.5}
```

at fixed `l`, together with the same order condition
`gamma_(eta_l) -> 1`.  Applying (R.2) just to `mu_m` gives

```math
\Phi(T_{A_m})
 \le \alpha+\eta_l
    +5\max\{b_m,C_{\eta_l}\}\sqrt{\delta_m^*}
    +\delta_m^*,
```

so the proof of convergence proceeds verbatim.  This drops both the
supremum over all source laws and the global operator norm.  EER implies
(R.5), because `delta_m^* <= delta_m` and `b_m <= D_m`.

An even weaker `UI` form replaces the second condition in (R.5) by uniform
integrability of the single sequence of products under `mu_m`.  Then weak
closeness to the envelope again passes the extremal energy integral.  This is
the weakest clean profile-level statement exposed by the present proof: one
matched extremizing law with enough tail control, on sufficiently covering
orders.

This extremizing-profile hypothesis is not equivalent to the desired scalar
convergence.  Knowing only the values of `Phi(T_{A_m})` does not force the
laws `mu_m` to approach the independently constructed envelope, and it does
not force their absolute product tails to be uniformly integrable.  It is
therefore a noncircular profile target, although it is deliberately tailored
to the energy observable.

## 6. Strictness, circularity, and the purported falsifier

The comparison with selected-phase recovery is valid as a one-way
implication, provided the selected cluster belongs to the envelope being
used:

```math
\partial_1(S,\mathcal E_\eta)
\le \partial_1(S,T)
\qquad(T\in K_\eta).
```

Together with the weaker covering condition, this confirms that EER imposes
fewer coherence obligations.  The prose should say exactly that.  The phrase
"strictly less information" is not established as a mathematical strictness
claim: an envelope may contain laws from several phases whereas a selected
object contains coherent profiles from one phase, so the two data packages
are not ordered without defining a forgetful map and exhibiting
noninjectivity.  "Forgets phase compatibility and all joint profiles" is
accurate and sufficient.

The circularity sentence at lines 164--166 is correct only with the
qualification already present in `minimal_all_order_action_recovery.md`:
conditional on deletion and the availability of exact minimizers, the
objective-only near-covering statement is equivalent to the missing
`limsup <= liminf` assertion.  It is methodologically circular if its
witnesses are selected by assuming that assertion.  EER itself is a
distributional condition, but formal existential quantifiers cannot encode
that its witnesses were constructed "independently"; that is a requirement
on a future proof, not on the logical statement.

Lines 168--171 do not give the exact falsifier of EER.  Infinitely many bad
orders can simply be omitted: for example, isolated bad orders can be
covered by the next good order while `s(N)/N -> 1`.  Moreover, failure of EER
may arise from the norm/error rate, from multiplicatively large gaps, or from
the absence of any favorable null tolerance sequence, without a fixed
positive LP separation.  The quantifier "every bounded-operator exact
signing" also does not match (EE.8), which permits varying and even growing
`D_m`.  That paragraph should be labeled a strong sufficient obstruction to
all-order recovery at a fixed tolerance, or deleted; it is not a logical
negation of EER.

## Recommended disposition

Retain (EE.3)--(EE.12) after adding the explicit tightness/UI argument (R.1)
and a one-line derivation of (R.2).  Change the title and section headings
from "weakest" to "weaker" unless the extremizing-profile and UI variants
are incorporated.  Qualify the strictness and circularity claims, and remove
or substantially rewrite the final falsifier paragraph.
