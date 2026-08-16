# Extremal-envelope recovery weakens phase-coherent AR

Date: 2026-08-16.

Status: **verified implication**, independently audited.  This note sharpens
the selected-phase statement in `minimal_all_order_action_recovery.md`.  It
does not prove the recovery hypothesis.

## 1. Normalization and purified cluster envelopes

Put

```math
p_n:=\frac{2M_n}{n^{3/2}},
\qquad
\alpha:=\liminf_{n\to\infty}p_n.
```

For every fixed purification tolerance `eta>0`, choose a sequence of exact
hollow signings `B_j` with orders tending to infinity such that

```math
\Phi(T_{B_j})\le\alpha+\eta,
\qquad
\|T_{B_j}\|_{2\to2}\le C_\eta.                         \tag{EE.1}
```

Let `K_eta` be the compact set of action cluster points of this sequence.
Every `T in K_eta` has

```math
\|T\|_{2\to2}\le C_\eta,
\qquad
\Phi(T)\le\alpha+\eta.                                  \tag{EE.2}
```

Let `S_1(T)` be the closed one-profile of `T`, and define the compact
**extremal profile envelope**

```math
\mathcal E_\eta
:=\overline{\bigcup_{T\in K_\eta}\mathcal S_1(T)}.       \tag{EE.3}
```

The closure is in Levy--Prokhorov distance.  Every law `nu` in this envelope
has first coordinate bounded by one, second-coordinate second moment at most
`C_eta^2`, and

```math
\left|\int xy\,d\nu(x,y)\right|\le\alpha+\eta.          \tag{EE.4}
```

Indeed, these statements hold on each `S_1(T)` by (EE.2).  Tightness and the
second-moment upper bound pass to the closure.  For the energy observable,
the required uniform integrability is explicitly

```math
\int |xy|\mathbf1_{\{|xy|>R\}}\,d\nu
\le\int |y|\mathbf1_{\{|y|>R\}}\,d\nu
\le\frac{C_\eta^2}{R}.                                  \tag{EE.4a}
```

Thus the integrals in (EE.4) also pass to the closure.

## 2. Directed distance to an envelope

For an operator `S`, put

```math
\partial_1(S,\mathcal E_\eta)
:=\sup_{\mu\in\mathcal S_1(S)}
  \inf_{\nu\in\mathcal E_\eta}d_{LP}(\mu,\nu).           \tag{EE.5}
```

The target law in (EE.5) may come from a different cluster object for every
source law.  No single phase `T`, even one depending on the order, has to
explain the complete profile of `S`.

The quantitative one-profile argument gives the following envelope form.

> **Envelope continuity.**  If `||S||_(2->2)<=D`, then, for
> `delta=partial_1(S,E_eta)`,
>
> ```math
> \boxed{
> \Phi(S)\le\alpha+\eta
>       +5\max\{D,C_\eta\}\sqrt\delta+\delta.}           \tag{EE.6}
> ```

To prove it, take an arbitrary law in `S_1(S)`, couple it to an envelope law
within `delta+o(1)`, and repeat the good-event/tail truncation proof of
directed one-profile continuity.  Only the two displayed second-moment bounds
enter.  Apply (EE.4), then take the supremum and let the auxiliary error tend
to zero.  Explicitly, if `K=max(D,C_eta)` and the LP error is `epsilon`,
truncation at `R` gives

```math
\left|\int xy\,d\mu-\int xy\,d\nu\right|
\le\frac{2K^2}{R}+(3R+1)\epsilon.                       \tag{EE.6a}
```

Taking `R=K/sqrt(epsilon)` yields (EE.6), with the zero cases obtained by
an auxiliary positive error.

## 3. A weakened order requirement

For an unbounded set of orders `N_eta`, define its upward covering function

```math
s_\eta(N):=\min\{m\in\mathcal N_\eta:m\ge N\},
\qquad
\gamma_\eta:=\limsup_{N\to\infty}\frac{s_\eta(N)}N.      \tag{EE.7}
```

Requiring `gamma_eta=1` at every fixed tolerance is unnecessary.  It is enough
to have a null sequence `eta_l` for which `gamma_(eta_l)->1`.

> **Extremal-envelope recovery (`EER`).**  There are `eta_l->0` and, for
> every `l`, exact hollow signings `A_m` on an unbounded order set
> `N_(eta_l)` such that
>
> ```math
> \delta_m:=\partial_1(T_{A_m},\mathcal E_{\eta_l})\to0,
> \qquad
> \max\{D_m,C_{\eta_l}\}\sqrt{\delta_m}\to0,             \tag{EE.8}
> ```
>
> where `||T_(A_m)||_(2->2)<=D_m`, and
>
> ```math
> \gamma_{\eta_l}\longrightarrow1.                       \tag{EE.9}
> ```

The limit in (EE.8) is taken at fixed `l` along `m in N_(eta_l)`; only after
that is `l` sent to infinity.

The exact order condition used below is that, along a favorable tolerance
subsequence,

```math
\gamma_{\eta_l}^{3/2}(\alpha+\eta_l)\longrightarrow\alpha. \tag{EE.9a}
```

The rigorous positive lower bound for `alpha` makes (EE.9) essentially
equivalent to (EE.9a) when `eta_l->0`.  Stating (EE.9) avoids a needless
case split.

## 4. EER implies convergence

By (EE.6)--(EE.8), at each fixed `l`,

```math
\limsup_{\substack{m\to\infty\\m\in\mathcal N_{\eta_l}}}
\Phi(T_{A_m})\le\alpha+\eta_l.                            \tag{EE.10}
```

For arbitrary `N`, take `m=s_(eta_l)(N)`.  Lossless principal deletion gives
an exact order-`N` signing with cap at most `Q(A_m)`, and therefore

```math
p_N\le\left(\frac mN\right)^{3/2}\Phi(T_{A_m}).          \tag{EE.11}
```

Equations (EE.7), (EE.10), and (EE.11) yield

```math
\limsup_{N\to\infty}p_N
\le\gamma_{\eta_l}^{3/2}(\alpha+\eta_l).                 \tag{EE.12}
```

Let `l->infinity`.  By (EE.9), the right side tends to `alpha`, which is the
liminf by definition.  Hence `p_n`, and equivalently `M_n/n^(3/2)`, converges.

## 5. Removing the source operator bound: `EER_UI`

The source `2 -> 2` bound and the quantitative rate in (EE.8) are sufficient,
not intrinsic.  The exact analytic issue is continuity of the unbounded
energy observable `(x,y) -> xy`.

For fixed `l`, define

```math
\mathfrak U_l(R):=
\limsup_{\substack{m\to\infty\\m\in\mathcal N_{\eta_l}}}
\ \sup_{\mu\in\mathcal S_1(T_{A_m})}
\int |xy|\mathbf1_{\{|xy|>R\}}\,d\mu.                  \tag{EE.13}
```

> **Uniform-integrability envelope recovery (`EER_UI`).**  Replace (EE.8)
> by
>
> ```math
> \partial_1(T_{A_m},\mathcal E_{\eta_l})\to0,
> \qquad
> \lim_{R\to\infty}\mathfrak U_l(R)=0,                 \tag{EE.14}
> ```
>
> and retain the order condition (EE.9a).

This still proves convergence.  Otherwise choose source laws whose energy
integrals violate (EE.10), and choose matching envelope laws.  Compactness of
the envelope gives a convergent target subsequence; directed LP error makes
the source laws converge weakly to the same limit.  Equation (EE.4a) and
(EE.14) give uniform integrability on the two sides, so Vitali convergence
passes the `xy` integrals and contradicts (EE.4).

Condition (EE.8) implies (EE.14) by a Strassen coupling and its
`D_m sqrt(delta_m)` exceptional-event bound.  Conversely, (EE.14) need not
control source output second moments or operator norms.  It is the weakest
natural tail condition for the observable actually used in the objective;
uniform integrability of the outputs `y` would be stronger.

## 6. A formally weaker extremizing-profile version

The convergence proof does not literally use every source profile.  For each
finite `A_m`, choose a Boolean `f_m` attaining `Phi(T_A_m)` and let

```math
\mu_m=\operatorname{Law}(f_m,T_{A_m}f_m),
\qquad
\delta_m^*=\inf_{\nu\in\mathcal E_{\eta_l}}d_{LP}(\mu_m,\nu).
```

It is enough that `delta_m^*->0` and that this single sequence of products is
uniformly integrable, together with (EE.9a).  Equivalently, a quantitative
version may require

```math
\max\left\{C_{\eta_l},
 \left(\int y^2\,d\mu_m\right)^{1/2}\right\}
\sqrt{\delta_m^*}\longrightarrow0.                      \tag{EE.15}
```

This is the formally weakest clean profile statement exposed by the proof.
It is not retained as the constructive target: identifying `mu_m` requires
solving the full target-order Boolean maximum or carrying its ground-state
layer.  Under the project's information criterion, it is target optimization
with an attached profile witness.  `EER_UI` is the weakest **optimizer-free**
profile target found.  Dropping profile matching as well leaves only the
scalar upper recovery statement, which is equivalent to convergence.

## 7. What has and has not been reduced

Compared with selected-cluster `AR_min^->`, EER removes two genuine
coherence obligations:

1. the recovered matrix need not approximate one fixed action phase, or even
   one phase per target order; each individual profile may match a different
   member of the liminf cluster envelope; and
2. target orders need only have asymptotic upward multiplicative gap
   `gamma_eta`, with `gamma_eta->1` as the objective tolerance vanishes.

The retained state is a compact unlabeled set of one-input/one-output laws.
It forgets vertex labels, phase compatibility between different test
functions, reverse inclusion, and every joint profile.  This is a genuine
loss of coherence obligations relative to selected-phase AR; no claim of a
formal information ordering is needed.

This is nevertheless not a proof that EER is easier to establish.  It still
requires excluding a separated outer profile for every bounded coloring of a
target-order signing.  Replacing (EE.5) only by the scalar upper bound in
(EE.6) gives the objective-only recovery statement, which is equivalent to
convergence and is circular as a target.

A strong sufficient obstruction would be an `epsilon>0` and multiplicatively
nonnegligible intervals of target orders on which every admissibly
bounded/tail-controlled signing has a profile law at LP distance `epsilon`
from the envelope.  Isolated bad orders do not falsify EER: the recovery set
may omit them and cover them from the next order.  Failure may also occur via
the tail condition or the covering ratios, so no single finite test is the
logical negation of EER.
