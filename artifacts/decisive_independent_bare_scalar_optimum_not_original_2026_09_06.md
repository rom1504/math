# The optimized bare scalar certificate cannot be the original limiting constant

Status: complete elementary conditional implication, pending independent audit.
Dependencies are the director's all-order partition-certificate realization
with deep value H=E, and existence of one bare certificate below 1/2.
This does NOT falsify H=E or the scalar rate-distortion formula itself.

## Statement

For the variance-one ternary source
`nu_p=(1-p)delta_0+(p/2)delta_(1/sqrt(p))+(p/2)delta_(-1/sqrt(p))`, let

```math
C(p,t)={t+p\log2+E_t(\nu_p)\over2t\sqrt p},\qquad
C_* =\inf_{0<p\le1,\ t>0} C(p,t).
```

Assume the known recursive signing construction provides, for every fixed
p,t, its partition-pressure upper certificate with numerator
`t+p log2+E_t(nu_p)`, at all asymptotically large orders as in the director's
construction. Suppose also `C_*<1/2`, as already supplied by a strict rational
certificate. Then

```math
\boxed{\limsup_n M_n/n^{3/2}<C_*.}                       \tag{1}
```

Consequently C_* cannot be a universal lower bound for original signings and
cannot be their limiting value. The proof constructs a strictly better
all-order bound from the SAME optimized pressure certificate by accounting
for the elementary entropy around one extremizer.

## 1. The bare optimum is attained away from all parameter boundaries

The constant label and full-revelation label give respectively

```math
E_t(\nu_p)\ge g_t(1),\qquad
E_t(\nu_p)\ge-H(\nu_p)=-h(p)-p\log2.                 \tag{2}
```

Also `E_t>=-t`, since `g_t(1)>=-t`. Hence

```math
C(p,t)\ge{t+g_t(1)+p\log2\over2t\sqrt p},\qquad
C(p,t)\ge{1-h(p)/t\over2\sqrt p},\qquad
C(p,t)\ge{\sqrt p\log2\over2t}.                     \tag{3}
```

If t tends to infinity, the middle bound gives liminf C>=1/2 uniformly in p.
If p stays bounded away from zero and t tends to zero, the last bound diverges.
If p tends to zero while t is bounded away from zero and infinity, the first
bound diverges, because `t+g_t(1)>0` for t>0. In the remaining joint boundary
p->0,t->0, the exact Gaussian formula gives

```math
t+g_t(1)=t^2+O(t^4).
```

Writing r=t/sqrt(p), the first bound and AM--GM give

```math
C(p,t)\ge\tfrac12[(1+o(1))r+(\log2)/r]
\ge\sqrt{(1+o(1))\log2}\longrightarrow\sqrt{\log2}>1/2
```

in the liminf sense. Finally, if p->1, the preceding bounds exclude t->0 or
t->infinity in a sublevel below 1/2; along the remaining bounded t range,
the middle bound has liminf at least 1/2.

Therefore every sublevel `C<=c_0<1/2` lies in a compact rectangle
`p in [epsilon,1-epsilon]`, `t in [epsilon,T]`.

For completeness, E_t is continuous there. Coupling two variance-one sources
X,Y at W_2 distance d and transporting a near-optimal label through the
coupling gives

```math
|E_t(X)-E_t(Y)|\le t(2d+d^2).
```

Information decreases under the transported channel, and the optimal
conditional mean squared error grows by at most 2d+d^2. The Gaussian reward
g_t is t-Lipschitz in variance. Moreover `|E_t(X)-E_s(X)|<=|t-s|`, since all
posterior average variances are at most one and
`partial_t g_t(v)=-v(1-rho)`. Ternary laws nu_p vary continuously in W_2 on
the indicated p interval. Compactness proves that C_* is attained at some
`0<p_*<1`, `0<t_*<infinity`.

## 2. Every finite-temperature bare certificate has a strict entropy improvement

For any signing on N vertices choose a signed extremizer x* with energy Q.
Flip its coordinates independently with probability delta, `0<delta<1/2`.
This law has entropy N h(delta), and its expected signed quadratic energy is
`(1-2delta)^2 Q`. The finite Gibbs variational principle therefore gives

```math
\log(Z_++Z_-)\ge N h(\delta)
+\beta(1-2\delta)^2 Q.                               \tag{4}
```

No random-sign assumption or independence of selected extrema is involved.
Apply (4) to the signing selected by the construction's partition bound.
In the same normalization that yields C(p,t), the retained vertex count
contributes p h(delta) to the numerator. Thus the all-order upper coefficient
improves to

```math
C_\delta(p,t)=
{t+p\log2+E_t(\nu_p)-p h(\delta)
 \over2t\sqrt p(1-2\delta)^2}
={C(p,t)-{\sqrt p\over2t}h(\delta)\over(1-2\delta)^2}.
                                                               \tag{5}
```

For fixed p,t and positive C, this is strictly smaller than C whenever

```math
{h(\delta)\over\delta}
>{8tC\over\sqrt p}(1-\delta).
```

Such a positive delta exists because `h(delta)/delta -> infinity` as
delta decreases to zero. Use p_*,t_* and then fix one such delta_*. Equation
(5) gives a fixed positive improvement below C_*, proving (1).

## Scope

This excludes identifying the UNCORRECTED optimized finite-temperature
certificate C_* with the original minimum. It does not exclude an appropriately
entropy-corrected, variationally optimized construction from being optimal.
It does not supply the actual limiting value or prove convergence. It also
does not prove the actual ensemble pressure equals its Finner/permanent upper
bound: no such tightness assumption is used.

The gap in (1) is an existence statement unless quantitative localization of
the minimizing parameters and a quantitative delta are additionally supplied.
The all-order sign families are supplied by the already established pressure
realization followed by the exact deterministic extraction (4), not by a
Gaussian or fractional-sign relaxation.
