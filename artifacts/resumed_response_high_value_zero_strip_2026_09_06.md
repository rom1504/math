# A high-value creation mask with an exact zero-slack strip

Date: 2026-09-06. Author: resumed response track.
Status: independently verified by the director's proof reconstruction and
exact rational replay. This is a counterexample to the proposed general-mask
zero-strip slack statement, not a new original Boolean lower bound.

## 1. Statement

In the actual marked-tree Gaussian creation space, there exists an even mask
`0<=H<=1` such that, with `W=UH`,

```
J(H)=E|W|(1-H) > .4301875 > .43,
H=1 whenever |W|<=1/1000.
```

Consequently **even the high-value hypothesis `J>=.43` does not imply any
positive zero-strip slack**, either pointwise for every target or uniformly
over those targets. The mask uses a contraction in the actual creation space;
it is not an arbitrarily prescribed Gaussian covariance model.

The exact executable certificate is
`computations/resumed_response_strip_certificate_2026_09_06.py`. It uses
the existing outward rational interval primitives, exact Gaussian Hermite
products, and literal Darboux upper sums. Its decisive bounds are

```
||(U^{-1}V)_-||_2 < .028006147076777 < .03,
sup ||T(W)-W0||_2^2 < .000467767224 < .000529,
sup_D ||D T(W)[D]||_2^2 / ||D||_2^2 < .893607225 < 1,
J(H) >= .4301875366539460547160317001386721753977010488234.
```

## 2. The fixed banked input and the new ingredient

Use exactly the 21-anchor instance from
`fresh_finite_anchor_fixed_point_2026_09_05.md`. Write

```
alpha=361/500,
H0=1{|V|<=alpha}, W0=UH0, g=U^{-1}V,
p0=EH0, c0=E[VW0].
```

The existing exact construction gives

```
Var V=1, ||g||_2=1,
.5297 < p0 < .5298, c0 > .7004,
J(H0)>.43065, 2 phi(alpha)<.62.
```

The improvement over the earlier attempted strip-filling argument is the
strong one-sided estimate

```
||g_-||_2 < .03.                                           (1)
```

For every nonnegative feedback increment `K` supported outside `H0`, it gives

```
E[V UK]=E[gK] >= -.03 ||K||_2.                            (2)
```

The much weaker covariance-only estimate uses
`sqrt(Var(V|W0))>.27` in place of `.03` and does not close the invariant
domain. The small negative part, not merely a small L2 residual, is essential.

### Proof and certificate of (1)

Let `G` denote the 21 independent Gaussian anchors and `Z` their independent
unit innovation from the banked fixed point. In that notation

```
V=rho.G+sZ, s^2=1-||rho||^2,
chi=1{|rho.G+sZ|<=alpha},
chi_D = total-Gaussian-Hermite-degree-200 projection of chi.
```

Let `R_a` be the resolvent in the innovation variable alone, with `a=17/5`:

```
R_a f(G,Z) = integral_0^1 t^(a-1)
               E[f(G,tZ+sqrt(1-t^2)Z') | G,Z] dt.
```

It acts on innovation Hermite degree `ell` by multiplication by
`1/(a+ell)` and preserves nonnegativity. Let `h_T` be the 21 normalized
removed anchor features, `c_T=E[chi h_T]`, and let `N` be the squared
norm of the unnormalized, anchor-deleted degree-200 resolvent. Then the
actual inverse image of `V` is exactly

```
g=P+A R_a chi_D,
A=s/sqrt(N)>0,
P=sum_T delta_T h_T,
delta_T=rho_T-(A/a)c_T.
```

All removed anchor features have innovation degree zero, which is why their
resolvent denominator is `a`. Thus

```
g=P+A R_a chi - A R_a(chi-chi_D),
||g_-||_2 <= ||P_-||_2 + A ||R_a(chi-chi_D)||_2.             (3)
```

Write `P=delta_0+R`, where `ER=0`. The exact coefficient calculation gives
`delta_0>.09659`. The elementary pointwise inequality

```
(-delta_0-R)_+^2 <= R^4/(16 delta_0^2)
```

follows by maximizing `(x-delta_0)^2/x^4` for `x>=delta_0`; its maximum
occurs at `x=2 delta_0`. Therefore

```
||P_-||_2 <= sqrt(E R^4)/(4 delta_0) < .010106746221.         (4)
```

The program computes `ER^4` exactly in outward intervals. It first expands
`R^2` in the normalized multivariate Hermite basis and then sums coefficient
squares. There are 207 resulting terms. The only product formula used is

```
h_m h_n = sum_(0<=k<=min(m,n))
 k! binom(m,k) binom(n,k)
 sqrt((m+n-2k)!/(m!n!)) h_(m+n-2k).
```

For the resolvent tail, let `beta_d=E[chi h_d(V)]`, and let
`L_d~Binomial(d,s^2)`. The multinomial identity gives

```
||R_a(chi-chi_D)||_2^2
 = sum_(d>200, d even) beta_d^2 E[(a+L_d)^(-2)]
 <= [p0-sum_(d<=200, d even) beta_d^2]
       E[(a+L_202)^(-2)].                               (5)
```

The inequality uses the coupling `L_(d+1)=L_d+Bernoulli(s^2)`, so the
expectation in (5) decreases with `d`. Every quantity in its upper bound
is a finite exact expression. The certificate gives

```
A ||R_a(chi-chi_D)||_2 < .017899400857.
```

Combining this with (3)--(4) proves (1). Positivity applies to the **full**
resolvent `R_a chi`; the finite Hermite truncation is not presumed positive.

## 3. A convex invariant domain for strip feedback

Set

```
e=.03, r=.023, tau=.001, b=.12,
q(w)=1                            for |w|<=tau,
q(w)=(b-|w|)/(b-tau)             for tau<|w|<b,
q(w)=0                            for |w|>=b.
```

On first Gaussian chaos define

```
T(W)=U[H0+(1-H0)q(W)].
```

Use the closed convex domain

```
D={W0+Delta:
       E[W0 Delta]=0,
       ||Delta||_2<=r,
       E[V Delta]>=-er }.
```

It is nonempty, since `W0` belongs to it, and complete in the L2 metric.
For any `W` in this domain, put `v=EW^2` and `c=EVW`. Then

```
p0<=v<=p0+r^2<=vmax=.530329,
c>=c0-er>=cmin=.69971.
```

Write `K=(1-H0)q(W)`. The feedback image satisfies

```
T(W)-W0=UK,
E[W0 UK]=E[H0 K]=0,
||UK||_2^2=EK^2,
E[V UK]>=-e ||K||_2.                                    (6)
```

Therefore an upper bound `EK^2<=r^2` makes the domain invariant.

## 4. Uniform Gaussian envelopes and every derivative direction

For `0<=w<=b`, define the conditional thresholds

```
sigma(v,c)=sqrt(1-c^2/v),
u(v,c,w)=(alpha-cw/v)/sigma(v,c),
a(v,c,w)=(alpha+cw/v)/sigma(v,c),
t0(v,c,w)=Phi(-u)+Phi(-a),
t2(v,c,w)=t0+u phi(u)+a phi(a).
```

Here `t0=P(|V|>alpha | W=w)`. If
`Y=(V-cW/v)/sigma(v,c)`, then
`t2=E[Y^2 1{|V|>alpha}|W=w]`.

Both thresholds are positive. Differentiating them shows that they increase
with `c` and decrease with `v` throughout the stated parameter domain.
For `u`, the relevant derivative signs reduce to

```
c alpha>w,
w(2v-c^2)<alpha v c,
```

which hold even under the conservative bounds
`b<cmin alpha` and `2b vmax<alpha(.5297)cmin`.
The threshold `a` has the same signs directly. Since
`Q(x)` and `Q(x)+x phi(x)` decrease on positive arguments, both `t0,t2`
are consequently bounded above by their values at `(vmax,cmin)`.

The density of `W` on this interval is bounded above by

```
d(w)=phi(w/sqrt(.5297))/sqrt(.5297),
```

because `w^2<.5297<=v` and the centered Gaussian density decreases in
variance there. Below, write `t0(w),t2(w)` for the worst-parameter envelopes.

These give

```
EK^2 <= A0:=2 integral_0^b d(w)t0(w)q(w)^2 dw,
E[|W|K] <= A1:=2 integral_0^b d(w)t0(w)w q(w) dw.           (7)
```

For the contraction, let `D1` be **any** first-chaos direction, not merely
an innovation independent of the mask. Decompose it in the orthonormal
Gaussian coordinates `(X,Y)=(W/sqrt(v),(V-cW/v)/sigma)` plus an independent
Gaussian remainder. Conditional Gaussian second moments show

```
E[(1-H0)q'(W)^2 D1^2] <= Lambda ||D1||_2^2,

Lambda := 2/(b-tau)^2 integral_tau^b
              d(w)[(w^2/.5297)t0(w)+t2(w)] dw.             (8)
```

For clarity, the two-dimensional weighted covariance matrix is positive
semidefinite, so its largest eigenvalue is at most its trace. That trace
is the integrand in (8) before the envelope bounds. The independent
remainder has coefficient `E[(1-H0)q'(W)^2]`, also dominated by the trace
because `t2>=t0`. Thus (8) includes directions correlated with the rare
unmasked conditional tail; none are discarded.

The piecewise linear `q` is Lipschitz. Along the segment between any two
points of the convex domain, Gaussian null boundary sets and the ordinary
fundamental theorem for Lipschitz functions give, using Jensen in the
segment parameter,

```
||T(W)-T(W')||_2^2 <= Lambda ||W-W'||_2^2.                  (9)
```

### Exact numerical enclosure of these bounds

The program uses one interval on `[0,tau]` and 128 equal intervals on
`[tau,b]`. On this range `t0(w),t2(w)` increase, while `d(w)` and `q(w)`
decrease. For `t2`, this follows since `u>=sqrt(2)` and the function
`x^2 phi(x)` decreases for `x>=sqrt(2)`; this endpoint condition is also
checked exactly. Every bin uses the increasing factors at its right endpoint
and decreasing factors at its left endpoint. Thus these are literal upper
Darboux sums, not approximate quadrature with an inferred error sign.

The exact replay gives

```
A0 < .000467767224 < r^2=.000529,
Lambda < .893607225 < 1,
A1 < .000034663347 < .00004.                              (10)
```

Equations (6), (9), and (10) prove that `T` is a strict contraction of its
complete invariant domain. Banach's theorem supplies `W*=T(W*)` there.

## 5. Value and exact vanishing of strip slack

Set `H=H0+(1-H0)q(W*)`. It is even, lies in `[0,1]`, and satisfies
`UH=W*` exactly. Moreover `q=1` on `[-tau,tau]`, so `H=1` there.

For a centered Gaussian pair `(V,W)` of variances `(1,v)` and covariance
`c>=0`, write

```
B(v,c)=E[|W|1{|V|>alpha}].
```

It increases with `v` at fixed `c`. Gaussian covariance differentiation,
or differentiation of the elementary folded-normal formula, gives

```
0<=partial_c B(v,c)
 =2 phi(alpha)[2 Phi(c alpha/sqrt(v-c^2))-1]
 <=2 phi(alpha)<.62.
```

Since `v>=p0` and `c>=c0-er`, (7) and (10) imply

```
J(H)=B(v,c)-E[|W*|K]
 >=J(H0)-.62 er-A1
 >.43065-.62*.03*.023-.000034663347
 >.4301875366539.
```

The construction has strictly smaller value than the banked original
mask, which is immaterial: its value remains above `.43` while its local
slack vanishes identically. It does not improve the original lower bound,
does not rule out special slack at actual marked-mask maximizers, and
does not settle convergence or all-order upper recovery.

## 6. Reproduction and independent audit obligations

Run

```
.venv/bin/python -B computations/resumed_response_strip_certificate_2026_09_06.py
```

The decisive mathematical checks are: the full-resolvent positive
decomposition of `g`; the normalized Hermite fourth-moment multiplication;
the binomial monotonicity in the tail; the convex one-sided invariant domain;
the **all-direction** derivative trace bound; and the parameter-monotone
Gaussian envelopes. A saved prior audit is not substituted for these checks.
The director has independently reconstructed the domain and parameter
derivatives, all-direction trace bound, positive-resolvent decomposition,
and exact Fraction replay.

For a single reproducible saved report, the program also accepts
`--output computations/results/resumed_response_strip_certificate_2026_09_06.json`.
