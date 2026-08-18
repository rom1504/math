# Actual-child sector min-entropy and the diffuse-retuning obstruction

Status: **rigorous actual-minimizer theorem and exact scope statement**.  The
full optimizer contraction (FC.8) gives more than an `O(N)` entropy ceiling:
at contracted temperature, every one-child Gibbs sector has exponentially
small atoms.  Consequently the exact two-child rank-one prior has linear
min-entropy and collision entropy, uniformly over the actual minimizing
children and the relative orientation.

This does not upper-bound posterior retuning.  It proves a different and
useful fact: any posterior retuning supported on a subexponential (or
sufficiently low-rate exponential) latent family necessarily costs linear
KL.  Thus any sublinear-retuning theorem must be genuinely diffuse; it cannot
come from selecting a small catalogue of child spin words.

All logarithms below are natural.

## 1. Optimizer contraction bounds the child ground cap

Let `A` be a pressure-minimizing signing of order `m` at raw temperature
`t>0`, and write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad K_A=\max_x|H_A(x)|,
\qquad d_m={m\choose2}.
```

The augmented normalized partition function is

```math
\overline Z_A(t)=2^{-m}\sum_x\cosh(tH_A(x)).
```

The optimizer contraction (FC.8), at the origin, gives

```math
\log\overline Z_A(t)\le d_m\log\cosh t.             \tag{ME.1}
```

If `x_*` realizes `K_A`, then `-x_*` realizes the same value and
`cosh(tK_A)>=e^(tK_A)/2`.  Their two terms show

```math
\overline Z_A(t)\ge2^{-m}e^{tK_A}.
```

Therefore

```math
\boxed{
{tK_A\over m}
\le \log2+{m-1\over2}\log\cosh t.}                  \tag{ME.2}
```

In the child application `t=beta/sqrt(N)` and `m<=N`, so

```math
\boxed{
{tK_A\over m}\le C_\beta,
\qquad C_\beta:=\log2+{\beta^2\over4}.}              \tag{ME.3}
```

Here only `log cosh t<=t^2/2` was used.  If `m/N->theta`, the sharper
asymptotic constant is `log2+beta^2 theta/4`.

## 2. A Hamming-sphere entropy lemma

For `s in {+-1}`, let

```math
\mu_{A,s}(x)
={e^{tsH_A(x)}\over\sum_z e^{tsH_A(z)}}.             \tag{ME.4}
```

Fix a maximizer `x_*` of `sH_A`, and let `x_*^R` be obtained by flipping the
vertices in `R`.  Uniformly over all `r`-subsets `R`, every edge has exactly
one endpoint flipped with probability `2r(m-r)/(m(m-1))`.  Hence

```math
{1\over{m\choose r}}\sum_{|R|=r}sH_A(x_*^R)
=\left(1-{4r(m-r)\over m(m-1)}\right)sH_A(x_*).
                                                               \tag{ME.5}
```

Jensen's inequality applied on this Hamming sphere gives the following
general lemma.

**Lemma ME.1 (quadratic Hamming-sphere min-entropy).**  For every quadratic
Boolean form, every `t>0`, every sector `s`, and every `1<=r<=m-1`,

```math
\boxed{
\|\mu_{A,s}\|_\infty
\le {m\choose r}^{-1}
 \exp\left\{{4tK_A r(m-r)\over m(m-1)}\right\}.}      \tag{ME.6}
```

Indeed, the denominator in (ME.4), restricted to the indicated sphere, is
at least

```math
{m\choose r}\exp\left\{
t\left(1-{4r(m-r)\over m(m-1)}\right)sH_A(x_*)
\right\},
```

while the numerator of the largest atom is `exp(tsH_A(x_*))`.

Already the one-flip sphere gives a useful exact finite-order estimate.
Combining `r=1` in (ME.6) with (ME.3) yields

```math
\boxed{
\|\mu_{A,s}\|_\infty
\le {e^{4C_\beta}\over m}
={16e^{\beta^2}\over m}.}                            \tag{ME.6a}
```

Thus no asymptotic choice of radius is needed merely to obtain a uniform
power atom bound.  The positive-density spheres below upgrade it to an
exponential bound.

## 3. Exponential sector and rank-one collision bounds

Let

```math
h(q)=-q\log q-(1-q)\log(1-q)
```

be binary entropy.  Combining (ME.3), (ME.6), and the type-class bound
`binom(m,r)>=exp(mh(r/m))/(m+1)` gives

```math
\begin{aligned}
-{1\over m}\log\|\mu_{A,s}\|_\infty
\ge{}&h(q_m)
-4C_\beta q_m(1-q_m){m\over m-1}
-{\log(m+1)\over m},\\
&q_m={r\over m}.                                    \tag{ME.7}
\end{aligned}
```

Define the positive constant

```math
\eta_\beta
=\sup_{0<q<1/2}\{h(q)-4C_\beta q(1-q)\}>0.           \tag{ME.8}
```

Positivity is elementary because `h(q)/q->infinity` as `q` decreases to
zero.  In fact, take `q_beta=exp(-4C_beta)`.  The elementary inequality

```math
h(q)\ge q\log(1/q)+q(1-q)
```

gives

```math
\boxed{
\eta_\beta\ge e^{-4C_\beta}={e^{-\beta^2}\over16}.}  \tag{ME.8a}
```

For a completely explicit eventually uniform exponent put

```math
\underline\eta_\beta={e^{-1-\beta^2}\over16}<\eta_\beta. \tag{ME.8b}
```

Taking `r=floor(qm)` in (ME.7) proves, uniformly in
the sector and in the actual minimizing child,

```math
\boxed{
\liminf_{m\to\infty}-{1\over m}
 \log\|\mu_{A,s}\|_\infty\ge\eta_\beta.}            \tag{ME.9}
```

In particular, for all sufficiently large `m`, uniformly in the actual
child and sector,

```math
\|\mu_{A,s}\|_\infty\le e^{-\underline\eta_\beta m}. \tag{ME.9a}
```

Now take actual minimizing children `A,D` of orders `m,n`, with `m+n=N`,
and fix the relative orientation `epsilon`.  Their exact zero-bridge law is

```math
\nu_\epsilon(s,x,y)
=\pi_s^{(\epsilon)}\mu_{A,s}(x)\mu_{D,\epsilon s}(y),
\qquad Q_{ij}=sx_iy_j.                                \tag{ME.10}
```

For fixed `Q` and fixed `s`, exactly two pairs `(x,y)` map to `Q`, related by
`(x,y)->(-x,-y)`.  Thus the induced actual-child prior `mu_epsilon` obeys

```math
\|\mu_\epsilon\|_\infty
\le2\max_s\|\mu_{A,s}\|_\infty
             \|\mu_{D,\epsilon s}\|_\infty.          \tag{ME.11}
```

In particular, (ME.6a) gives the fully finite estimate

```math
\boxed{
\|\mu_\epsilon\|_\infty
\le {512e^{2\beta^2}\over mn}.}                      \tag{ME.11a}
```

This factor has been checked directly.  The map
`(s,x,y)->Q=sxy^T` has exactly four preimages: two choices of `s`, and for
each `s` two simultaneous gauges `(x,y),(-x,-y)`.  The two fixed-`s` weights
are equal.  Summing them gives the factor `2` in (ME.11), while summing the
sector weights uses `sum_s pi_s=1`; there is no missing factor `2` or `4`.
For the one-child augmented projective law, quotienting `x~-x` similarly
multiplies a conditional sector atom by exactly `2`, which does not change
the exponential rate.

For any sequence of comparable splits this yields

```math
\boxed{
\liminf_{N\to\infty}-{1\over N}
 \log\|\mu_\epsilon\|_\infty\ge\eta_\beta,
\qquad
\sum_Q\mu_\epsilon(Q)^2
\le\|\mu_\epsilon\|_\infty
\le e^{-\eta_\beta N+o(N)}.}                         \tag{ME.12}
```

Equivalently, the actual rank-one prior has

```math
H_\infty(\mu_\epsilon),\ H_2(\mu_\epsilon)
\ge\eta_\beta N-o(N).                                \tag{ME.13}
```

Consequently, after increasing the finite threshold if necessary,

```math
\boxed{
\|\mu_\epsilon\|_\infty,
\ \sum_Q\mu_\epsilon(Q)^2
\le e^{-\underline\eta_\beta N}.}                   \tag{ME.13a}
```

No conference, Paley, or generic row surrogate is used.

This settles the previously open literal two-word test (PA.13): the
projective-atom observable `eta_epsilon(A,D;t)` of PA.12 is not merely
`o(1)` but `exp[-Omega_beta(N)]`.  It also upgrades the one-child
projective ceiling PA.7.  Conditional on either augmented sector, projective
atoms have exponentially small mass, so the unconditional augmented
projective law does as well.  There is no contradiction with Proposition
PA.3: that construction satisfies the signs of AC.32 but is explicitly not
the Gibbs law of a quadratic signing.  The Hamming-sphere identity (ME.5) is
the additional Gibbs/quadratic structure missing from that abstract cone.

It also controls the projective caps in PA.13a--PA.13c.  A cap

```math
\left\{[x]:{|\langle x,u\rangle|\over m}\ge1-\delta\right\}
```

is a projective Hamming ball of radius at most `delta m/2`.  For fixed
`0<delta<1`, its number of projective classes is at most

```math
\sum_{j\le\delta m/2}{m\choose j}
=\exp\{m h(\delta/2)+o(m)\}.                         \tag{ME.13b}
```

The projective sector atom is exactly twice the corresponding ordinary
sector atom, because `[x]={x,-x}` and the two weights agree.  Therefore

```math
\boxed{
\kappa_{A,s}(\delta)
\le\exp\{-[\eta_\beta-h(\delta/2)]m+o(m)\}}
                                                               \tag{ME.13c}
```

whenever `h(delta/2)<eta_beta`.  Substitution in the upper half of PA.13c
gives, on comparable splits,

```math
\boxed{
\Xi_\epsilon(\delta)
\le\exp\{-[\eta_\beta-h(\delta/2)]N+o(N)\}.}          \tag{ME.13d}
```

There is no lost orientation factor: the upper comparison is a convex sum
over the two sector weights.  Thus actual optimality excludes not only a
fixed rank-one word but every fixed-radius common-sign projective cluster
whose cap entropy is below the min-entropy rate.  It still does not exclude
a diffuse union of exponentially many such caps.

Explicitly, let `C_delta(u,v)` be any rank-one `delta`-cap in PA.13b.  If a
posterior average `bar mu` assigns it mass at least `a`, binary data
processing and (ME.13d) give

```math
\boxed{
D(\bar\mu\Vert\mu_\epsilon)
\ge a[\eta_\beta-h(\delta/2)]N-o(N)-h(a).}            \tag{ME.13e}
```

Thus one coherent geometric direction cannot carry fixed posterior mass at
sublinear retuning cost whenever `h(delta/2)<eta_beta`.

## 4. What this proves about posterior retuning

Let `bar mu` be any retuned latent law, including the posterior average
arising from any point on the actual negative-disorder path.  If
`bar mu(R)>=a`, data processing to the indicator of `R` gives

```math
D(\bar\mu\Vert\mu_\epsilon)
\ge a\log{1\over\mu_\epsilon(R)}-h(a).               \tag{ME.14}
```

If `|R|<=exp(xi N)` with `xi<eta_beta`, (ME.12) therefore implies

```math
\boxed{
D(\bar\mu\Vert\mu_\epsilon)
\ge a(\eta_\beta-\xi)N-o(N)-h(a).}                   \tag{ME.15}
```

In particular, a retuned law with `D(bar mu||mu_epsilon)=o(N)` places
vanishing mass on every subexponential latent catalogue.  This is a strict
actual-child consequence: the complete rank-one support bound alone allows
arbitrarily heavy atoms and does not imply (ME.15).

The same conclusion has a compact entropy form.  Since

```math
E_{\bar\mu}[-\log\mu_\epsilon]
=H(\bar\mu)+D(\bar\mu\Vert\mu_\epsilon)
\ge H_\infty(\mu_\epsilon),
```

(ME.13) gives

```math
\boxed{
H(\bar\mu)
\ge\eta_\beta N-D(\bar\mu\Vert\mu_\epsilon)-o(N).}  \tag{ME.16}
```

Thus sublinear retuning KL forces a posterior of exponential effective
support, rather than merely ruling out one large atom.

There is also a geometric cluster version.  Work on the latent bit space
`(s,x,y)` before its fourfold gauge quotient, and call the image of a
Hamming ball there a latent Hamming cluster.  A union `R_N` of at most
`exp(kappa N)` such clusters of radius at most `delta N`, with
`0<=delta<1/2`, has

```math
|R_N|\le
4\exp\{[\kappa+h(\delta)+o(1)]N\}.                  \tag{ME.17}
```

Therefore, whenever `kappa+h(delta)<underline eta_beta`,

```math
\boxed{
\mu_\epsilon(R_N)
\le\exp\{-[\underline\eta_\beta-\kappa-h(\delta)-o(1)]N\}.}
                                                               \tag{ME.18}
```

If `bar mu(R_N)>=a`, (ME.14) turns (ME.18) into the same coefficient times
`aN`, up to `o(N)+h(a)`.  In particular, subexponentially many
`o(N)`-radius latent clusters cannot carry fixed posterior mass at
sublinear retuning cost.

The direction is also the precise limitation.  Linear min-entropy does not
upper-bound `D(bar mu||mu_epsilon)`: a posterior may retune diffusely across
an exponentially large family.  Nor does (ME.12) reduce the weighted
transport cross-entropy to `o(N)`.  Therefore AC.32/FC.8 exclude
low-complexity atomic retuning, but a new theorem is still needed to control
**diffuse exponential-rate retuning** under the actual inverse bridge tilt.

## 5. Research consequence

The prior-complexity question has a rigorous answer:

```math
\text{actual minimizing child}
\Longrightarrow
\text{linear sector min-entropy and collision entropy}.
```

Accordingly, no finite, polynomial, or subexponential catalogue of planted
rank-one words can be the missing posterior direction at sublinear KL cost.
The remaining smallest lemma is narrower: determine whether the actual
negative-disorder posterior retunes diffusely over an exponential-rate
family, and if so whether that diffuse retuning has an optimizer-specific
low-dimensional direction.  Atom or collision estimates alone cannot decide
that question.
