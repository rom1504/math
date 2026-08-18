# Interaction-curvature interpolation for the actual child escort

Status: **rigorous actual-law structural theorem**.  The canonical
row-erased product `r` and the full negative-disorder escort `q` are joined
by an exact one-parameter family.  Every law on this path retains the same
uniform conditional row Renyi-two bound as the endpoints.  The directed
error of the canonical product is exactly the integrated variance of the
collision--cavity interaction along this regular path.

This theorem does not decide the asymptotic branch.  It turns a linear
**canonical-certificate** error into an explicit extensive-fluctuation
alternative inside the actual child channel, without invoking conference
surrogates or a parent-order optimizer.  This distinction matters:
`D(r||q)` upper-bounds, but need not equal, the directed distance from `q`
to the best row-product law.

## 1. The hybrid path

Use the notation of CR.1.  Thus `p` is the full forward bridge likelihood,
`p_i` are its exact row marginals,

```math
h(B)=\log p(B)-\sum_i\log p_i(B_i),                   \tag{IC.1}
```

and

```math
{dq\over dU}\propto p^{-\lambda},
\qquad
{dr\over dU}\propto\prod_i p_i^{-\lambda}.           \tag{IC.2}
```

For `0<=s<=lambda`, define

```math
{dq_s\over dr}(B)
={e^{-s h(B)}\over E_r e^{-s h}},                    \tag{IC.3}
```

or equivalently

```math
{dq_s\over dU}(B)
\propto p(B)^{-s}\prod_i p_i(B_i)^{-(\lambda-s)}.    \tag{IC.4}
```

Then `q_0=r` and `q_lambda=q`.

## 2. Exact curvature identity

**Theorem IC.1 (canonical error is integrated interaction variance).**  Put

```math
K(s)=\log E_r\exp\{-s(h-E_rh)\}.                     \tag{IC.5}
```

For every finite actual child pair,

```math
K(0)=K'(0)=0,
\qquad
K''(s)=\operatorname{Var}_{q_s}(h),                  \tag{IC.6}
```

and hence

```math
\boxed{
\mathcal J=D(r\Vert q)
=K(\lambda)
=\int_0^\lambda(\lambda-s)
        \operatorname{Var}_{q_s}(h)\,ds.}            \tag{IC.7}
```

Consequently:

```math
\sup_{0\le s\le\lambda}\operatorname{Var}_{q_s}(h)=o(N)
\quad\Longrightarrow\quad \mathcal J=o(N),          \tag{IC.8}
```

whereas

```math
\mathcal J\ge\eta N
\quad\Longrightarrow\quad
\exists s_N\in[0,\lambda]:
\operatorname{Var}_{q_{s_N}}(h)
\ge {2\eta\over\lambda^2}N.                         \tag{IC.9}
```

*Proof.*  If `psi(s)=log E_r e^(-s h)`, differentiation on the finite cube
gives `psi'(s)=-E_(q_s)h` and
`psi''(s)=Var_(q_s)(h)`.  Since `K(s)=psi(s)+sE_rh`, (IC.6) follows.
The density relation in (IC.3) gives

```math
D(r\Vert q)=\lambda E_rh+\log E_re^{-\lambda h}
=K(\lambda).
```

Twice integrating (IC.6) proves (IC.7).  Equations (IC.8)--(IC.9) follow
because `int_0^lambda(lambda-s)ds=lambda^2/2`. `square`

Here and below `N=m+n`, and the asymptotic statements take `lambda` fixed.
Equation (IC.7) also records the exact logical strength of the curvature
criterion:

```math
\mathcal J=o(N)
\quad\Longleftrightarrow\quad
\int_0^\lambda(\lambda-s)\operatorname{Var}_{q_s}(h)\,ds=o(N).
                                                               \tag{IC.7a}
```

The uniform condition in (IC.8) is a convenient **stronger sufficient
condition**, not an equivalent characterization.

## 3. The whole path has tight conditional row complexity

**Theorem IC.2 (uniform regularity of the interaction path).**  Let the
channel amplitude be `u`.  For every `s in [0,lambda]`, every row order,
every prefix value, and every row `i` of length `n`,

```math
\boxed{
D_2(q_s(R_i\mid R_{<i})\Vert U_n)
\le n\log(1+\tanh^2(\lambda u))
\le\lambda^2u^2n.}                                  \tag{IC.10}
```

The same conclusion holds after arbitrary conditioning outside a retained
part of a row and arbitrary marginalization of unmentioned coordinates.
At `u=beta/sqrt(N)` it is uniformly `O_(beta,lambda)(1)` for every split
(indeed the bound is `lambda^2 beta^2 n/N`).

*Proof.*  Flipping a bridge bit changes `log p` by at most `2u`: this is the
likelihood ratio of a binary channel with amplitude `u`.  It changes the
corresponding `log p_i` by at most `2u` and leaves every other row marginal
unchanged.  Therefore the negative log density in (IC.4) has one-bit
oscillation at most

```math
2u\{s+(\lambda-s)\}=2\lambda u.                      \tag{IC.11}
```

Marginalization and conditioning preserve this oscillation bound.  Apply
the conditional cube lemma AC.1 with half-log-odds bound `lambda u` and sum
over the retained row coordinates. `square`

Combining IC.1 and IC.2 gives an exhaustive structural alternative for the
explicit canonical product certificate:

```text
sublinear collision--cavity curvature along the regular hybrid path,
or an extensive interaction fluctuation under an actual-child law whose
every conditional row still has tight Renyi-two complexity.
```

The second branch is an extensive mismatch from the explicit canonical
row product.  It witnesses nonproduct structure, and it cannot be
reclassified as escaping conditional component complexity.  It does **not**
by itself prove that

```math
\inf_{a=\otimes_i a_i}D(a\Vert q)=\Omega(N),
```

because a different row product could be much closer to `q` than `r` is.

## 4. A tilted average-influence criterion

The supremum projective diameter in CC.14 is often too costly.  There is a
different, potentially weaker certificate which averages row influence
under the same hybrid laws.  It is lower-information in form, but without
an additional uniform increment bound its asymptotic criterion is not
logically implied by the projective-diameter criterion.  Let

```math
\bar h_i(B_{-i})=E_{R_i\sim r_i}h(R_i,B_{-i}),
\qquad \tau(a)=e^a-a-1,                              \tag{IC.12}
```

and, for `s>0`, define

```math
\mathcal A_s
={1\over s^2}E_{q_s}\sum_{i=1}^m
 \tau\{s(h(B)-\bar h_i(B_{-i}))\}.                  \tag{IC.13}
```

At `s=0`, use the continuous value

```math
\mathcal A_0={1\over2}E_r\sum_i
 (h-\bar h_i)^2.                                    \tag{IC.14}
```

**Theorem IC.3 (tilted row-influence bound).**  One has

```math
\boxed{
\mathcal J\le\lambda\int_0^\lambda\mathcal A_s\,ds.} \tag{IC.15}
```

Therefore

```math
\sup_{0\le s\le\lambda}\mathcal A_s=o(N)
\quad\Longrightarrow\quad\mathcal J=o(N),           \tag{IC.16}
```

while `J>=eta N` forces `A_(s_N)>=eta N/lambda^2` for
some `s_N`.

*Proof.*  Apply entropy tensorization on the product row law `r`, followed
by the modified logarithmic Sobolev inequality, to `Z=-h` and the
row-deleted comparison `Z_i=-bar h_i`:

```math
\operatorname{Ent}_r(e^{-s h})
\le E_r e^{-s h}\sum_i
 \tau\{s(h-\bar h_i)\}.                              \tag{IC.17}
```

This is the standard product-space modified log-Sobolev inequality of
Boucheron--Lugosi--Massart, *Annals of Probability* **31** (2003),
doi:`10.1214/aop/1055425791`; it allows any comparison `Z_i` independent of
row `i`.

Dividing (IC.17) by `E_r e^(-s h)` and using (IC.5) gives

```math
sK'(s)-K(s)\le s^2\mathcal A_s.                     \tag{IC.18}
```

Since `(K(s)/s)'=(sK'(s)-K(s))/s^2` and
`K(s)/s ->0` as `s downarrow0`, integration proves (IC.15).  The remaining
claims follow by averaging over an interval of length `lambda`. `square`

IC.3 is not a bounded-difference restatement.  It retains the sign and size
of each complete collision--cavity row replacement through `tau`, averages
under the relevant inverse tilt, and needs no worst-case posterior range.
It is still a sufficient criterion, not a proof that actual minimizers have
small influence.  More precisely, if `delta_i` denotes the row range in
CC.14, then

```math
|h-\bar h_i|\le\delta_i(B_{-i}),
\qquad
\mathcal A_s\le {1\over s^2}
 \sum_i\tau\{s\sup_{B_{-i}}\delta_i(B_{-i})\}.       \tag{IC.19}
```

Thus projective synchronization plus a common bounded increment gives the
average-influence criterion (with a constant depending on that bound).
Without the common bound, `sum_i delta_i^2=o(N)` permits a sparse growing
increment whose exponential `tau` cost is not controlled; the two
sufficient hypotheses are then not ordered.

## 5. Relation to the smallest missing lemma

The identity (IC.7) alone is a cumulant reformulation and is not claimed as
a strict reduction.  Its value is that the random variable `h` now has the
child-only collision--cavity representation CC.8--CC.11, while the entire
interpolating family satisfies the uniform structural hypothesis (IC.10).
Thus the next proof obligation can be stated without a product variational
oracle:

> **Actual-child collision-curvature lemma.**  Prove the weighted integrated
> curvature in (IC.7a) is `o(N)` for contracted-temperature minimizing
> children.  A stronger, potentially easier-to-localize sufficient version
> is `sup_(0<=s<=lambda) Var_(q_s)(log G)=o(N)`; a different sufficient
> route is `sup_s A_s=o(N)` in (IC.13).  Alternatively, exhibit a
> child-derived statistic forcing either `Omega(N)` alternative.

This obligation has a smaller **declared output** than evaluating the
optimal row-product shadow: it asks for scalar variances of the explicit
child-only partition `G` under a one-parameter family with known local
regularity, and it removes the product variational oracle.  But (IC.7a) is
exactly equivalent to `J=o(N)`, while `J=o(N)` is only a sufficient (and
possibly stronger) route to the desired optimal-shadow conclusion.  Thus
IC.1 is not yet a strict mathematical reduction.  It becomes one only if
the variance can be controlled from demonstrably coarser optimizer-specific
child data.  Generic weak-coordinate laws can have linear collision
curvature, as the rank-one example CC.18--CC.24 shows.
