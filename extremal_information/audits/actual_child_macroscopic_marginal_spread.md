# Every macroscopic child marginal has exponential spread

Status: **rigorous actual-child theorem with uniform conditional and marginal
scope**.  The only optimizer input is the scalar cap contraction already
proved in Theorem 37.61.  The theorem rules out frozen or low-rate common
patterns on every prescribed positive-density vertex block under an actual
child sector.  It does not control posterior retuning, approximate row
lifetime at large Hamming radius, or the full balanced-product target.

All logarithms are natural.

## 1. A subset-flip lemma

Let

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad K_A=\max_x|H_A(x)|,
```

and split the vertices as `U union V`, with `|U|=k>=2`.  Write

```math
H_A(u,v)=H_U(u)+H_{U,V}(u,v)+H_V(v).                \tag{MS.1}
```

For `R` uniform among the `r`-subsets of `U`, let `u^R` be obtained by
flipping the coordinates in `R`, and put

```math
a_{k,r}={4r(k-r)\over k(k-1)},
\qquad b_{k,r}={2r\over k}.                         \tag{MS.2}
```

An edge internal to `U` has exactly one endpoint flipped with probability
`2r(k-r)/(k(k-1))`, while a cross edge has its `U` endpoint flipped with
probability `r/k`.  Therefore, exactly,

```math
E_R H_A(u^R,v)
=(1-a_{k,r})H_U(u)+(1-b_{k,r})H_{U,V}(u,v)+H_V(v).
                                                               \tag{MS.3}
```

Both pieces affected in (MS.3) are controlled by the *global* cap:

```math
H_U(u)=E_{v'\text{ uniform}}H_A(u,v'),
\qquad
H_{U,V}(u,v)={H_A(u,v)-H_A(u,-v)\over2},             \tag{MS.4}
```

so

```math
|H_U(u)|\le K_A,
\qquad |H_{U,V}(u,v)|\le K_A.                       \tag{MS.5}
```

The first identity uses that every nonconstant quadratic monomial in `v'`
has uniform mean zero.  The second uses that `H_V(-v)=H_V(v)` and the cross
term changes sign.  Equations (MS.3)--(MS.5) imply, for either sector sign
`s`, every `u,v`, and `t>0`,

```math
E_R\,stH_A(u^R,v)
\ge stH_A(u,v)-t(a_{k,r}+b_{k,r})K_A.                \tag{MS.6}
```

## 2. Uniform conditional and marginal min-entropy

Let

```math
\mu_{A,s}(u,v)={e^{stH_A(u,v)}\over
                        \sum_{u',v'}e^{stH_A(u',v')}}.
```

For fixed `v`, restrict the conditional partition sum over `u'` to the
Hamming sphere `{u^R:|R|=r}`.  Jensen and (MS.6) give

```math
{1\over {k\choose r}}\sum_{|R|=r}e^{stH_A(u^R,v)}
\ge \exp\{E_R stH_A(u^R,v)\}
\ge e^{stH_A(u,v)-t(a_{k,r}+b_{k,r})K_A}.            \tag{MS.7}
```

Hence, uniformly in the exterior configuration,

```math
\boxed{
\sup_v\|\mu_{A,s}(X_U\in\cdot\mid X_V=v)\|_\infty
\le {k\choose r}^{-1}
     e^{t(a_{k,r}+b_{k,r})K_A}.}                    \tag{MS.8}
```

Every marginal is a mixture of these conditional laws, so the same bound
holds for the `U`-marginal:

```math
\boxed{
\|\mu_{A,s}^{U}\|_\infty
\le {k\choose r}^{-1}
     e^{t(a_{k,r}+b_{k,r})K_A}.}                    \tag{MS.9}
```

Equivalently, define

```math
W(u)=\sum_v e^{stH_A(u,v)}.
```

Convexity of log-sum-exp and (MS.6) give directly

```math
E_R\log W(u^R)
\ge\log W(u)-t(a_{k,r}+b_{k,r})K_A,                 \tag{MS.10}
```

and arithmetic-geometric mean on the Hamming sphere proves (MS.9).
Thus (MS.8) is genuinely a conditional strengthening, while (MS.9) is the
requested vertex-marginal statement.

## 3. Uniform exponential rate for actual optimizing children

Suppose now that `A` is an actual pressure-minimizing child of order `m` at
the contracted temperature

```math
t={\beta\over\sqrt N},\qquad m\le N.
```

The optimizer contraction from Theorem 37.61 states

```math
{tK_A\over m}\le C_\beta,
\qquad C_\beta=\log2+{\beta^2\over4}.               \tag{MS.11}
```

Fix `0<theta<=1`, take any `U` with `k=|U|>=theta m`, choose
`r=floor(qk)`, and write `q_k=r/k`.  The type-class bound and
(MS.8)--(MS.11) give the fully finite estimate

```math
\begin{aligned}
-{1\over k}\log
 \sup_v\|\mu_{A,s}(X_U\in\cdot\mid X_V=v)\|_\infty
\ge{}&h(q_k)-{\log(k+1)\over k}\\
&-{C_\beta\over\theta}
 \left\{4q_k(1-q_k){k\over k-1}+2q_k\right\},       \tag{MS.12}
\end{aligned}
```

and the identical lower bound for the marginal min-entropy rate.  Here
`h(q)=-q\log q-(1-q)\log(1-q)`.

Define

```math
\eta_{\beta,\theta}^{\rm marg}
=\sup_{0<q<1/2}
 \left\{h(q)-{C_\beta\over\theta}
                   (6q-4q^2)\right\}.              \tag{MS.13}
```

This constant is strictly positive.  Indeed, with

```math
D={C_\beta\over\theta},\qquad q_*=e^{-6D},
```

the elementary inequality
`h(q)>=q log(1/q)+q(1-q)` yields

```math
\boxed{
\eta_{\beta,\theta}^{\rm marg}
\ge q_*(1-q_*)
=e^{-6C_\beta/\theta}
 (1-e^{-6C_\beta/\theta})>0.}                       \tag{MS.14}
```

Consequently, uniformly over the actual child, sector, subset
`|U|>=theta m`, and exterior condition `v`,

```math
\boxed{
\liminf_{m\to\infty}
\inf_{\substack{A\ {\rm actual},\ s,\ U:\ |U|\ge\theta m\\v\in\{\pm1\}^{V}}}
 -{1\over |U|}\log
 \|\mu_{A,s}(X_U\in\cdot\mid X_V=v)\|_\infty
\ge\eta_{\beta,\theta}^{\rm marg}.}                \tag{MS.15}
```

The same statement holds without conditioning.  In particular, with

```math
\underline\eta_{\beta,\theta}
={1\over2}e^{-6C_\beta/\theta}
 (1-e^{-6C_\beta/\theta}),                          \tag{MS.16}
```

all sufficiently large actual children satisfy

```math
\boxed{
\sup_{\substack{U:\ |U|\ge\theta m\\v,s}}
\|\mu_{A,s}(X_U\in\cdot\mid X_V=v)\|_\infty
\le e^{-\underline\eta_{\beta,\theta}|U|}
\le e^{-\theta\underline\eta_{\beta,\theta}m}.} \tag{MS.17}
```

The marginal has the same bound.

## 4. Exact and approximate common-block consequences

Let `X,X'` be independent draws from arbitrary conditional sector laws on
the same block `U`; the two exterior conditions may differ.  By (MS.17),

```math
P\{X_U=X'_U\}\le e^{-\underline\eta_{\beta,\theta}k},
\qquad
P\{X_U=\mathord\pm X'_U\}
\le2e^{-\underline\eta_{\beta,\theta}k}.            \tag{MS.18}
```

More generally, for `0<=delta<1/2`, a Hamming-ball union bound gives

```math
P\{d_H(X_U,\mathord\pm X'_U)\le\delta k\}
\le2(k+1)e^{-[\underline\eta_{\beta,\theta}-h(\delta)]k}.  \tag{MS.19}
```

Thus whenever `h(delta)<underline eta_(beta,theta)`, actual child sectors
cannot contain a positive-density frozen, common, or anticommon row block,
even after conditioning on all exterior child spins.  Equivalently, every
catalogue of `exp(xi k)` block patterns has conditional and marginal mass at
most `exp[-(underline eta_(beta,theta)-xi)k]`.

## 5. Exact scope

The theorem is optimizer-specific through (MS.11), uses only a single
global cap scalar, and is strictly lower-information than the full Gibbs
landscape.  It upgrades whole-word spread to every macroscopic coordinate
marginal and, in fact, every exact exterior conditioning.

It does **not** prove approximate independence, small total correlation,
small connected-cluster mass, or persistence along the negative bridge
path.  A posterior may remain diffuse over exponentially many patterns and
still retune their correlations coherently.  Therefore (MS.15)--(MS.19)
exclude macroscopic frozen/common blocks in the actual zero-bridge sector
law but do not close row lifetime, the balanced-product SML, or Level 6.
