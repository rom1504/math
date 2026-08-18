# Bounded-degree coherent witnesses for actual-child factor retuning

Status: **rigorous representation theorem; not yet a child-only branch
decision**.  This note starts from the exact optimal row-product shadow of
the actual optimized-child bridge law.  It proves that extensive factor
retuning cannot be hidden exclusively in exponentially detailed row
tables: at physical scale it is separated from the canonical product by
one additive Walsh polynomial of bounded row degree.  The degree depends
only on the fixed thermodynamic parameters and on the size of the claimed
retuning gap, not on the child orders.

The theorem is genuinely about the actual child law, but its coefficients
are extracted from a globally optimal product shadow.  Thus it compresses
the *representation* of branch (iii); it does not yet decide that branch
from the children without solving the product variational problem.

## 1. A uniform low-degree separation lemma

Let `U_k` be the fair law on `{+-1}^k`.  For a function `f` on the cube,
write

```math
f=\sum_{S\subseteq[k]}\widehat f(S)\chi_S,
\qquad
\Pi_{[1,d]}f=\sum_{1\le |S|\le d}\widehat f(S)\chi_S.
```

**Lemma CRW.1 (bounded-degree likelihood-ratio separator).**  Fix
`A,C,eta>0`.  There are an integer

```math
d=d(A,C,\eta)<\infty
```

and `delta=delta(A,C,eta)>0` with the following property.  Let `r,p` be
strictly positive probability laws on `{+-1}^k` such that

```math
D_2(r\Vert U_k)\le C,
\qquad
D_2(p\Vert U_k)\le C,                              \tag{CRW.1}
```

and every one-bit flip changes each of `log(dr/dU_k)` and
`log(dp/dU_k)` by at most `A/sqrt(k)`.  Put

```math
\ell=\log{dp\over dr},
\qquad
g=\Pi_{[1,d]}\ell.
```

If

```math
\chi^2(p\Vert r)\ge\eta,                            \tag{CRW.2}
```

then

```math
\boxed{\mathbb E_p g-\mathbb E_r g\ge\delta.}       \tag{CRW.3}
```

In particular, a fixed reverse-chi-square separation between two regular
weak-coordinate laws is visible in bounded Walsh degree.  The assertion is
uniform in the row length.

*Proof.*  We first record a quantitative consequence of the hypotheses.
There is a number

```math
j=j(A,C,\eta)>0                                    \tag{CRW.4}
```

such that every admissible pair obeying (CRW.2) satisfies

```math
\mathbb E_p\ell-\mathbb E_r\ell
=D(p\Vert r)+D(r\Vert p)\ge j.                     \tag{CRW.5}
```

The bit oscillation of `ell` is at most `2A/sqrt(k)`, so bounded differences
gives dimension-free exponential moments for `ell-E_U ell`.  The centering
constants are uniformly bounded as well.  Indeed, if `f=log(dq/dU)` is
either log density, then Jensen gives `E_Uf<=0`, while Hoeffding and
`E_Ue^f=1` give `E_Uf>=-A^2/8`.  Thus `E_Uell` stays in a fixed interval.
The `L^2(U)` bounds in (CRW.1) and Cauchy--Schwarz now give, for every
fixed `s>0`, a finite `K_s=K_s(A,C)` such that

```math
\mathbb E_re^{s|\ell|}+\mathbb E_pe^{s|\ell|}\le K_s. \tag{CRW.5a}
```

Put `Z=e^ell-1`.  Equation (CRW.5a) makes `Z^2` uniformly integrable under
`r`: for example, its positive-`ell` tail is bounded by an exponential
moment of order strictly larger than two, while on the negative tail
`Z^2<=1`.  Choose `L=L(A,C,eta)` so that

```math
\mathbb E_r[Z^2\mathbf1_{\{|\ell|>L\}}]\le{\eta\over2}. \tag{CRW.5b}
```

Pointwise, `(e^ell-1)ell>=0`, and on `|ell|<=L`,

```math
(e^\ell-1)\ell\ge e^{-L}(e^\ell-1)^2.              \tag{CRW.5c}
```

Indeed, for positive `ell`, `e^ell-1<=ell e^L`; for negative `ell`, the
inequality is even stronger.  Since

```math
D(p\Vert r)+D(r\Vert p)
=\mathbb E_r(e^\ell-1)\ell,
```

(CRW.2), (CRW.5b), and (CRW.5c) prove (CRW.5) with
`j=e^{-L}eta/2`.  To justify (CRW.5a) explicitly, Hoeffding gives
`sup_k E_U exp(s|ell|)<infty` for each fixed `s`; multiplication by either
`dr/dU` or `dp/dU` and Cauchy--Schwarz transfers this bound to the two
laws.  This also shows why weak-coordinate control, rather than `D_2`
alone, rules out chi-square mass hidden on a vanishing event.

Now use the Walsh Dirichlet identity.  Since the bit oscillation of `ell`
is at most `2A/sqrt(k)`,

```math
4\sum_S |S|\widehat\ell(S)^2
=\sum_{a=1}^k\mathbb E_U
  \{\ell(B)-\ell(B^{(a)})\}^2
\le4A^2.                                             \tag{CRW.6}
```

Consequently

```math
\|\ell-\widehat\ell(\varnothing)-g\|_{L^2(U)}
\le {A\over\sqrt{d+1}}.                             \tag{CRW.7}
```

The constant coefficient vanishes against `p-r`.  Also (CRW.1) gives

```math
\left\|{dp\over dU}-{dr\over dU}\right\|_{L^2(U)}
\le2e^{C/2}.                                         \tag{CRW.8}
```

Therefore

```math
\begin{aligned}
\mathbb E_pg-\mathbb E_rg
&=\mathbb E_p\ell-\mathbb E_r\ell\\
&\quad-
 \mathbb E_U\left[
 \left({dp\over dU}-{dr\over dU}\right)
 \{\ell-\widehat\ell(\varnothing)-g\}\right]\\
&\ge j-{2Ae^{C/2}\over\sqrt{d+1}}.                  \tag{CRW.9}
\end{aligned}
```

Choose `d` so that the last error is at most `j/2`, and take
`delta=j/2`. `square`

The use of both `D_2` bounds in CRW.1 is essential to the short projection
argument: it converts `L^2(U)` Walsh approximation into a directional
expectation estimate.  Merely bounding `D_2(p||r)` would not do this.

## 2. Exact score equation for the actual product shadow

Fix actual contracted-temperature minimizing children `A,D`, a relative
orientation, and a comparable split `m+n=N`.  Put

```math
u={\beta\over\sqrt N}.
```

Let `p_for(B)` be the actual forward bridge likelihood, let `p_i` be its
erased-row marginals, and write

```math
h(B)=\log p_{\rm for}(B)-\sum_{i=1}^m\log p_i(B_i).
                                                               \tag{CRW.10}
```

The canonical inverse row product is `r=tensor_i r_i`, where

```math
{dr_i\over dU_n}(b)\ \propto\ p_i(b)^{-\lambda}.
```

Let `p^*=tensor_i p_i^*` be any globally optimal row-product shadow.  The
product objective differs by a constant from

```math
\mathcal G(P)=\mathbb E_Ph
 +{1\over\lambda}\sum_iD(P_i\Vert r_i).             \tag{CRW.11}
```

Strict convexity in one factor gives the exact relative score equation

```math
\boxed{
\ell_i(b):=\log{dp_i^*\over dr_i}(b)
=-\lambda\mathbb E_{p_{-i}^*}h(b,B_{-i})-c_i.}      \tag{CRW.12}
```

Thus for every nonempty Walsh set `S subseteq[n]`,

```math
\boxed{
\widehat\ell_i(S)
=-\lambda\mathbb E_{p_{-i}^*}
 \mathbb E_{U_n}\{h(B)\chi_S(B_i)\}.}              \tag{CRW.13}
```

Equation (CRW.13) is an exact optimizer/minimality identity for the actual
children.  It uses the collision--cavity interaction `h`, not a surrogate
row law.

Both terms in `ell_i` are weak-coordinate scores.  AC.17 and the erased-row
formula give

```math
\operatorname {osc}_{b_a}\log{dp_i^*\over dU_n}\le2\lambda u,
\qquad
\operatorname {osc}_{b_a}\log{dr_i\over dU_n}\le2\lambda u,
```

and hence

```math
\operatorname {osc}_{b_a}\ell_i\le4\lambda u.       \tag{CRW.14}
```

The erased-row instance of AC.1 gives the first bound below, while the
coordinate best-response equation AC.17 followed by the same lemma gives
the second (this is AC.18 / Theorem 37.19):

```math
D_2(r_i\Vert U_n)\le\lambda^2u^2n,
\qquad
D_2(p_i^*\Vert U_n)\le\lambda^2u^2n.                \tag{CRW.15}
```

## 3. One coherent bounded-degree direction

**Theorem CRW.2 (coherent low-degree witness for extensive actual
retuning).**  Fix `beta,lambda` and a comparable-split window
`theta N<=m,n<=(1-theta)N`.  Suppose along a sequence of actual optimizing
children

```math
\mathcal J-\mathcal I^{\leftarrow}\ge\eta N.         \tag{CRW.16}
```

Then there are constants

```math
d=d(\beta,\lambda,\theta,\eta)<\infty,
\qquad c=c(\beta,\lambda,\theta,\eta)>0,             \tag{CRW.17}
```

and, for every member of the sequence, a set `I subseteq[m]` with
`|I|>=cN` such that the single additive row polynomial

```math
G(B)=\sum_{i\in I}g_i(B_i),
\qquad
g_i=\Pi_{[1,d]}\log{dp_i^*\over dr_i},               \tag{CRW.18}
```

obeys

```math
\boxed{
\mathbb E_{p^*}G-\mathbb E_rG\ge cN.}                \tag{CRW.19}
```

The witness contains at most

```math
m\sum_{a=1}^d{n\choose a}=N^{O_{\beta,\lambda,\theta,\eta}(1)} \tag{CRW.20}
```

real coefficients.  It is therefore a strict representation reduction
from the `m2^n` optimal row tables and, a fortiori, from the full
`2^{mn}` bridge landscape.

This remains a finite-information reduction after quantization.  By
(CRW.6), every coefficient vector in (CRW.18) lies in a fixed `ell_2`
ball.  If `g_i` is changed by `epsilon` in `L^2(U_n)`, then (CRW.15) and
Cauchy--Schwarz change its directional gap by at most
`2e^{C/2}epsilon`.  Taking a fixed sufficiently small `epsilon` therefore
preserves a positive fraction of (CRW.19).  Coordinatewise quantization at
mesh `epsilon/sqrt(sum_(a<=d) binom(n,a))` uses at most
`N^{O(d)}log N` bits in total.

*Proof.*  Theorem 37.27, ES.35, supplies constants `a_0,c_0>0`, depending
only on the parameters in (CRW.17), and a set `I` of at least `c_0N` rows
on which

```math
\chi^2(p_i^*\Vert r_i)\ge a_0^2.                    \tag{CRW.21}
```

At a comparable split, (CRW.14) can be written as an `A/sqrt(n)` bound
with fixed `A`, and (CRW.15) has a fixed right side.  Apply CRW.1 to every
row in `I`, with the same degree and separation constant.  Independence of
both products gives

```math
\mathbb E_{p^*}G-\mathbb E_rG
=\sum_{i\in I}(\mathbb E_{p_i^*}g_i-\mathbb E_{r_i}g_i)
\ge\delta|I|.
```

This proves (CRW.19), after decreasing `c`.  The coefficient count is
immediate. `square`

**Corollary CRW.3 (one-scalar distinguishability).**  The witness in
CRW.2 may be viewed as a single scalar observable, not a list of separately
paid row tests.  There is a constant `C_G` depending only on the fixed
parameters such that

```math
\operatorname {Var}_{p^*}G+\operatorname {Var}_rG\le C_GN.    \tag{CRW.22}
```

Consequently, thresholding `G` halfway between its two expectations
distinguishes `p^*` from `r` with error `O(1/N)` under either product law.

*Proof.*  The Dirichlet bound (CRW.6) bounds `||g_i||_(L^2(U))` uniformly.
Bonami hypercontractivity for a degree-`d` polynomial gives

```math
\|g_i\|_{L^4(U)}\le3^{d/2}\|g_i\|_{L^2(U)}.
```

The `D_2(U)` bounds in (CRW.15) and Cauchy--Schwarz therefore bound the
second moment of `g_i` under either `p_i^*` or `r_i` by a fixed constant.
The factors are independent, so variances add.  Equation (CRW.19) and
Chebyshev prove the threshold claim. `square`

Combining CRW.13 and CRW.18, every coefficient of the coherent direction
is a bounded-degree Walsh moment of the exact actual-child
collision--cavity response, averaged against the other optimal factors.
The theorem does not split the channels or pay absolute values row by row:
all retained modes are evaluated together in the one scalar response
`G`.

## 4. Operational scope and the remaining lemma

CRW.2 rules out one possible hiding place for branch (iii): extensive
retuning cannot live only in degree tending to infinity, nor does one need
the complete row tables merely to *represent a separator*.  At fixed
accuracy, one bounded-degree additive observable carries a linear
directional separation.

There are, however, two precise limitations.

1. The functions `g_i` are obtained from `p^*`.  Although their total
   description is polynomial, (CRW.13) still averages against the unknown
   optimal product factors.  The theorem is therefore not an algorithm for
   deciding branch (iii) directly from the children.
2. Nothing in weak-coordinate regularity forces `d=2`.  The sector--Gram
   state controls the quadratic interaction tangent, but CRW.2 by itself
   does not identify the physical retuning witness with that state.

The retuning part of the SML can nevertheless be stated more narrowly:

> **Bounded-degree child closure.**  Show that the coefficients in
> (CRW.13), for `|S|<=d(beta,lambda,theta,eta)`, are determined or
> uniformly controlled by an optimizer-specific child state which does not
> require the optimal row tables; or prove directly that the
> bounded-degree exponential row family generated by these modes captures
> a fixed fraction of (CRW.16).

Either result would turn the representation theorem into a genuine
low-information branch decision.  Without it, CRW.2 is a strict
information reduction but not a Level-6 recurrence.
