# Adversarial audit: effective row support dichotomy

**Object audited:**
[`../drafts/actual_child_effective_row_support_dichotomy.md`](../drafts/actual_child_effective_row_support_dichotomy.md).

**Verdict:** **PASS.**  The conditional Renyi-two comparison ES.0, the
conditional exponential path, the `D_2`-to-fourth-moment argument, every
constant in ES.5, the support count ES.18, and the dual-total-correlation
identity are correct.  The source also
states the essential scope limitation correctly: this localizes the
canonical error `D(r||q)`, not the optimal reverse projection onto all row
products.

## 1. Conditional Renyi-two comparison

The logarithmic moment lemma ES.3c has the stated constants.  If
`g=log f`, Hoeffding gives

```math
 \log E e^{a(g-Eg)}\le a^2V.
```

Using `E e^g=1` first yields `-V<=Eg<=0`.  Thus the centered bound gives
`log E f^a<=a^2V` for `a>=0`; for `a<0`, multiplication reverses which
endpoint of the interval for `Eg` is relevant and contributes exactly
`|a|V`.  In particular, the negative third moment used below is bounded by
`12V`, not merely asserted finite.

The one-bit oscillation `2lambda u` really is retained by each conditional
row law in ES.0.  In the global density IC.4, flipping a bit changes the
`-s log p` term by at most `2su` and its sole row-marginal term by at most
`2(lambda-s)u`.  Fixing coordinates preserves this pointwise ratio bound.
Marginalizing unmentioned coordinates also preserves it, because summing
the two pointwise inequalities over the same fibres bounds the ratio of the
two sums.  Hence it survives both an earlier-row prefix and arbitrary
conditioning on all other rows.  Independently, `r_i propto p_i^{-lambda}`
has the same `2lambda u` oscillation.  Both associated moment constants are
therefore at most

```math
 V_0={1\over8}n(2\lambda u)^2
 ={1\over2}\lambda^2u^2n.
```

Finally, Holder with exponents `3/2` and `3` gives

```math
 E_U{p^2\over r_i}
 <=(E_Up^3)^{2/3}(E_Ur_i^{-3})^{1/3}.
```

ES.3c contributes `(2/3)9V_0=6V_0` to the first factor and
`(1/3)(9+3)V_0=4V_0` to the second.  Thus
`D_2(p||r_i)<=10V_0=5lambda^2u^2n`, exactly ES.3a.  The comparison remains
uniform at `u=beta/sqrt(N)` and does not rely on an unstated balance
condition.

## 2. Conditional path and variance ceiling

After fixing `B_(-i)`, IC.3 gives exactly

```math
 {d\nu_t\over dr_i}
 ={e^{-tX}\over E_{r_i}e^{-tX}},
 \qquad X(R_i)=h(R_i,B_{-i}).
```

There is no residual marginal factor: it is already contained in `r_i`.
For this exponential family,

```math
 {d\over dt}D(\nu_t||\nu_0)
 =t\operatorname{Var}_{\nu_t}(X),
```

so ES.13 has the correct sign and endpoint.

A bit flip in row `i` changes `log p` by at most `2u` and the sole affected
`log p_i` by at most `2u`; hence `X` has difference constant `4u` per bit.
McDiarmid therefore has variance proxy

```math
 v=\sum_{j=1}^n(4u)^2=16u^2n
```

and tail `2exp(-2a^2/v)`.  Layer-cake integration gives

```math
 E_U(X-E_UX)^4
 <=8\int_0^\infty a^3e^{-2a^2/v}\,da=v^2.
```

With `f=dnu_t/dU`, Cauchy--Schwarz yields

```math
 \operatorname{Var}_{\nu_t}(X)
 <=(E_Uf^2)^{1/2}(E_UZ^4)^{1/2}
 <=e^{C/2}\,16u^2n,
```

because IC.2 gives `D_2(nu_t||U)<=C=lambda^2u^2n`.  Integrating
`t Var_(nu_t)(X)` contributes the factor `s^2/2`; after division by `s^2`
this is exactly

```math
 e_i(s)<=8e^{C/2}u^2n=L.
```

The zero endpoint is `Var_(r_i)(X)/2`, so the same bound survives when the
witness parameter tends to zero.  Substitution `u=beta/sqrt(N)` gives ES.6
without an omitted balance assumption.

## 3. Support count

Let `a=eta/(2lambda^2)` and `k=K_s(a)`.  Rows above threshold contribute at
most `kL`; all other rows contribute at most `a`, and `m<=N`.  Hence

```math
 {\eta\over\lambda^2}N
 <=\sum_i e_i(s)
 <=kL+{\eta\over2\lambda^2}N,
```

which gives ES.18 with exactly the displayed constant.  If the resulting
lower bound exceeds the available row count, the hypothesis `J>=eta N` is
itself impossible; this is not a defect.  For physical fixed parameters,
`L=O(1)`, so the conclusion really is a positive-density row set.

## 4. Entropy-production and DTC identities

For the centered cumulant `K`, direct evaluation gives

```math
 D(q_s||r)=sK'(s)-K(s).
```

Since `(K(s)/s)'=D(q_s||r)/s^2` and `K(s)/s->0`, ES.21 follows with the
outer factor `lambda` in the correct place.

For row-product `r`, expansion into entropies gives

```math
 \sum_iE_qD(q(R_i|B_{-i})||r_i)-D(q||r)
 =H(q)-\sum_iH_q(R_i|B_{-i}),
```

which is dual total correlation and is nonnegative.  Therefore the exact
slack of IC.4 is ES.23.  The ordinary total-correlation split ES.25 and the
ordered chain rule ES.26 are also in the correct KL direction; integrating
them gives ES.27--ES.29.

## 5. Scope

The result is uniform for arbitrary finite child priors, so it is valid for
actual optimizing children but does not use their minimality.  It proves
that a **linear canonical mismatch** cannot live on `o(N)` rows.  It does
not prove the best row-product distance is linear, because

```math
 \inf_{a=\otimes_i a_i}D(a||q)\le D(r||q)=\mathcal J.
```

Likewise, ES.21 and the integrated forms are exact reformulations of the
canonical cumulant, not strict reductions by themselves.  The theorem file
states both limitations explicitly, so no scope correction is needed.

## 6. Canonical-to-best-product retuning alternative

Theorem ES.3 also passes.  The comparison `D_2(p_i^*||r_i)<=5C` implies

```math
a_i^2=\chi^2(p_i^*\Vert r_i)\le e^{5C}-1=X.
```

When the product factors are changed one at a time, all other factors only
enter through the average `f_i(b)=E[L(B)|B_i=b]`.  A bridge-bit flip changes
`L` by at most `2u`, so averaging against any intermediate product retains
the same `2u` bound.  Its fair-cube bounded-difference parameter is therefore
`4u^2n`.  Combining the fourth-moment estimate with
`D_2(r_i||U_n)<=C` gives exactly

```math
\operatorname {Var}_{r_i}(f_i)
\le4e^{C/2}u^2n.
```

Chi-square Cauchy--Schwarz then proves ES.38 with coefficient
`2e^{C/4}u sqrt(n)`.

The potentially delicate absolute entropy-difference step admits no
counterexample: it is an immediate exact change-of-reference identity.  If
`g_i=log(dr_i/dU_n)`, then

```math
D(p_i^*\Vert U_n)-D(r_i\Vert U_n)
=D(p_i^*\Vert r_i)+(E_{p_i^*}-E_{r_i})g_i.
```

The first term is at most
`D_2(p_i^*||r_i)=log(1+a_i^2)<=a_i^2`.  Since `g_i` has bit oscillation
`2lambda u`, the same fourth-moment argument gives

```math
\operatorname {Var}_{r_i}(g_i)
\le4e^{C/2}\lambda^2u^2n,
```

and chi-square Cauchy--Schwarz bounds the absolute second term by
`2e^{C/4}lambda u sqrt(n) a_i`.  This proves ES.40 in both signs; no
unspoken ordering of the two entropies is used.

There are two copies of the linear expectation-change coefficient after
division by `lambda`: one from energy and one from entropy.  Together they
give `4e^{C/4}u sqrt(n) a_i`.  The remaining entropy cost obeys
`a_i^2<=sqrt(X)a_i`, yielding exactly

```math
K=4e^{C/4}u\sqrt n+{\sqrt X\over\lambda}.
```

The Gibbs variational identity gives
`F(r)-F(p^*)=(J-I^leftarrow)/lambda>=0`, because `r` is itself a row
product.  Thus ES.33--ES.34 follow.  Finally, if `k` rows exceed
`eta/(2lambda K)`, the uniform ceiling `a_i<=sqrt(X)` gives

```math
\sum_i a_i\le k\sqrt X+{\eta\over2\lambda K}N,
```

and comparison with ES.34 proves the constant in ES.35.  As usual, if the
displayed lower bound exceeds the number of available rows, the antecedent
is impossible.  The degenerate `u=0` case should be read separately (all
three gaps vanish); at the stated positive physical temperature the
displayed denominators are nonzero.

This is a strict structural localization of a large
`J-I^leftarrow`: a better product shadow cannot repair the canonical one by
retuning only `o(N)` factors.  It still does not establish that this gap is
large or small for optimizing children, a limitation the theorem states.
