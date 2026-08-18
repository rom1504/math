# Effective row support in the actual child interaction path

Status: **rigorous actual-law structural theorem**.  This note sharpens the
conditional-entropic influence alternative in IC.4.  It proves that a
linear canonical interaction cumulant cannot be hidden in a vanishing
number of exceptional bridge rows, even if the witnessing hybrid parameter
approaches zero.

The theorem applies to the negative-disorder law induced by arbitrary
finite children; in particular it applies to every pair of
contracted-temperature optimizing children.  It does not decide which side
of the dichotomy optimizing children occupy asymptotically.

## 1. Setup

Use the notation of IC.1--IC.4.  There are `m` bridge rows of length `n`,
`N=m+n`, the channel amplitude is `u`,

```math
r=\bigotimes_{i=1}^m r_i,
\qquad {dq_s\over dr}\propto e^{-s h},
\qquad 0\le s\le\lambda,
```

and

```math
\mathcal J=D(r\Vert q_\lambda).
```

For `s>0` define the scaled conditional influence of row `i` by

```math
e_i(s)={1\over s^2}E_{q_s(B_{-i})}
D\bigl(q_s(R_i\mid B_{-i})\Vert r_i\bigr).          \tag{ES.1}
```

At zero use the continuous value

```math
e_i(0)={1\over2}E_{r(B_{-i})}
\operatorname {Var}_{r_i}\bigl(h(R_i,B_{-i})\bigr).\tag{ES.2}
```

Thus the optimized conditional-entropic influence in IC.4 is

```math
\mathcal E_s=\sum_{i=1}^m e_i(s).                  \tag{ES.3}
```

For each row, the chain rule gives the exact operational split

```math
\boxed{
e_i(s)={1\over s^2}\left\{
D(q_{s,i}\Vert r_i)
+I_{q_s}(R_i;B_{-i})\right\}.}                    \tag{ES.3-}
```

Thus `e_i` is precisely canonical marginal drift plus irreducible
row-versus-rest information, at the quadratic scale of the hybrid path.

Before localizing this KL response, one can close the conditional
Renyi-complexity question relative to the *natural* canonical factor, not
only relative to a fair row.

**Theorem ES.0 (tight conditional Renyi-two relative to the canonical
row).**  For every `s in [0,lambda]`, every row order, every prefix value,
and also after conditioning on all other rows,

```math
\boxed{
D_2\bigl(q_s(R_i\mid\mathcal C)\Vert r_i\bigr)
\le5\lambda^2u^2n.}                                \tag{ES.3a}
```

Here `mathcal C` can be a row prefix or an arbitrary fixed configuration of
the other rows.  At `u=beta/sqrt(N)` this is uniformly
`5lambda^2beta^2n/N=O(1)`.  Hence no mass of conditional components can
escape to growing Renyi-two complexity anywhere on the actual hybrid path.

*Proof.*  We use the following finite-cube moment lemma.  If `f=dP/dU_n`,
`E_U f=1`, and every bit flip changes `log f` by at most `c_j`, put

```math
V={1\over8}\sum_jc_j^2.                              \tag{ES.3b}
```

Writing `g=log f`, bounded differences and normalization give

```math
-V\le E_Ug\le0,
\qquad
\log E_Uf^a\le
\begin{cases}
a^2V,&a\ge0,\\
(a^2+|a|)V,&a<0.
\end{cases}                                         \tag{ES.3c}
```

Indeed Hoeffding's lemma bounds
`log E exp{a(g-Eg)}` by `a^2V`; take `a=1` and use `E f=1` to get the lower
bound on `Eg`, while Jensen gives the upper bound.

Let `p` be any conditional row law displayed in (ES.3a), with densities
still denoted `p,r_i` relative to `U_n`.  The bit-oscillation proof of IC.2
gives `c_j=2lambda u` for `log p`; the same bound for `log r_i` follows
directly from (CC.5).  Thus both moment constants in (ES.3b) are at most

```math
V_0={1\over2}\lambda^2u^2n.                         \tag{ES.3d}
```

Holder with conjugate exponents `3/2` and `3`, followed by (ES.3c), gives

```math
\begin{aligned}
\exp D_2(p\Vert r_i)
 &=E_{U_n}{p^2\over r_i}\\
 &\le(E_{U_n}p^3)^{2/3}(E_{U_n}r_i^{-3})^{1/3}\\
 &\le\exp\{6V_0+4V_0\}.
\end{aligned}                                       \tag{ES.3e}
```

Since `10V_0=5lambda^2u^2n`, this proves (ES.3a). `square`

The identical argument applies to every factor `p_i^*` of a globally best
row-product shadow from AC.3, because its mean-field equation AC.17 also
has one-bit log oscillation at most `2lambda u`.  Therefore

```math
\boxed{D_2(p_i^*\Vert r_i)\le5\lambda^2u^2n.}       \tag{ES.3f}
```

Thus the canonical factors, every hybrid conditional factor, and every
factor of the optimal row-product variational competitor all lie in one
order-one Renyi neighborhood at physical scaling.  The unresolved resource
is the accumulation and compatibility of those bounded row changes, not an
escaping component.

## 2. A uniform per-row ceiling

**Theorem ES.1 (no sparse-row concentration).**  Put

```math
C=\lambda^2u^2n,
\qquad
L=8e^{C/2}u^2n.                                    \tag{ES.4}
```

For every `0<=s<=lambda` and every row `i`,

```math
\boxed{0\le e_i(s)\le L.}                         \tag{ES.5}
```

At physical amplitude `u=beta/sqrt(N)`, this is the order-independent
bound

```math
e_i(s)\le
8\beta^2{n\over N}
\exp\left\{{\lambda^2\beta^2n\over2N}\right\}.
                                                               \tag{ES.6}
```

*Proof.*  Fix `i` and `B_{-i}`, and write

```math
X(R_i)=h(R_i,B_{-i}),
\qquad
\nu_t=q_t(R_i\mid B_{-i}).                          \tag{ES.7}
```

The conditional hybrid identity gives

```math
{d\nu_t\over dr_i}\propto e^{-tX}.                 \tag{ES.8}
```

Changing one bridge bit changes `log p` and `log p_i` by at most `2u`
each, and hence changes `X=log p-sum_j log p_j` by at most `4u`.  If `U_n`
is uniform on the row cube and `Z=X-E_{U_n}X`, bounded differences gives

```math
P_{U_n}(|Z|\ge a)\le2\exp\{-2a^2/v\},
\qquad v=16u^2n.                                   \tag{ES.9}
```

Integrating the tail yields

```math
E_{U_n}Z^4\le v^2.                                 \tag{ES.10}
```

The uniform conditional regularity theorem IC.2 gives

```math
D_2(\nu_t\Vert U_n)\le C                           \tag{ES.11}
```

for every `t` and every fixed outside row configuration.  Cauchy--Schwarz
therefore implies

```math
\operatorname {Var}_{\nu_t}(X)
\le E_{\nu_t}Z^2
\le e^{C/2}(E_{U_n}Z^4)^{1/2}
\le16e^{C/2}u^2n.                                  \tag{ES.12}
```

For the exponential family (ES.8), direct differentiation gives

```math
D(\nu_s\Vert\nu_0)
=\int_0^s t\operatorname {Var}_{\nu_t}(X)\,dt.     \tag{ES.13}
```

Since `nu_0=r_i`, equations (ES.12)--(ES.13) give the pointwise bound

```math
D\bigl(q_s(R_i\mid B_{-i})\Vert r_i\bigr)
\le8e^{C/2}s^2u^2n.                                \tag{ES.14}
```

Averaging over `B_{-i}`, dividing by `s^2`, and taking the zero limit prove
(ES.5). `square`

## 3. The actual-child dichotomy has extensive effective support

Define the threshold support at scale `a>0` by

```math
K_s(a)=\#\{i:e_i(s)\ge a\}.                         \tag{ES.15}
```

**Corollary ES.2 (extensive-row alternative).**  If for some fixed
`eta>0`

```math
\mathcal J\ge\eta N,                               \tag{ES.16}
```

then there is `s_N in [0,lambda]` for which

```math
\sum_i e_i(s_N)\ge {\eta\over\lambda^2}N           \tag{ES.17}
```

and

```math
\boxed{
K_{s_N}\left({\eta\over2\lambda^2}\right)
\ge {\eta\over2\lambda^2L}\,N.}                 \tag{ES.18}
```

Thus a linear canonical row-product error forces a fixed positive fraction
of bridge rows to have nonzero order-one **scaled conditional KL response**.
The statement remains nondegenerate if `s_N` tends to zero, because the
normalization in (ES.1) converges to the conditional variance in (ES.2).

By (ES.3-), at least half of the rows counted in (ES.18) have one of

```math
{D(q_{s_N,i}\Vert r_i)\over s_N^2}
\ge {\eta\over4\lambda^2},
\qquad
{I_{q_{s_N}}(R_i;B_{-i})\over s_N^2}
\ge {\eta\over4\lambda^2}.                        \tag{ES.18a}
```

Consequently the extensive support is forced into an explicit
marginal-retuning or irreducible-row-information channel; a small hybrid
parameter cannot conceal it.

*Proof.*  IC.4 gives (ES.17).  If `k=K_s(eta/(2lambda^2))`, then (ES.5) and
`m<=N` give

```math
\sum_i e_i(s)
\le kL+(m-k){\eta\over2\lambda^2}
\le kL+{\eta\over2\lambda^2}N.                    \tag{ES.19}
```

Compare (ES.17) and (ES.19). `square`

Equations (ES.5) and (ES.18) make the phrase **effective coordinate
support** quantitative.  They exclude a sparse exceptional-row explanation
of a failed canonical product certificate.  They do not say that the raw
conditional KL at the witness is linear: it is `s_N^2 sum_i e_i(s_N)`, and
that distinction is necessary when `s_N` tends to zero.

## 4. Exact entropy-production identity and its limitation

There is also an exact identity which clarifies what IC.4 does and does not
reduce.  Let

```math
K(s)=\log E_r\exp\{-s(h-E_rh)\}.
```

Then

```math
D(q_s\Vert r)=sK'(s)-K(s),                          \tag{ES.20}
```

so

```math
\boxed{
\mathcal J
=\lambda\int_0^\lambda {D(q_s\Vert r)\over s^2}\,ds.}
                                                               \tag{ES.21}
```

Moreover, if `\mathrm{DTC}(q_s)` denotes dual total correlation across
bridge rows, then

```math
\sum_iE_{q_s}D(q_s(R_i\mid B_{-i})\Vert r_i)
=D(q_s\Vert r)+\mathrm{DTC}(q_s).                  \tag{ES.22}
```

Consequently IC.4 has the exact slack

```math
\lambda\int_0^\lambda
 {\mathrm{DTC}(q_s)\over s^2}\,ds.                 \tag{ES.23}
```

*Proof.*  Equation (ES.20) is differentiation of the exponential family,
and integrating `(K(s)/s)'=D(q_s||r)/s^2` proves (ES.21).  For (ES.22),
expand both sides into Shannon entropies and use that `r` is a row product:

```math
\sum_i\{ -H_{q_s}(R_i\mid B_{-i})-E_{q_s}\log r_i(R_i)\}
-D(q_s\Vert r)
=H(q_s)-\sum_iH_{q_s}(R_i\mid B_{-i}).             \tag{ES.24}
```

The last expression is dual total correlation. `square`

The same forward divergence also separates ordinary row total correlation
from drift of the one-row channels.  If `q_{s,i}` is the row-`i` marginal,
then

```math
D(q_s\Vert r)
=\operatorname {TC}(q_s)+\sum_iD(q_{s,i}\Vert r_i),              \tag{ES.25}
```

and, for every fixed ordering of the rows,

```math
D(q_s\Vert r)
=\sum_iE_{q_s}D\bigl(q_s(R_i\mid R_{<i})\Vert r_i\bigr).         \tag{ES.26}
```

Combining these identities with (ES.21) gives the exact filtration and
row-correlation representations

```math
\mathcal J
=\lambda\int_0^\lambda {1\over s^2}
 \sum_iE_{q_s}D\bigl(q_s(R_i\mid R_{<i})\Vert r_i\bigr)\,ds,
                                                               \tag{ES.27}
```

```math
\mathcal J
=\lambda\int_0^\lambda {1\over s^2}
 \left\{\operatorname {TC}(q_s)
       +\sum_iD(q_{s,i}\Vert r_i)\right\}\,ds.                  \tag{ES.28}
```

In particular, `J>=eta N` forces at least one of

```math
\lambda\int_0^\lambda {\operatorname {TC}(q_s)\over s^2}\,ds
\ge{\eta\over2}N,
\qquad
\lambda\int_0^\lambda {\sum_iD(q_{s,i}\Vert r_i)\over s^2}\,ds
\ge{\eta\over2}N.                                      \tag{ES.29}
```

The first alternative is irreducible cross-row dependence under the
standard forward information projection; the second says that the
canonical erased-row factors drift collectively and leaves open a better
product shadow.  Equations (ES.25)--(ES.29) are KL chain rules, but their
combination with the exact actual-child interpolation identifies precisely
which row resource must carry a failed canonical certificate.  Every
conditional row appearing in (ES.27) nevertheless has the order-one
Renyi-two bound IC.10.

This identity is an important scope warning.  The optimized conditional
influence is a useful *localization* of a linear obstruction, but by itself
it is not a strict reduction: it contains the forward hybrid divergence
whose integral is exactly the canonical error.  Progress now requires an
optimizer-specific child theorem controlling the individual row responses
in (ES.1), or a coarser high-transport statistic which certifies that a
positive fraction persist.

## 5. Canonical error, optimal dependence, or extensive row retuning

The same local regularity also prevents a better product shadow from
repairing the canonical product through only a few exceptional factors.
Write

```math
\mathcal F(P)=E_PL+{1\over\lambda}D(P\Vert U_B),     \tag{ES.30}
```

let `p^*=tensor_i p_i^*` minimize this functional over row products, and
put

```math
\mathcal I^{\leftarrow}
=\inf_{P\ {\rm row\ product}}D(P\Vert q_\lambda).
```

For each row define

```math
\chi_i^2=\chi^2(p_i^*\Vert r_i),
\qquad a_i=\sqrt{\chi_i^2}.                          \tag{ES.31}
```

**Theorem ES.3 (extensive retuning alternative).**  With `C` as in (ES.4),
put

```math
X=e^{5C}-1,
\qquad
K=4e^{C/4}u\sqrt n+{\sqrt X\over\lambda}.          \tag{ES.32}
```

Then

```math
\boxed{
0\le {\mathcal J-\mathcal I^{\leftarrow}\over\lambda}
=\mathcal F(r)-\mathcal F(p^*)
\le K\sum_{i=1}^m a_i.}                            \tag{ES.33}
```

In particular, if `J-I^leftarrow>=eta N`, then

```math
\sum_i a_i\ge{\eta\over\lambda K}N,               \tag{ES.34}
```

and a positive density of rows is separated from the canonical factor:

```math
\boxed{
\#\left\{i:a_i\ge{\eta\over2\lambda K}\right\}
\ge {\eta\over2\lambda K\sqrt X}\,N.}            \tag{ES.35}
```

When `u=beta/sqrt(N)` and the split is comparable, `C,X,K` are fixed
constants depending only on `beta,lambda` and the aspect ratio.
The case `u=0` is the trivial law `q=r=U` and is read separately; the
displayed positive-gap alternatives are then vacuous.

*Proof.*  Equation (ES.3f) gives `a_i^2<=X`.  Telescope from `r` to `p^*`
one product factor at a time.  With all other factors fixed at the
corresponding intermediate product, set

```math
f_i(b)=E[L(B)\mid B_i=b].                            \tag{ES.36}
```

A bit flip changes `f_i` by at most `2u`.  Under the fair row law its
bounded-difference parameter is `4u^2n`, and the fourth-moment argument of
ES.1, together with `D_2(r_i||U_n)<=C`, gives

```math
\operatorname {Var}_{r_i}(f_i)
\le4e^{C/2}u^2n.                                    \tag{ES.37}
```

The chi-square Cauchy--Schwarz inequality therefore yields

```math
|E_{p_i^*}f_i-E_{r_i}f_i|
\le2e^{C/4}u\sqrt n\,a_i.                           \tag{ES.38}
```

For the entropy term put `g_i=log(dr_i/dU_n)`.  Its bit oscillation is at
most `2lambda u`, so the same argument gives

```math
\operatorname {Var}_{r_i}(g_i)
\le4e^{C/2}\lambda^2u^2n.                           \tag{ES.39}
```

The exact change of reference identity and
`D(p_i^*||r_i)<=log(1+chi_i^2)<=a_i^2` imply

```math
\begin{aligned}
|D(p_i^*\Vert U_n)-D(r_i\Vert U_n)|
&\le a_i^2
 +2e^{C/4}\lambda u\sqrt n\,a_i.                  \tag{ES.40}
\end{aligned}
```

After division by `lambda`, (ES.38)--(ES.40), telescoping, and
`a_i^2<=sqrt(X)a_i` prove the upper bound in (ES.33).  Its identity and
nonnegativity follow from AC.3 and the Gibbs variational identity:

```math
\mathcal F(r)=V_\lambda+{\mathcal J\over\lambda},
\qquad
\mathcal F(p^*)=V_\lambda
 +{\mathcal I^{\leftarrow}\over\lambda}.           \tag{ES.41}
```

Equation (ES.34) follows immediately.  If `k` rows have
`a_i>=eta/(2lambda K)`, then

```math
\sum_i a_i
\le k\sqrt X+{\eta\over2\lambda K}N;
```

comparison with (ES.34) gives (ES.35). `square`

Theorems ES.0, ES.2, and ES.3 give a finite resource classification for the
actual child law.  Asymptotically, after passing to any subsequence on which
the relevant normalized gaps have definite zero/positive behavior, the
alternatives are:

1. `J=o(N)`: the explicit one-child iid product is asymptotically accurate;
2. `I^leftarrow=Omega(N)` with `J-I^leftarrow=o(N)`: irreducible reverse
   product dependence carries the canonical obstruction; or
3. `J-I^leftarrow=Omega(N)`: a positive density of individually regular
   row factors must retune by order one.

Mixtures of the last two alternatives are allowed, and different
subsequences may occupy different branches.  No branch can be
attributed to non-tight conditional Renyi components or to a vanishing set
of bad rows.  The theorem does not determine which branch contracted-
temperature minimizing children occupy.

Equivalently, no regularity of normalized gaps is needed for the following
pointwise version.  Whenever `J>=eta N`, either

```math
\mathcal I^{\leftarrow}\ge{\eta\over2}N             \tag{ES.42}
```

or `J-I^leftarrow>=eta N/2`, in which case (ES.34)--(ES.35) hold with
`eta/2`.  Thus every infinite linear-canonical-error subsequence contains
an infinite irreducible-dependence subsequence or an infinite
positive-density-retuning subsequence.

## 6. Smallest missing lemma

For contracted-temperature minimizing children, prove one of:

1. the integrated scaled row response in (ES.21) is `o(N)` from a child
   statistic strictly smaller than the full external-field landscape; or
2. `I^leftarrow>=cN` from a high-transport aggregate child observable; or
3. the extensive-retuning branch in ES.35 occurs, from an optimizer-specific
   statistic which identifies the common direction of the regular row
   shifts without solving the full product variational problem.

The last two statements separate irreducible dependence from collective
factor retuning.  ES.1 and ES.3 show that neither a vanishing set of bad rows
nor a near-zero hybrid parameter can evade this formulation.
