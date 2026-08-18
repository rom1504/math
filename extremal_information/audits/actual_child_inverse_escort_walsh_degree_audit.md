# Walsh degree under the actual inverse bridge escort

**Status.**  Rigorous degree audit.  The complete bridge likelihood has an
exact two-replica formula for its Walsh-level masses.  At strong fixed
temperature, every actual rank-one child law has exponentially large
ordinary `L^2(U)` mass above every degree `o(N)`.  Thus an argument which
first approximates the likelihood in ordinary `L^2`, or controls it through
positive replicas and Cauchy--Schwarz, is forced to degree `Omega(N)`.

There is nevertheless a universal degree-`O(N)` approximation theorem in
the **correct inverse-escort norm**, including simultaneous approximation of
all exact edge-cavity responses.  Its proof explicitly pays an extensive
negative-moment cost.  The resulting linear-degree object is not a strict
compression: on rank-one words it already contains the complete projective
latent law.

The ordinary-`L^2` lower bound is not an inverse-escort lower bound.  An
antipodal one-atom channel has the same extensive high-degree `L^2` tail,
while its inverse escort makes the log likelihood constant up to `O(1)` and
all cavities functions of one linear statistic.  Hence the unresolved
question is genuinely a **negative-tail/conditional response** question,
not a Walsh-tail question under the fair bridge law.

## 1. Setup

Let

```math
d=mn,qquad N=m+n,qquad \beta>0,qquad t={\beta\over\sqrt N},\qquad
\rho=\tanh t,qquad U=\operatorname {Unif}\{\pm1\}^{d}.
                                                               \tag{EW.1}
```

Let `mu` be an antipodally symmetric law on rank-one bridge words

```math
Q=(Q_{ij})=(\tau X_iY_j)\in\{\pm1\}^{m\times n}.
```

The normalized bridge likelihood and its inverse escort are

```math
P(B)=E_\mu\prod_{e=1}^{d}(1+\rho B_eQ_e),
\qquad {dq_\lambda\over dU}={P^{-\lambda}\over Z_\lambda},
\qquad Z_a=E_UP^{-a}.                              \tag{EW.2}
```

Thus `E_UP=1`.  The additive constant between `log P` and the bridge
pressure is immaterial here.  For an actual pair of contracted-temperature
children, Theorem 37.52 supplies precisely (EW.2) with

```math
|\operatorname {supp}\mu|\le 2^{N-1}.             \tag{EW.3}
```

For `S subseteq[d]`, put

```math
m(S)=E_\mu\chi_S(Q).
```

Theorem 37.71 gives

```math
\widehat P(S)=\rho^{|S|}m(S)                       \tag{EW.4}
```

at even `|S|` and zero at odd `|S|`.  Write `P_(<=K)` for the Walsh
projection of `P` onto degrees at most `K`, and

```math
T_K=\|P-P_{\le K}\|_{L^2(U)}^2.                   \tag{EW.5}
```

## 2. The exact two-replica level enumerator

### Theorem EW.1 (Walsh mass is a two-replica overlap polynomial)

Define

```math
W_k=\sum_{|S|=k}|\widehat P(S)|^2.
```

For every scalar `z`, one has the polynomial identity

```math
\boxed{
\sum_{k=0}^{d}W_kz^k
=E_{Q,Q'\sim\mu}\prod_{e=1}^{d}
       (1+z\rho^2Q_eQ'_e).}                       \tag{EW.6}
```

If `C=<Q,Q'>`, then

```math
\boxed{
\sum_kW_kz^k
=E_{Q,Q'}(1+z\rho^2)^{(d+C)/2}
          (1-z\rho^2)^{(d-C)/2}.}                 \tag{EW.7}
```

For rank-one replicas,

```math
C=(\tau\tau')\langle X,X'\rangle\langle Y,Y'\rangle.          \tag{EW.8}
```

Thus the complete ordinary Walsh-level profile is exactly a two-replica
child-overlap object; no higher replicas are needed for this particular
`L^2` question.

*Proof.*  By (EW.4),

```math
\begin{aligned}
\sum_kW_kz^k
&=\sum_{S\subseteq[d]}(z\rho^2)^{|S|}
  E_{Q,Q'}\chi_S(QQ')\\
&=E_{Q,Q'}\prod_e(1+z\rho^2Q_eQ'_e).
\end{aligned}
```

There are `(d+C)/2` coordinates on which `Q_eQ'_e=1` and `(d-C)/2`
on which it is `-1`, proving (EW.7).  Equation (EW.8) is the rank-one
factorization.  `square`

In particular,

```math
T_K=\sum_{k>K}W_k
\le\sum_{k>K}{d\choose k}\rho^{2k}.               \tag{EW.9}
```

This upper bound will give the linear-degree recovery theorem.  The next
section shows that it cannot be replaced by a useful `o(N)` ordinary-`L^2`
bound at the temperatures relevant to the zero-temperature program.

## 3. An actual-law `Omega(N)`-degree obstruction in ordinary `L^2`

### Theorem EW.2 (support collision forces an extensive high-degree tail)

Assume only the actual-law support bound (EW.3).  If

```math
\liminf_N{d\over N^2}\ge\gamma_0>0,
\qquad \gamma_0\beta^2>\log2,                     \tag{EW.10}
```

then, for every sequence `K_N=o(N)`, there is `c>0` such that

```math
\boxed{
\inf_{\deg f\le K_N}\|P-f\|_{L^2(U)}^2
=T_{K_N}\ge e^{cN}}                               \tag{EW.11}
```

for all sufficiently large `N`.  More precisely, the exponential rate in
(EW.11) may be taken arbitrarily below
`gamma_0 beta^2-log 2`.

At a balanced split, `gamma_0=1/4`, so the threshold is

```math
\boxed{\beta>\sqrt{4\log2}=1.665109\ldots .}       \tag{EW.12}
```

This applies uniformly to the bridge laws of all actual minimizing-child
pairs and both orientations.

*Proof.*  Let

```math
k_Q(B)=\prod_e(1+\rho B_eQ_e).
```

For every pair of latent words,

```math
\langle k_Q,k_{Q'}\rangle_U
=\prod_e(1+\rho^2Q_eQ'_e)>0.                      \tag{EW.13}
```

Keeping only diagonal replica pairs and using collision at least inverse
support size gives

```math
\begin{aligned}
\|P\|_2^2
&=E_{Q,Q'}\prod_e(1+\rho^2Q_eQ'_e)\\
&\ge\Pr(Q=Q')(1+\rho^2)^d\\
&\ge2^{-(N-1)}(1+\rho^2)^d\\
&\ge\exp\{(\gamma_0\beta^2-\log2-o(1))N\}.        \tag{EW.14}
\end{aligned}
```

On the other hand, `|m(S)|<=1`, and hence

```math
\|P_{\le K}\|_2^2
\le\sum_{j\le K}{d\choose j}\rho^{2j}.           \tag{EW.15}
```

Put `a_N=d rho^2=Theta(N)`.  For `K=o(N)`,

```math
\log\sum_{j\le K}{d\choose j}\rho^{2j}
\le O(\log K)+K\log{e a_N\over K}=o(N).           \tag{EW.16}
```

Orthogonality says that `P_(<=K)` is the best degree-`K` approximant.
Subtracting (EW.15) from (EW.14) proves (EW.11).  `square`

The obstruction uses the **actual** rank-one support size, not a conference
surrogate.  It also pinpoints the collision with Theorems 37.39--37.44:
those results compress each candidate row-product *factor* to fixed row
degree.  They do not approximate the complete child likelihood in global
bridge Walsh degree.  The latter necessarily has extensive ordinary-`L^2`
degree at strong temperature.

## 4. Negative moments and why they cost an exponential factor

Let `X=log P`.  Conditional insertion gives

```math
P(B)=P_{-e}(B_{-e})\{1+\rho B_er_e(B_{-e})\},
\qquad |r_e|\le1.                                 \tag{EW.17}
```

Consequently the oscillation of `X` when one bridge bit is flipped is at
most

```math
\log{1+\rho\over1-\rho}=2t.                      \tag{EW.18}
```

### Lemma EW.3 (uniform negative-moment envelope)

For every `a>=0`,

```math
\boxed{
1\le Z_a=E_UP^{-a}
\le\exp\left\{{a(a+1)\over2}dt^2\right\}.}       \tag{EW.19}
```

*Proof.*  Bounded differences on the fair cube gives

```math
\log E_Ue^{s(X-E_UX)}\le{s^2dt^2\over2}.          \tag{EW.20}
```

Since `E_Ue^X=E_UP=1`, (EW.20) at `s=1` gives
`E_UX>=-dt^2/2`, while Jensen gives `E_UX<=0`.  Apply (EW.20) at `s=-a`
to obtain the upper bound.  The lower bound is Jensen.  `square`

The scale `Theta(N)` in (EW.19) is genuine even inside the rank-one
channel.  If `mu` is uniform on one antipodal pair `{q,-q}`, then, with
`S=<B,q>`,

```math
P(B)=(1-\rho^2)^{d/2}\cosh(tS),                   \tag{EW.21}
```

and a local central-binomial estimate gives, for fixed `a>0` and comparable
splits,

```math
\boxed{
\log E_UP^{-a}
={a\over2}d\rho^2+O_{a,\beta}(\log N)=\Theta(N).} \tag{EW.22}
```

Indeed, the first factor in (EW.21) contributes
`a d rho^2/2+O(d rho^4)`, and
`E cosh(tS)^(-a)=Theta((t sqrt d)^(-1))=Theta(N^(-1/2))`.

Thus Cauchy--Schwarz transfer from `U` to the inverse escort must pay an
`exp{Theta(N)}` negative-moment factor.  This is why polynomial ordinary
`L^2` accuracy is useless and why the Poisson Walsh tail must be pushed past
linear degree in the universal argument below.

## 5. Universal linear-degree recovery in the inverse-escort norm

The following is a positive theorem, but not a compression theorem.

### Theorem EW.4 (degree `O(N)` recovers pressure and all cavities)

Fix `beta,lambda` and comparable splits.  There is a finite constant
`c=c(beta,lambda,gamma_0)>0` such that, for

```math
K_N=\lceil cN\rceil,                              \tag{EW.23}
```

one can construct from the degree-`K_N` Walsh coefficients of `P`:

1. a positive clipped decoder `widetilde P_K`; and
2. decoders `widetilde r_(e,K) in[-1,1]` for every edge cavity,

such that, uniformly over every rank-one latent law and hence every actual
minimizing-child law,

```math
\boxed{
E_{q_\lambda}|\log P-\log\widetilde P_K|
\le e^{-c_1N},}                                   \tag{EW.24}
```

and

```math
\boxed{
\sum_{e=1}^{d}E_{q_\lambda}
 |r_e-\widetilde r_{e,K}|\le e^{-c_1N}}           \tag{EW.25}
```

for some `c_1>0`.  The same conclusion holds for the summed cavity fields
`atanh(rho r_e)`.

*Proof.*  Let `R_K=P-P_(<=K)` and retain the exact tail `T_K=||R_K||_2^2`.
By Cauchy--Schwarz and (EW.19),

```math
\begin{aligned}
A_K
&:=E_{q_\lambda}{|R_K|\over P}\\
&={E_U|R_K|P^{-(\lambda+1)}\over Z_\lambda}\\
&\le {\sqrt{T_KZ_{2\lambda+2}}\over Z_\lambda}\\
&\le\sqrt{T_K}\,
 \exp\left\{{(2\lambda+2)(2\lambda+3)\over4}dt^2\right\}.
                                                               \tag{EW.26}
\end{aligned}
```

Let `a_N=d rho^2`.  If `K+1>=2a_N`, then (EW.9), the factorial bound, and
a geometric tail give

```math
T_K\le2\left({e a_N\over K+1}\right)^{K+1}.       \tag{EW.27}
```

Since `a_N=Theta(N)` and `dt^2=Theta(N)`, a sufficiently large fixed `c`
in (EW.23) makes the right side of (EW.26) `exp{-Omega(N)}`.  An explicit
sufficient asymptotic condition is

```math
{c\over2}\log{c\over e\bar a}
>{(2\lambda+2)(2\lambda+3)\over4}\bar\sigma,     \tag{EW.28}
```

where

```math
\bar a=\limsup_N{d\rho^2\over N},
\qquad \bar\sigma=\limsup_N{dt^2\over N},
```

with a strict margin and `c>2 bar a`.

For deterministic positivity repair, put

```math
p_-=(1-\rho)^d,qquad p_+=(1+\rho)^d
```

and define

```math
\widetilde P_K
=\min\{2p_+,\max\{p_-/2,P_{\le K}\}\}.           \tag{EW.29}
```

Every exact kernel, and hence `P`, lies in `[p_-,p_+]`.  On
`|R_K|<=P/2`, clipping is inactive and

```math
|\log P-\log\widetilde P_K|\le2{|R_K|\over P}.
```

The complementary event has `q_lambda`-probability at most `2A_K`, while
the clipped logarithmic difference is at most

```math
M_N=\log4+2dt.                                     \tag{EW.30}
```

Therefore

```math
E_{q_\lambda}|\log P-\log\widetilde P_K|
\le2(1+M_N)A_K=e^{-\Omega(N)},                    \tag{EW.31}
```

proving (EW.24).

For cavities define on the cube with edge `e` deleted

```math
D_e=E_{B_e}P,qquad
C_e=\rho^{-1}E_{B_e}(B_eP),qquad r_e={C_e\over D_e},           \tag{EW.32}
```

and define `D_(e,K),C_(e,K)` by replacing `P` with `P_(<=K)`.  Conditional
expectation is an `L^2` contraction, so

```math
\|D_e-D_{e,K}\|_2\le\sqrt{T_K},
\qquad
\|C_e-C_{e,K}\|_2\le\rho^{-1}\sqrt{T_K}.         \tag{EW.33}
```

Moreover

```math
D_e\in[(1-\rho)^{d-1},(1+\rho)^{d-1}],
\qquad P=D_e(1+\rho B_er_e).                      \tag{EW.34}
```

The marginal of `q_lambda` on `B_(-e)` is bounded above by

```math
{(1-\rho)^{-\lambda}\over Z_\lambda}D_e^{-\lambda}U_{-e}.
                                                               \tag{EW.35}
```

The negative-moment proof of Lemma EW.3 applies to `D_e` on `d-1`
coordinates.  Hence

```math
E_{q_\lambda}
 { |D_e-D_{e,K}|+|C_e-C_{e,K}|\over D_e}
\le(1-\rho)^{-\lambda}(1+\rho^{-1})\sqrt{T_K}
 e^{C_\lambda(d-1)t^2},                           \tag{EW.36}
```

where

```math
C_\lambda={(2\lambda+2)(2\lambda+3)\over4}.
```

Set

```math
\widetilde D_{e,K}
=\max\{(1-\rho)^{d-1}/2,D_{e,K}\},
\qquad
\widetilde r_{e,K}
=\operatorname {clip}_{[-1,1]}
 {C_{e,K}\over\widetilde D_{e,K}}.               \tag{EW.37}
```

On the event where the numerator in (EW.36) is at most `D_e/2`, elementary
ratio perturbation bounds the cavity error by a fixed multiple of that
numerator divided by `D_e`; off the event both clipped responses are bounded
by one, and Markov gives the same bound.  Equations (EW.27)--(EW.28), with
the harmless factors `d` and `rho^(-1)=O(sqrt N)`, prove (EW.25).  Finally,
`r -> atanh(rho r)` is `rho/(1-rho^2)`-Lipschitz.  `square`

The theorem supplies much more than the requested `o(N)` total error, but
at the wrong information scale.  Once
`K>=2 floor((N-1)/2)`, the retained coefficients
include all even subsets of a spanning tree of `K_(m,n)`.  By Theorem 37.71
these characters recover the complete `2^(N-2)`-atom projective latent law.
Thus (EW.24)--(EW.25) are a stable finite theorem, not a strict quotient.

## 6. Why the `L^2` obstruction is not an escort obstruction

Return to the antipodal one-atom law (EW.21), with fixed `lambda>0`.  Under
its inverse escort,

```math
{dq_\lambda\over dU}(B)
\propto\cosh(tS)^{-\lambda}.                      \tag{EW.38}
```

For comparable splits, `t sqrt d=Theta(sqrt N)`.  The central-binomial
upper and lower bounds, followed by a geometric sum on the lattice of `S`,
give

```math
\boxed{
E_{q_\lambda}\log\cosh(tS)=O_{\beta,\lambda}(1).} \tag{EW.39}
```

Indeed, the denominator in (EW.38) is
`Theta((t sqrt d)^(-1))`, and the numerator with the additional integrable
factor `log cosh(tS)` has the same order.  Consequently the degree-zero
quantity

```math
{d\over2}\log(1-\rho^2)                           \tag{EW.40}
```

approximates `log P` in inverse-escort `L^1` to `O(1)`.  Moreover every
exact cavity is a function of the single linear statistic

```math
r_e(B_{-e})=Q_e\tanh\left(t\sum_{f\ne e}B_fQ_f\right).         \tag{EW.41}
```

Thus this channel simultaneously has an extensive ordinary-`L^2` tail
above every `o(N)` degree and an exact one-feature cavity presentation under
the inverse escort.  It is not an actual diffuse child law, but it proves
that Theorem EW.2 cannot logically be promoted to an escort-weighted degree
lower bound.

There is also a universal actual-law warning that pressure-only
approximation is too weak.  Let `q_u` have density proportional to `P^{-u}`
for `u in[0,lambda]`.  The exact one-bit entropy estimate used in Theorem
37.52 gives

```math
D(q_u\Vert U)\le d\,\kappa(ut)\le{u^2dt^2\over2}.              \tag{EW.42}
```

If `phi(u)=log E_UP^{-u}`, then

```math
D(q_u\Vert U)=u\phi'(u)-\phi(u)
=\int_0^u v\operatorname {Var}_{q_v}(\log P)\,dv.              \tag{EW.43}
```

For every fixed `0<delta<lambda`, Cauchy--Schwarz therefore yields

```math
\boxed{
\int_\delta^\lambda
 E_{q_u}|\log P-E_{q_u}\log P|\,du
\le\left{{(\lambda-\delta)\lambda^2dt^2\over2\delta}
     \right}^{1/2}
=O_{\beta,\lambda,\delta}(\sqrt N).}             \tag{EW.44}
```

So a scalar already approximates the pressure to `o(N)` after averaging in
inverse temperature.  This gives neither a fixed-endpoint theorem nor a
cavity-response approximation.  It shows that the useful rare-event state
must preserve conditional/directional responses, not merely the value of
`log P` under its own escort.

## 7. Research conclusion

The degree audit has a sharp three-part answer.

1. **Degree `O(1)` or `o(N)` is impossible for ordinary `L^2(U)` likelihood
   approximation** at strong temperature, uniformly on the actual child
   laws.  The obstruction already follows from actual rank-one support
   collision and has the balanced threshold (EW.12).
2. **Degree `O(N)` is universally sufficient in inverse-escort `L^1`**, for
   both pressure and every edge cavity, after deterministic positivity
   repair.  It pays the unavoidable `exp{O(N)}` negative-moment factor and
   is full-information rather than compression.
3. **No escort-weighted `Omega(N)` degree lower bound follows.**  Inverse
   escorting can erase the high-pressure atoms that create the ordinary
   `L^2` obstruction.  The one-atom channel makes this exact, and (EW.44)
   shows that pressure-only self-escort approximation is generically a weak
   target.

Accordingly, a strict rare-event theorem must work directly with a
`q_lambda`-weighted **relative** Walsh/cavity error, or prove an
optimizer-specific negative-tail synchronization law.  Positive-replica
`L^2` control, even with the exact overlap enumerator (EW.7), cannot by
itself produce a sublinear-degree quotient at the strong temperatures
needed by the project.
