# An extensive Jensen-certificate barrier for the actual child bridge escort

Status: **task-local rigorous theorem note**.  This note analyzes the most
natural attempt to transfer the forward child-spin product-channel latent to
the negative-disorder bridge escort.  The attempt admits an exact
product-mixture proposal with uniformly bounded row Renyi-two complexity,
but the direct convexity domination certificate pays `Theta(N)` at every
comparable split.  Thus this canonical *Jensen transfer* does not establish
the desired tight-component closure.  The theorem is not a lower bound on
the optimal rejection constant or on the distance between the two laws.
It applies to the actual optimizing
children, although its proof only needs the global-spin symmetry shared by
all children.

The normalization is that of
[`actual_child_negative_escort_structure.md`](actual_child_negative_escort_structure.md).

## 1. Conditional forward channel and inverse escort

Fix children `A,D`, raw temperature `t>0`, and an orientation
`epsilon in {+-1}`.  Put `L=mn`.  Conditional on `epsilon`, let `mu_epsilon`
be the law of the augmented child Gibbs variables and define their rank-one
bridge word

```math
Q_{ij}=\tau_1x_iy_j.
```

With `rho=tanh t`, the forward bridge-channel density relative to the fair
bridge law `U_B` is

```math
k_Q(B)=\prod_{i,j}(1+\rho Q_{ij}B_{ij})
=\prod_{i,j}{e^{tQ_{ij}B_{ij}}\over\cosh t}.           \tag{IM.1}
```

Consequently

```math
p_\epsilon(B):={d\Pi_\epsilon\over dU_B}(B)
=\mathbb E_{\mu_\epsilon}k_Q(B),                       \tag{IM.2}
```

and the negative escort has density

```math
q_{\lambda,\epsilon}(B)
={p_\epsilon(B)^{-\lambda}\over Z_{\lambda,\epsilon}},
\qquad
Z_{\lambda,\epsilon}=\mathbb E_{U_B}p_\epsilon^{-\lambda},
\qquad\lambda>0.                                      \tag{IM.3}
```

## 2. The exact reversed-channel product proposal

Define

```math
c_\lambda(t)
=\mathbb E_{b\sim U_1}(1+\rho b)^{-\lambda}
=(\cosh t)^\lambda\cosh(\lambda t).                   \tag{IM.4}
```

For every latent word `Q`, normalized coordinatewise inversion of its
forward channel gives

```math
\widetilde k_Q(B)
={k_Q(B)^{-\lambda}\over c_\lambda(t)^L}
=\prod_{i,j}{e^{-\lambda tQ_{ij}B_{ij}}\over\cosh(\lambda t)}.
                                                               \tag{IM.5}
```

Thus `widetilde k_Q` is again an exact product channel, now with coordinate
mean `-tanh(lambda t)Q_ij`.  Let

```math
r_{\lambda,\epsilon}(B)
=\mathbb E_{\mu_\epsilon}\widetilde k_Q(B).            \tag{IM.6}
```

This is a common-latent mixture of bit-product, hence row-product, bridge
laws.  Every component row has the exact Renyi-two complexity

```math
D_2(\widetilde k_{Q,i}\Vert U_n)
=n\log(1+\tanh^2(\lambda t))
\le n\lambda^2t^2.                                     \tag{IM.7}
```

It is therefore uniformly tight at fixed `beta,lambda` and comparable
splits when `t=beta/sqrt(N)`.

Convexity of `u mapsto u^{-lambda}` gives pointwise

```math
p_\epsilon(B)^{-\lambda}
\le \mathbb E_{\mu_\epsilon}k_Q(B)^{-\lambda}
=c_\lambda(t)^Lr_{\lambda,\epsilon}(B).                \tag{IM.8}
```

Define the inverse-mixture Jensen gap

```math
J_{\lambda,\epsilon}
=L\log c_\lambda(t)-\log Z_{\lambda,\epsilon}.         \tag{IM.9}
```

Then

```math
q_{\lambda,\epsilon}(B)
\le e^{J_{\lambda,\epsilon}}
 r_{\lambda,\epsilon}(B),
\qquad
D(q_{\lambda,\epsilon}\Vert r_{\lambda,\epsilon})
\le J_{\lambda,\epsilon}.                             \tag{IM.10}
```

In particular, the pointwise Jensen certificate supplies the following
(possibly very inefficient) exact rejection scheme: sample `r`, accept `B`
with probability

```math
e^{-J}{q(B)\over r(B)},                                \tag{IM.11}
```

and obtain `q` conditional on acceptance; this particular scheme has
acceptance probability exactly `e^{-J}`.  The optimal rejection scheme can
be much better because (IM.10) only proves
`D_infinity(q||r)<=J`; it does not prove equality.

## 3. Central symmetry forces an extensive gap

**Theorem IM.1 (extensive Jensen-certificate barrier).**  For every finite
child pair, every orientation, and every `t,lambda>0`,

```math
\boxed{
L\log\cosh(\lambda t)
\le J_{\lambda,\epsilon}
\le L\{\lambda\log\cosh t+\log\cosh(\lambda t)\}.}    \tag{IM.12}
```

In particular, if `m/N -> theta in (0,1)` and `t=beta/sqrt(N)`, then

```math
{\lambda^2\beta^2\over2}\theta(1-\theta)
\le\liminf {J_{\lambda,\epsilon}\over N}
\le\limsup {J_{\lambda,\epsilon}\over N}
\le{(\lambda+\lambda^2)\beta^2\over2}\theta(1-\theta).
                                                               \tag{IM.13}
```

Thus the direct Jensen domination of the negative escort by the canonical
reversed-channel proposal supplies only an `exp[-Theta(N)]` acceptance
certificate.  This statement alone does **not** show that every transfer or
even the optimal rejection coupling has exponential cost.

*Proof.*  Conditional on `epsilon`, the map `x mapsto -x` preserves the
child Gibbs weight and sends `Q mapsto -Q`.  Hence `mu_epsilon` is centrally
symmetric.  Decompose it into orbits `{Q,-Q}`.  For one such orbit, put

```math
S_Q(B)=\sum_{i,j}Q_{ij}B_{ij}.
```

Its forward mixture density is

```math
p_Q^{\rm pair}(B)
={k_Q(B)+k_{-Q}(B)\over2}
={\cosh(tS_Q(B))\over(\cosh t)^L}.                     \tag{IM.14}
```

Under `U_B`, `S_Q(B)` has the law of a sum `S_L` of `L` fair signs,
independently of `Q`.  Therefore every orbit has the same inverse moment

```math
Z_{\rm pair}
=(\cosh t)^{\lambda L}
 \mathbb E\cosh(tS_L)^{-\lambda}.                      \tag{IM.15}
```

If `p_epsilon=sum_a w_a p_a^pair`, convexity gives

```math
Z_{\lambda,\epsilon}
=\mathbb E(p_\epsilon)^{-\lambda}
\le\sum_aw_a\mathbb E(p_a^{\rm pair})^{-\lambda}
=Z_{\rm pair}.                                        \tag{IM.16}
```

Combining (IM.4), (IM.9), and (IM.15),

```math
J_{\lambda,\epsilon}
\ge L\log\cosh(\lambda t)
 -\log\mathbb E\cosh(tS_L)^{-\lambda}
\ge L\log\cosh(\lambda t),                           \tag{IM.17}
```

because the final expectation is at most one.  Conversely
`E_U p_epsilon=1`, so Jensen gives
`Z_(lambda,epsilon)>=1`.  Equations (IM.4) and (IM.9) give the upper bound.
Taylor expansion at `t=beta/sqrt(N)` proves (IM.13). `square`

## 4. A genuine separation at sufficiently low physical temperature

The preceding Jensen gap is only a certificate cost.  Rank-one support does,
however, give a genuine separation in a definite parameter regime.

**Theorem IM.2 (large-`beta` separation from the canonical proposal).**  Let
`K_epsilon` be the number of rank-one words in the support of
`mu_epsilon`.  At finite temperature the latent law has full support, so
`K_epsilon=2^(m+n-1)`.  In particular,

```math
\boxed{
D(U_B\Vert\Pi_\epsilon)
\ge L\log\cosh t-t\sqrt{2L\log K_\epsilon},
\qquad K_\epsilon=2^{m+n-1}.}                          \tag{IM.18}
```

At escort parameter `lambda=1`, central symmetry gives

```math
r_{1,\epsilon}=p_\epsilon.                             \tag{IM.19}
```

Therefore

```math
\boxed{
D(q_{1,\epsilon}\Vert r_{1,\epsilon})
\ge D(U_B\Vert p_\epsilon),
\qquad
\|q_{1,\epsilon}-r_{1,\epsilon}\|_{\rm TV}
\ge1-\exp\{-D(U_B\Vert p_\epsilon)/2\}.}              \tag{IM.20}
```

If `m/N->theta`, `t=beta/sqrt(N)`, and

```math
\gamma(\beta,\theta)
={\beta^2\over2}\theta(1-\theta)
-\beta\sqrt{2\theta(1-\theta)\log2}>0,                \tag{IM.21}
```

then the KL in (IM.20) is at least
`(gamma(beta,theta)-o(1))N` and the total variation tends to one
exponentially.  Equivalently, this conclusion holds when

```math
\beta>\sqrt{8\log2\over\theta(1-\theta)}.              \tag{IM.22}
```

*Proof.*  Every supported word is rank one.  Conversely, every rank-one word
`uv^T` occurs by taking, for example, `tau_1=1,x=u,y=v,tau_2=epsilon`;
all such states have positive finite-temperature Gibbs weight.  The usual
simultaneous sign redundancy leaves exactly `2^(m+n-1)` distinct words.
For a fair bridge put
`S_Q=sum_eQ_eB_e`.  Each `S_Q` is a centered `L`-subgaussian variable, and
the standard exponential-max bound gives

```math
\mathbb E_{U_B}\max_{Q\in\operatorname{supp}\mu}S_Q
\le\sqrt{2L\log K_\epsilon}.                           \tag{IM.23}
```

Equation (IM.2) implies pointwise

```math
\log p_\epsilon(B)
\le t\max_QS_Q(B)-L\log\cosh t.
```

Taking the negative uniform expectation proves (IM.18).

For `lambda=1`, (IM.5) is exactly
`widetilde k_Q=k_(-Q)`.  Central symmetry of `mu_epsilon` proves (IM.19).
Writing `Z=E_Up_epsilon^{-1}`, the Hellinger affinity between `q_1` and
`p_epsilon` is exactly

```math
\mathbb E_U\sqrt{q_1p_\epsilon}=Z^{-1/2}.              \tag{IM.24}
```

Jensen gives `log Z>=E_U[-log p]=D(U||p)`.  Monotonicity of Renyi
divergence gives
`D(q_1||p)>=-2log(Z^(-1/2))=log Z`, while
`1-TV<=` the Hellinger affinity.  This proves (IM.20).  Substitution of
`L=theta(1-theta)N^2+o(N^2)` and
`log K<=(N-1)log2` proves (IM.21)--(IM.22). `square`

This is an actual distance theorem, not just a large Jensen certificate.
Its threshold is deliberately not claimed sharp.  Below (IM.22), deciding
whether rank-one channel resolvability makes `q` close to a product mixture
is a covering/resolvability problem rather than a consequence of the present
entropy bound.

## 5. What the barrier does and does not prove

The theorem is a sharp falsifier of the **canonical** latent-product route:

- its proposed components have uniformly bounded row `D_2` by (IM.7);
- nevertheless the direct Jensen domination of the actual inverse escort
  pays a linear log-density/rejection *certificate* by (IM.12);
- this obstruction holds for the actual optimized children and cannot be
  removed by selecting a better child minimizer.

The gap `J` is not ordinary row total correlation and does not by itself
lower-bound `D(q||r)`, total variation, the optimal rejection constant, or
either KL projection onto row-product laws.  That limitation is real.  If
the latent distribution were uniform on all bridge words, both the forward
and reversed mixtures would be uniform, even though the same Jensen gap is
extensive.  The actual latent words are rank one, but turning that support
restriction into a genuine separation between `q` and `r`, or into an
extensive row-dependence lower bound, is a separate theorem.

Together with the exact product-shadow identity AC.24, the result narrows
the structural question to the following alternatives:

1. control the best row-product variational term directly by actual-child
   minimality; or
2. prove that the directed product-projection cost in AC.16 is extensive.

The naive convexity comparison to the forward Gibbs latent and its
reversed-channel image cannot decide this dichotomy at sublinear certified
cost.  It remains possible that a sharper comparison shows the two marginal
laws close.  In particular, tight sequential row `D_2` is true but
insufficient, while tight common-latent product complexity remains unproved.
