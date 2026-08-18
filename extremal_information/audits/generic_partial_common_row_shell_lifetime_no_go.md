# Partial-common-row obstruction to a generic shell-to-lifetime inequality

Status: **rigorous scalable generic counterexample**. This is not an actual
minimizing-child channel. Its right projective factor is fixed, so it
violates factorwise quadratic spread even though the global rank-one prior
has linear min-entropy, exponentially small collision, and exponentially
decaying narrow rank-one caps.

The example shows that global latent spread and the averaged-posterior
shell/geometric split cannot generically control `J-I^leftarrow`.

## 1. A partially shared rank-one prior

Let `m=k+ell`, let `m+n=N`, and assume

```math
 {k\over N}\to\kappa>0,\qquad
 {\ell\over N}\to\zeta>0,\qquad
 {n\over N}\to\nu>0.                                 \tag{PC.1}
```

Take independent fair signs `sigma,xi_1,...,xi_ell` and put

```math
 X=(\sigma\mathbf1_k,\xi_1,\ldots,\xi_\ell),
 \qquad Q=X\mathbf1_n^{\mathsf T}.                    \tag{PC.2}
```

The latent support is the row-sign subgroup

```math
 G=\{(\sigma\mathbf1_k,\xi):
       \sigma\in\{+-1\},\ \xi\in\{+-1\}^\ell\}
 \le\{+-1\}^m.                                       \tag{PC.3}
```

The map `X -> Q` is injective, so the uniform prior `mu` has

```math
 \boxed{
 \|\mu\|_\infty=\sum_Q\mu(Q)^2=2^{-(\ell+1)}
 =e^{-\zeta(\log2)N+o(N)}.}                          \tag{PC.4}
```

It also has geometric spread. For every fixed rank-one word `uv^T`,

```math
 { |\langle Q,uv^T\rangle|\over mn}
 ={ |\langle X,u\rangle|\over m}
  { |\langle\mathbf1_n,v\rangle|\over n}.            \tag{PC.5}
```

If the left side is at least `1-delta`, then
`|<X,u>|>=(1-delta)m`. The common block contributes at most `k`, so the
`ell` independent signs must have absolute sum at least `ell-delta m`.
Consequently, whenever `0<delta<ell/m`,

```math
 \boxed{
 \mu\left\{Q:
 { |\langle Q,uv^T\rangle|\over mn}\ge1-\delta\right\}
 \le {2\over2^\ell}
   \sum_{r\le\lfloor\delta m/2\rfloor}{\ell\choose r}.} \tag{PC.6}
```

If `ell/m->a in (0,1)`, the exponent is at most

```math
 -\ell\left\{\log2-h\left({\delta\over2a}\right)\right\}
 +o(m)<0                                               \tag{PC.7}
```

for every fixed `delta<a`.

## 2. Exact likelihood factorization

Let `u=beta/sqrt(N)`, let `U` be the fair bridge law, and write

```math
 S_i(B)=\sum_{j=1}^nB_{ij},
 \qquad S_C(B)=\sum_{i=1}^kS_i(B).                    \tag{PC.8}
```

Independence of the latent signs gives the exact forward likelihood

```math
 \boxed{
 p(B)=p_C(B_C)\prod_{i=k+1}^mp_i(B_i),}               \tag{PC.9}
```

where

```math
 p_C={\cosh(uS_C)\over(\cosh u)^{kn}},
 \qquad p_i={\cosh(uS_i)\over(\cosh u)^n}.            \tag{PC.10}
```

For fixed `lambda>0`, let

```math
 {dq\over dU}={p^{-\lambda}\over E_Up^{-\lambda}}.   \tag{PC.11}
```

Then

```math
 \boxed{
 q=q_C\otimes\bigotimes_{i=k+1}^mq_i,}               \tag{PC.12}
```

where `q_C` is the inverse escort of the `k`-row common-sign channel and
every `q_i` is the inverse escort of one isolated row. The canonical
inverse row product is

```math
 r=\bigotimes_{i=1}^mq_i.                             \tag{PC.13}
```

Therefore product additivity gives exactly

```math
 \boxed{
 J(q,r)=J(q_C,r_C),\qquad
 I^\leftarrow(q)=I^\leftarrow(q_C).}                 \tag{PC.14}
```

Indeed, for every row product `P`,

```math
 D(P\Vert q)
 =D(P_C\Vert q_C)+\sum_{i>k}D(P_i\Vert q_i),          \tag{PC.15}
```

and the nuisance minima occur at `P_i=q_i`.

## 3. Linear canonical retuning survives global spread

Let `G_0~N(0,beta^2 nu)`, `f(x)=log cosh x`, and put

```math
 z_0=Ee^{-\lambda f(G_0)},
 \qquad
 d_0=-\lambda{E[f(G_0)e^{-\lambda f(G_0)}]\over z_0}
      -\log z_0>0.                                    \tag{PC.16}
```

The proof of Theorem 37.58 applies verbatim with ambient scale `N` and
dimensions `k,n`. Explicitly,

```math
 J(q_C,r_C)
 =kD(q_i\Vert U_n)
  +\lambda E_{r_C}\log\cosh\left(u\sum_{i=1}^kS_i\right)
  +\log E_{U_{kn}}\cosh(uS_C)^{-\lambda}.             \tag{PC.16a}
```

The first term is `k d_0+o(N)` by the CLT, while each remaining term is
`O(sqrt(N))` by the second-moment and central-event estimates in that proof.
Together with (PC.14), this gives

```math
 \boxed{J=\kappa d_0N+o(N).}                          \tag{PC.17}
```

Take the fair row product on the common block and the exact factors `q_i`
on nuisance rows. Then

```math
 I^\leftarrow(q)
 \le D(U_{kn}\Vert q_C)
 \le\lambda u\sqrt{kn}=O(\sqrt N).                  \tag{PC.18}
```

The last estimate follows from

```math
 D(U\Vert q_C)
 =\lambda E_U\log\cosh(uS_C)
  +\log E_U\cosh(uS_C)^{-\lambda},
```

whose second term is nonpositive, together with
`log cosh z<=|z|` and `E|S_C|<=sqrt(kn)`. Consequently

```math
 \boxed{J-I^\leftarrow=\kappa d_0N+o(N).}             \tag{PC.19}
```

## 4. Averaged posterior retuning is exactly zero

For `g in G`, let `T_g` switch bridge row `i` by `g_i`. The fair bridge law
is invariant under `T_g`. Uniformity of the subgroup prior gives

```math
 p(T_gB)=p(B),\qquad k_X(T_gB)=k_{Xg}(B),             \tag{PC.20}
```

where `k_X` is the channel density conditional on `X`.

Let `w` be any normalized bridge law whose density relative to `U` is a
function of `p(B)`; this includes every raw disorder tilt and the inverse
escort (PC.11). Its averaged forward Bayes posterior is

```math
 \bar\mu(X)=\mu(X)E_w{k_X(B)\over p(B)}.              \tag{PC.21}
```

Changing variables `B -> T_gB` in (PC.21) shows that the multiplier is the
same for `X` and `Xg`. The group acts transitively on its support, so the
multiplier is constant. Normalization proves

```math
 \boxed{\bar\mu=\mu.}                                 \tag{PC.22}
```

It follows that for every deterministic quotient `F(Q)`, including any
energy-shell analogue,

```math
 \boxed{
 D(\bar\mu\Vert\mu)=0,
 \quad D(F_\#\bar\mu\Vert F_\#\mu)=0,
 \quad
 \sum_f\bar\mu(F=f)
 D(\bar\mu(\cdot|f)\Vert\mu(\cdot|f))=0.}            \tag{PC.23}
```

Thus both terms of every shell/geometric retuning split vanish while
`J-I^leftarrow=Theta(N)`.

## 5. Exact scope

The example proves the generic no-go

```math
 \boxed{
 \text{global exponential prior spread + zero averaged-posterior retuning}
 \centernot\Longrightarrow J-I^\leftarrow=o(N).}       \tag{PC.24}
```

It rules out every generic inequality controlling `J-I^leftarrow` by the
shell KL plus within-shell deficit of the averaged latent posterior, even
after assuming exponential maximum-atom, collision, and fixed-cap decay.

This does not falsify an actual-child theorem. The right factor in
`Q=X1_n^T` is a point mass and violates Theorem 37.61 factorwise. The
counterexample says that Theorem 37.62 is a classification of latent
retuning, not a proxy for canonical row-product regret. Any viable
`L_shell-to-row-lifetime` theorem must use factorwise child spread or a
stronger optimizer identity linking posterior geometry to the canonical
row factors; global rank-one spread is insufficient.
