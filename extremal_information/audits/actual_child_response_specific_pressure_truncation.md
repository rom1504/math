# Response-specific truncation of the actual-child bridge pressure

Status: **rigorous actual-child truncation theorem; strict removal of the
ultra-rare external-pressure tail, not a complete response quotient**.

## 1. Target

The all-tilt superconcentration route fails because every actual child pair
has bridge atoms with pressure `Theta(N^(3/2))`.  Those atoms have probability
`exp(-Theta(N^2))`, however, while the product-factor problem only tests the
pressure against row products whose row `L^2` norms are uniformly bounded.
This note proves that the entire pressure above a fixed multiple of `N` is
uniformly negligible for every such product, including every exact
square-polynomial factor product from Theorem 37.39.

The proof uses only two optimizer-specific scalar facts:

1. actual child minimality bounds each zero-bridge pressure by its annealed
   value;
2. the exact product-cavity identities give a dimension-free row `L^2`
   bound for every optimal factor, and Theorem 37.39 preserves such a bound
   for the square carrier.

Thus the result applies to the actual optimizing children and no surrogate
child model.  It is response-specific: it controls `E_P L`, rather than the
all-parameter MGF of `L` under an external base law.

## 2. One annealed moment of the actual pressure

Let `A,D` be exact contracted-temperature pressure minimizers of orders
`m,n`, let `N=m+n`, and put `t=beta/sqrt(N)`.  In either relative
orientation define the raw bridge pressure

```math
L(B)=\log E_{x,z}\cosh\left(
t\{H_A(x)+\epsilon H_D(z)+x^{\mathsf T}Bz\}\right).
\tag{RT.1}
```

Write

```math
p_A(t)=\log E_x\cosh(tH_A(x)),
\qquad p_D(t)=\log E_z\cosh(tH_D(z)).
```

Actual minimality and averaging over all internal edge signings give

```math
p_A(t)\le {m\choose2}\log\cosh t,
\qquad
p_D(t)\le {n\choose2}\log\cosh t.               \tag{RT.2}
```

### Lemma RT.1 (replica annealed moment)

For every positive integer `k`, under the fair bridge law `U_(mn)`,

```math
\boxed{
\log E_Ue^{kL}
\le k\{\log2+p_A(t)+p_D(t)\}
   +{k^2t^2mn\over2}.}                           \tag{RT.3}
```

*Proof.*  Expand the `k`th power using independent replicas and introduce
one auxiliary sign `tau_l` in each replica.  Averaging a bridge edge gives

```math
E_{B_{ij}}\exp\left{
tB_{ij}\sum_{l=1}^k\tau_lx_i^lz_j^l\right}
=\cosh\left(t\sum_l\tau_lx_i^lz_j^l\right)
\le e^{k^2t^2/2}.
```

After all `mn` bridge edges are removed this way, the internal expectation
factorizes over replicas.  For one replica,

```math
E_{\tau,x,z}e^{t\tau(H_A+\epsilon H_D)}
\le2e^{p_A(t)+p_D(t)},
```

by `cosh(a+b)<=2 cosh(a)cosh(b)`.  Multiplying the `k` replica bounds proves
(RT.3). `square`

Taking `k=1`, using (RT.2), `log cosh t<=t^2/2`, and `mn<=N^2/4`, gives

```math
\boxed{
\log E_Ue^L\le\log2+{3\beta^2\over8}N.}          \tag{RT.4}
```

No external-field stability or all-tilt estimate is used.

## 3. Uniform clipping over every bounded-`L^2` row product

Let

```math
P=\bigotimes_{i=1}^m q_iU_n,
\qquad E_Uq_i=1,
\qquad \|q_i\|_2\le K,                           \tag{RT.5}
```

where `K` is fixed independently of `m,n` and the children.  Put

```math
L^{(C)}(B)=L(B)\wedge CN.
```

### Theorem RT.2 (uniform moderate-pressure recovery)

For every fixed `beta,K`, every

```math
C>{3\beta^2\over8}+2\log K                       \tag{RT.6}
```

admits `c=c(beta,K,C)>0` such that, uniformly over the actual children,
their orientation, and every row product (RT.5),

```math
\boxed{
0\le E_PL-E_PL^{(C)}
\le C_{\beta}N^{3/2}e^{-cN}=o(1).}               \tag{RT.7}
```

In particular the error is much smaller than the required `o(N)` response
scale.

*Proof.*  Since the total internal-plus-bridge Hamiltonian has at most
`N choose 2` unit terms,

```math
0\le L(B)\le t{N\choose2}\le {\beta\over2}N^{3/2}. \tag{RT.8}
```

Let `E_C={L>CN}`.  Equations (RT.4) and Markov's inequality give

```math
U(E_C)\le2\exp\left{-\left(C-{3\beta^2\over8}\right)N\right\}. \tag{RT.9}
```

The product likelihood `Q=prod_i q_i` satisfies
`||Q||_2<=K^m<=K^N`.  Therefore Cauchy--Schwarz and (RT.8)--(RT.9) imply

```math
\begin{aligned}
E_P(L-L^{(C)})
&\le E_P[L1_{E_C}]\\
&\le\|Q\|_2\{E_U[L^21_{E_C}]\}^{1/2}\\
&\le C_\beta N^{3/2}
 \exp\left{N\log K
 -{N\over2}\left(C-{3\beta^2\over8}\right)\right\}.
\end{aligned}                                    \tag{RT.10}
```

The exponent is negative by (RT.6), proving (RT.7). `square`

The theorem applies simultaneously to:

- every globally optimal product factor from the exact cavity fixed-point
  equations, because the coordinate best-response equation gives directly
  `D_2(q_iU_n||U_n)<=lambda^2 beta^2 n/N`;
- every factor in the exact degree-`2d` square carrier used by Theorem
  37.39, because that carrier has one fixed row `L^2` bound `K_1`.

Consequently, if `V^row` and `V^(d,sq)` are the exact product variational
values and `V_C^row,V_C^(d,sq)` are obtained by replacing `L` by
`L^(C)`, then

```math
\boxed{
|V^{\rm row}-V_C^{\rm row}|
+|V^{(d,{\rm sq})}-V_C^{(d,{\rm sq})}|
\le C_\beta N^{3/2}e^{-cN}.}                    \tag{RT.11}
```

The entropy terms are unchanged.  Thus every extensive reverse-product or
coherent-retuning phase is already present in the moderate-pressure
landscape `L<=CN`.

For completeness, the clipped unrestricted optimum also has the same row
collision bound: composition with `x mapsto x wedge CN` does not increase a
one-edge oscillation, so its coordinate best-response equation gives the
identical Renyi-two estimate.  Evaluating the clipped objective at the
original optimizer and the original objective at the clipped optimizer,
and applying (RT.7) in both directions, proves (RT.11).  The square-carrier
case is immediate because its entire feasible class already satisfies
(RT.5).

## 4. Exact scope and remaining assumption

The low-information **certificate** for (RT.7) consists only of the two
child pressure bounds (RT.2), the row collision bound `K`, and the cutoff
constant `C`.  All are verifiable from proved child minimality and cavity
identities.  In particular, the `Theta(N^(3/2))` rank-one atom which
falsifies all-tilt superconcentration contributes exponentially little to
every admissible product response.  The rare-atom falsifier does not survive
the actual product-factor optimization.

What RT.2 does **not** provide is a low-information representation of the
retained function `L^(C)`.  Its moderate level sets may still carry the full
cross-row response tensor.  Quantizing those level sets would merely rename
the pressure table.  Therefore the next response-specific statement must
concern only the clipped actual pressure, for example:

> **Moderate-pressure actual-child quotient.**  Construct from polynomially
> many child observables a function `L_tilde` such that, uniformly over the
> fixed-degree square carrier products,
> `|E_P L^(C)-E_P L_tilde|=o(N)`.

This is strictly narrower than the prior target: ultra-rare external
pressures, arbitrary all-tilt MGFs, high row Walsh degree, positivity, and
row entropy have all been removed.  It is not proved by RT.2, and no current
minimality identity controls the moderate-pressure cross-row level sets.

The theorem is therefore a response-specific truncation reset, not a
Level-6 closure.
