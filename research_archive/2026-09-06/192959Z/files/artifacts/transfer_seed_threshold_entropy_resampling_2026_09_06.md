# Flat Gaussian threshold entropy by smoothing and coordinate resampling

Date: 2026-09-06. Seed-transfer track. This proves a uniform discrete
entropy upper bound with explicit parameters improving the clipped-law
credit. No Gaussian-replica universality or higher-Hermite truncation is used.

## 1. Uniform finite theorem

Let `A` be an order-`N` hollow symmetric signing, `B=A/sqrt(N)`, and
`||B||op<=C`. Fix `delta in (0,1/2)` and write

```math
t=\Phi^{-1}(\delta),\quad m=1-2\delta,\quad
v=4\delta(1-\delta),\quad a=2\phi(t),\quad \eta=a^2/v.
```

Choose real `s` and smoothing parameter `r` with `0<|s|C<r<1`.
Let `p_s` be the law of `X_i=sign(G_i-t)` for
`G~N(0,I+sB)`, and let `p_0` be the product law with mean `m`. Define

```math
\kappa=\eta\int_0^r
 \frac{\exp(t^2z/(1+z))}{\sqrt{1-z^2}}\,dz,\qquad
\lambda=\frac{s^2}{2(r^2-s^2C^2)},                        (1)
```

and

```math
F_B(\lambda,\kappa)=\frac1{(1-\kappa)^2}
 \exp\left(\frac{4\lambda\kappa e^{2\lambda}}{1-\kappa}\right).
```

Then uniformly over all eligible actual signings,

```math
\boxed{D(p_s\Vert p_0)\le\log(1+\chi^2(p_s\Vert p_0))
\le N\lambda\kappa^2F_B(\lambda,\kappa)+\log(N+1).}       (2)
```

All constants are explicit and finite at fixed parameters. The exact
intermediate comparison is

```math
1+\chi^2(p_s\Vert p_0)
\le\mathbb E_{S\sim\operatorname{Ber}(\kappa)^{\otimes N}}
 \det(I-(s/r)^2B_S^2)^{-1/2}
\le\mathbb E_{K\sim\operatorname{Bin}(N,\kappa)}
 e^{\lambda K(K-1)/N}.                                  (3)
```

The second inequality uses flat entries. The operator comparison preceding
it is valid more generally and is not a pointwise comparison of pair laws.

## 2. Gaussian smoothing and a rank-one operator bound

Let `gamma_N` be standard Gaussian measure and
`L_u=dN(0,I+uB)/d gamma_N`. The condition `|s|C/r<1` ensures
`L_(s/r) in L^2(gamma_N)`. Define the Gaussian smoothing operator

```math
(T_\alpha f)(x)=\mathbb E_Z f(\alpha x+\sqrt{1-\alpha^2}Z).
```

It is self-adjoint in Gaussian `L^2`. Convolution changes a Gaussian
covariance `R` to `alpha^2 R+(1-alpha^2)I`; consequently

```math
L_s=T_{\sqrt r}^{\otimes N}L_{s/r}.                       (4)
```

Let `P` be conditional expectation onto one threshold bit. The functions
`1` and `u(x)=(sign(x-t)-m)/sqrt(v)` are orthonormal, so
`P=|1><1|+|u><u|`. Thresholding a likelihood is conditional expectation,
which gives the exact identity

```math
1+\chi^2(p_s\Vert p_0)
=\|P^{\otimes N}L_s\|_2^2
=\langle L_{s/r},(T_{\sqrt r}PT_{\sqrt r})^{\otimes N}L_{s/r}\rangle.
                                                               (5)
```

Put `w=T_(sqrt(r))u`. It is centered and its squared norm is the
normalized covariance of threshold bits at Gaussian correlation `r`.
The derivative of the bivariate threshold probability is
`phi(t)^2 exp[t^2r/(1+r)]/sqrt(1-r^2)`. Integration from zero gives
`||w||_2^2=kappa` from (1), in particular `0<kappa<1`.
For a direct calculus verification, its Gaussian density satisfies
`partial_r phi_r(x,y)=partial_x partial_y phi_r(x,y)`; integrating on
`(-infinity,t]^2` leaves the density at `(t,t)`.
The same identity follows from the convergent full Hermite expansion.

Rank-one Cauchy--Schwarz on the centered subspace now gives

```math
0\le T_{\sqrt r}PT_{\sqrt r}
\le K_\kappa:=|1\rangle\langle1|
                  +\kappa(I-|1\rangle\langle1|).          (6)
```

These are positive-semidefinite OPERATOR inequalities. Their tensor powers
preserve the ordering: expand the difference as a sum of positive tensor
products with one factor `K_kappa-TPT`. Every higher Hermite is retained.

## 3. Coordinate resampling and the exact use of flatness

The operator `K_kappa` keeps a coordinate with probability `kappa`
and averages it away otherwise. Its product is the average of orthogonal
projections conditioning on subsets `S~Ber(kappa)^(tensor N)`. Thus

```math
1+\chi^2(p_s\Vert p_0)
\le\mathbb E_S\|\mathbb E[L_{s/r}\mid G_S]\|_2^2.
```

The conditional density is exactly the Gaussian marginal likelihood
with covariance `I+(s/r)B_S`. Direct Gaussian integration gives

```math
\|L_{I+H}\|_2^2=\det(I-H^2)^{-1/2},\qquad\|H\|_{op}<1. (7)
```

Indeed its two determinant factors are `det(I+H)^(-1)` and
`det(2(I+H)^(-1)-I)^(-1/2)`, whose product is (7).
This proves the first inequality of (3), including the precise `L^2`
domain and square-root power.

Put `rho=|s|C/r<1`. Principal operator monotonicity and
`-log(1-x)<=x/(1-rho^2)` for `0<=x<=rho^2` give a determinant
upper bound in terms of `tr B_S^2`. Flatness is used exactly here:

```math
\operatorname{tr}B_S^2=|S|(|S|-1)/N.
```

The logarithm of (7) is therefore at most
`lambda |S|(|S|-1)/N`, proving the second inequality of (3).
An operator-only bound does not give this cardinality identity. This
does not rule out other, weaker nonflat bounds with explicit dependence
on the threshold probability.

## 4. An explicit Bernoulli count bound

The elementary binomial entropy bound gives

```math
\mathbb E e^{\lambda K(K-1)/N}
\le(N+1)\exp\left[N\sup_{0\le x\le1}
 \{\lambda x^2-D(\operatorname{Ber}(x)\Vert
                         \operatorname{Ber}(\kappa))\}\right]. (8)
```

The maximizer is interior: the derivative tends to opposite infinities
at the endpoints. Its stationary equation is

```math
\frac{x}{1-x}=\frac\kappa{1-\kappa}e^{2\lambda x}.
```

First `x<=b:=kappa e^(2lambda)/(1-kappa)`, and then the same equation
gives `x<=kappa exp(2lambda b)/(1-kappa)`. Discarding the nonnegative
binary divergence bounds the supremum by `lambda kappa^2 F_B`.
This proves (2). Finally `D<=log(1+chi^2)` is Jensen under `p_s`.
There is no saddle-point approximation; at rare `kappa` this elementary
factor tends to one.

## 5. A rational coefficient at the current flip density

Choose

```math
\delta=2^{-24},\quad C=33/32,\quad s=1/32,\quad r=1/25.
```

Then the following convenient consequence holds for every eligible matrix:

```math
\boxed{D(p_s\Vert p_0)\le22\frac{\eta^2s^2}{4}N+\log(N+1).} (9)
```

All constants have rational certificates. Put `z=-t`. Mills bounds
`phi(z)z/(1+z^2)<=Phi(-z)<=phi(z)/z`, together with
`8/3<e<11/4` and `3<pi<22/7`, give

```math
\Phi(-5)>\frac5{78(11/4)^{13}}>2^{-24},\qquad
\Phi(-\sqrt{29})<\frac5{104(8/3)^{14}}<2^{-24}.
```

The bounds on `e` follow from its positive series and a geometric tail:
the first four terms equal `8/3`, and the sum from degree three onward
is at most `(1/6)/(1-1/4)=2/9`, giving `e<49/18<11/4`.
The Mills lower bound follows by integration by parts and
`integral_z^infinity phi(u)/u^2 du<=Phi(-z)/z^2`; the upper bound
is the same integration-by-parts identity with its positive remainder.
Thus `5<z<sqrt(29)`. With `x=a^2/v^2`, Mills also gives

```math
25/4<x<\frac{225}{29(1-\delta)^2}<8,\qquad \eta=xv<32\delta. (10)
```

Since `z^2<=29`,

```math
\frac\kappa{\eta r}
\le\frac{e^{29r}-1}{29r\sqrt{1-r^2}}<19/10.               (11)
```

A rational Taylor sum with geometric remainder proves
`e^(29/25)<16/5`; then `(55/29)^2(625/624)<(19/10)^2`
proves (11). Moreover `lambda=320000/367951<1` and `kappa<10^-6`.
Thus `b=kappa e^(2lambda)/(1-kappa)<9*10^-6` and

```math
F_B<[(1-10^{-6})^2(1-36\cdot10^{-6})]^{-1}<1001/1000.
```

Here `e^y<=1/(1-y)` for `0<=y<1`. The coefficient of
`eta^2s^2/4` in (2) is less than

```math
\frac{2(19/10)^2(1001/1000)}{1-(sC/r)^2}
=\frac{23682154496}{1149846875}<22.                       (12)
```

## 6. Direct pressure consequence, with all error terms

Take physical inverse temperature `b_phys=beta/sqrt(N)`, where
`8<=beta<=33/4`. Switch an absolute maximizing spin to all ones and
choose its energy sign. The resulting signing `A'` has unchanged
flatness and norm and satisfies `q_(A')(1)=Q(A)`. Apply the threshold
law to covariance `I+s A'/sqrt(N)`.

The marginals have mean `m`. The full threshold Hermite expansion and
Parseval give, for distinct coordinates of correlation `rho_ij`,

```math
\mathbb E X_iX_j=m^2+a^2\rho_{ij}+R_{ij},\qquad
|R_{ij}|\le v\rho_{ij}^2.
```

Here `rho_ij=s A'_(ij)/sqrt(N)`. After multiplying energy by
`beta/sqrt(N)`, the total remainder is at most
`beta v s^2 sqrt(N)/2`. Also `H(X)=N h(delta)-D(p_s||p_0)` exactly.
Gibbs variational therefore proves

```math
\log(Z_+(b_{phys})+Z_-(b_{phys}))
\ge\frac{\beta m^2 Q(A)}{\sqrt N}
 +N\{h(\delta)+g_{new}\}
 -\frac{\beta a^2s}{2}-\frac{\beta v s^2\sqrt N}{2}
 -\log(N+1),                                            (13)
```

where

```math
g_{new}=\frac{\beta a^2s}{2}-22\frac{\eta^2s^2}{4}.
```

All errors are explicit and subextensive at the fixed parameters.
The law is built after switching the actual maximizer, not assumed
independent of that maximizer. Equation (13) is deterministic for each
eligible signing.

The previous clipped-law credit is `g_clip=beta^2 v^2/128`.
Using (10), `s=1/32`, and `8<=beta<=33/4`,

```math
\frac{g_{new}}{v^2}
=\frac{\beta x}{64}-\frac{11x^2}{2048}
\ge\frac{x}{8}-\frac{11x^2}{2048}
\ge\frac{18725}{32768},\qquad
\frac{g_{clip}}{v^2}\le\frac{1089}{2048}.
```

The quadratic is increasing on `[25/4,8]`. Hence the certified strict
improvement is

```math
g_{new}-g_{clip}\ge\frac{1301}{32768}v^2>0.               (14)
```

Retained weaves at `p=31/32,t_weave=4` satisfy
`beta=8/sqrt(p)+o(1)` and normalized operator norm
`1/sqrt(p)+o(1)`, within the strict rational ranges above for large orders.
Substituting (13) into the existing annealed extraction strictly improves
its clipped-law upper coefficient. The director owns that final integration
and its ledger value; no new unaudited decimal endpoint is asserted here.

## Reproduction and scope

The proof retains every rare-threshold and localized-Gram effect through
a PSD operator comparison to coordinate resampling. It does not compare
the replica law to a Gaussian of correlation `eta`.

The checker is
`computations/transfer_seed_threshold_entropy_resampling_verify_2026_09_06.py`.
It replays rational constants and finite low-dimensional threshold and
determinant tests. No input or output files are required.
