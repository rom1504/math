# Double-partial rank-one obstruction: whole-factor spread still does not control row lifetime

Status: **rigorous scalable generic counterexample**.  This is not an
actual minimizing-child channel.  Both projective factors have exponential
whole-word min-entropy and the averaged posterior is exactly unchanged by
every likelihood-dependent disorder tilt.  Nevertheless the canonical
row-product lifetime exceeds the best reverse row-product projection by a
linear amount.  The construction is still excluded by the actual-child
macroscopic-marginal theorem: each factor contains a positive-density
two-point common block.

All logarithms are natural.

## 1. The double-partial product group

Let

```math
m=k+\ell,\qquad n=r+s,\qquad m+n=N,
```

and suppose that all four ratios `k/N`, `ell/N`, `r/N`, and `s/N`
converge to positive constants.  Let

```math
X=(\sigma\mathbf1_k,\xi_1,\ldots,\xi_\ell),
\qquad
Y=(\tau\mathbf1_r,\eta_1,\ldots,\eta_s),             \tag{DP.1}
```

where all displayed signs are independent and fair, and put `Q=XY^T`.
The two factor supports are multiplicative groups.  Modulo the simultaneous
global flip, the rank-one support is a multiplicative subgroup `H` of the
bridge cube of size

```math
|H|=2^{\ell+s+1}.                                    \tag{DP.2}
```

Thus the uniform rank-one prior satisfies

```math
\|\mu\|_\infty=\sum_Q\mu(Q)^2=2^{-(\ell+s+1)},       \tag{DP.3}
```

and the two projective factor laws have maximum atoms `2^(-ell)` and
`2^(-s)`, respectively.  In particular, both factors have exponential
whole-word spread.  The standard binomial-tail argument used for the
partial-common-row construction also gives exponentially small sufficiently
narrow caps around every fixed factor word, and hence around every fixed
rank-one word.

This spread is not local.  On the first `k` coordinates, `X` is supported on
the two words `+-mathbf1_k`; on the first `r` coordinates, `Y` has the same
defect.  Therefore the example violates the uniform macroscopic conditional
and marginal spread of actual minimizing children.

## 2. Likelihood and the canonical row product

Let `u=beta/sqrt(N)`, let `U` be the fair law on the `mn` bridge bits, and
write the forward likelihood as

```math
p(B)={E_{X,Y}\exp\{u\langle B,XY^T\rangle\}
       \over(\cosh u)^{mn}}.                         \tag{DP.4}
```

For every row, irrespective of whether its index lies in the common or
nuisance part of `X`, the latent row word has the distribution

```math
(a\mathbf1_r,c_1,\ldots,c_s),
```

with `a,c_1,...,c_s` independent and fair.  If
`S_D(b)=sum_(j<=r)b_j`, the erased-row likelihood is therefore exactly

```math
p_i(b)={\cosh\{uS_D(b)\}\over(\cosh u)^r}.           \tag{DP.5}
```

For fixed `lambda>0`, let

```math
{dq\over dU}={p^{-\lambda}\over E_Up^{-\lambda}},
\qquad
{d\nu\over dU_r}
 ={\cosh(uS_D)^{-\lambda}\over Z_{r,N}},
\qquad
Z_{r,N}=E_{U_r}\cosh(uS_D)^{-\lambda}.              \tag{DP.6}
```

The canonical inverse row product is

```math
R=\bigotimes_{i=1}^m(\nu\otimes U_s).                \tag{DP.7}
```

Define the comparison row product

```math
P=U_n^{\otimes k}\otimes
  (\nu\otimes U_s)^{\otimes\ell}.                   \tag{DP.8}
```

Thus `P` removes the canonical retuning only on the positive-density common
row block.

## 3. A conditional one-dimensional identity

Let `C=[k]`, `D=[r]`, and

```math
T(B)=\sum_{i\in C,j\in D}B_{ij}.                    \tag{DP.9}
```

Fix every edge outside `C times D`.  Conditional on
`a=sigma tau`, all bridge energy inside that block is `aT`.  Summing the
remaining positive latent weights separately for `a=+-1` yields positive
numbers `W_+` and `W_-`, depending only on the fixed outside edges.  Hence

```math
\boxed{
\log p(B)=c(B_{(C\times D)^c})
 +\log\cosh\{uT(B)+h(B_{(C\times D)^c})\},}          \tag{DP.10}
```

where `h=(1/2)log(W_+/W_-)`; all constants independent of `T` have been
absorbed into `c`.  This identity is exact and does not require a
high-temperature expansion of the nuisance rectangle.

Under both `R` and `P`, the block `C times D` is independent of its
complement, and the complement has exactly the same law.  Under `P`, `T` is
a sum of `kr` fair signs.  Under `R`, it is a sum of `k` independent,
centered row sums with one-row law `nu`.  Since `z mapsto log cosh z` is
one-Lipschitz,

```math
\left|E_R\log p-E_P\log p\right|
\le u\{E_R|T|+E_P|T|\}.                             \tag{DP.11}
```

The elementary central-window bound

```math
Z_{r,N}\ge {1\over2}\cosh(\sqrt2\,\beta)^{-\lambda}
 =:z_*(\beta,\lambda)>0                             \tag{DP.12}
```

follows from `U_r{|S_D|<=sqrt(2r)}>=1/2` and `r<=N`.  Therefore

```math
E_\nu S_D^2\le {r\over Z_{r,N}}\le {r\over z_*},
```

and (DP.11) gives the uniform estimate

```math
\boxed{
\left|E_R\log p-E_P\log p\right|
\le {\beta\sqrt{kr}\over\sqrt N}
       (1+z_*^{-1/2})=O_{\beta,\lambda}(\sqrt N).}   \tag{DP.13}
```

The large nuisance rectangle can therefore change the comparison only at
the square-root scale.  There is no hidden linear cancellation.

## 4. Linear canonical-minus-best-product gap

Put

```math
d_N=D(\nu\Vert U_r)
=-\lambda E_\nu\log\cosh(uS_D)-\log Z_{r,N}.        \tag{DP.14}
```

If `r/N->rho>0`, the triangular-array CLT and bounded convergence give

```math
d_N\longrightarrow d_0
=D\!\left(
 {e^{-\lambda\log\cosh G}\over
  E e^{-\lambda\log\cosh G}}\,\gamma_{\beta^2\rho}
 \middle\Vert \gamma_{\beta^2\rho}
 \right)>0,                                        \tag{DP.15}
```

where `G~N(0,beta^2 rho)`.  Positivity holds because `beta,rho,lambda>0`
make the tilt nonconstant.

Let

```math
J=D(R\Vert q),
\qquad
I^\leftarrow=\inf_{A=\otimes_iA_i}D(A\Vert q).
```

Using `P` as a competitor and cancelling the common normalizer of `q`,

```math
\begin{aligned}
J-I^\leftarrow
&\ge D(R\Vert q)-D(P\Vert q)\\
&=D(R\Vert U)-D(P\Vert U)
 +\lambda(E_R-E_P)\log p\\
&=k d_N+\lambda(E_R-E_P)\log p.
\end{aligned}                                      \tag{DP.16}
```

Combining (DP.13)--(DP.16), if `k/N->kappa>0`, gives

```math
\boxed{
J-I^\leftarrow\ge \kappa d_0N-O(\sqrt N)+o(N)
=\Omega(N).}                                       \tag{DP.17}
```

Thus making *both* factors globally diffuse does not cure the coherent
row-retuning obstruction.  The exact conditional `log cosh` identity is
what prevents the two nuisance sides from producing a cancelling linear
term.

## 5. Averaged posterior invariance

The rank-one support `H` is a multiplicative subgroup.  For `g in H`,
entrywise switching `B mapsto gB` preserves `U`, satisfies `p(gB)=p(B)`,
and permutes the latent channel components transitively.  Consequently, for
every bridge law whose density relative to `U` is a function of `p(B)`, the
averaged Bayes posterior on `H` is exactly uniform:

```math
\boxed{\bar\mu=\mu.}                                \tag{DP.18}
```

Every deterministic latent quotient, including every energy-shell or
geometric split, consequently has zero averaged-posterior retuning, while
(DP.17) is linear.

## 6. Exact scope

The example proves the generic no-go

```math
\boxed{
\begin{gathered}
\text{exponential whole-word spread of both factors}\ +
\text{zero posterior retuning}\\
\centernot\Longrightarrow J-I^\leftarrow=o(N).
\end{gathered}}                                    \tag{DP.19}
```

It does **not** survive the actual-child macroscopic-marginal theorem:
the `k`-coordinate marginal of `X` and the `r`-coordinate marginal of `Y`
are each supported on two antipodal words.  Hence the next generic
falsifier must replace the hard common blocks by a diffuse synchronization
mechanism, rather than merely adding diffuse nuisance coordinates outside
them.
