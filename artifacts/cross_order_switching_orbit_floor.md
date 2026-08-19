# Switching-orbit bridge averaging has a low-temperature linear floor

Status: **proved direct cross-order reduction and a sharply scoped
method-class obstruction for actual own-scale pressure minimizers**.  The
reduction averages the parent *partition function* over all row and column
switches of one bridge.  Even after optimizing the bridge orbit, its
orientation-blind certificate has a positive linear defect at every fixed
`beta>4.05160...` on large balanced splits.  The exact orientation gain is
quantified separately.  Evading the floor through that scalar requires an
exponentially one-sided signed partition, but the best presently available
bound on its rate is large enough that this escape cannot be excluded.

This does not lower-bound the true parent optimum.  It rules out obtaining
an `o(N)` defect by replacing the rare-bridge minimum with a full
row/column-switch-orbit average and then retaining only pure bipartite
pressure and radial child payments.

## 1. Exact switching-orbit reduction

For a hollow sign matrix `A`, write

```math
Z_A(t)=\mathbb E_x\cosh(tH_A(x)),
\qquad
S_A(t)=\mathbb E_x\sinh(tH_A(x)),
\qquad
u_A(t)={S_A(t)\over Z_A(t)}.
```

For an `m` by `n` sign matrix `B`, define its pure bipartite pressure

```math
\psi_B(t)=\log\mathbb E_{p\in\{\pm1\}^m,
q\in\{\pm1\}^n}e^{t p^{\mathsf T}Bq},
\qquad
\Psi_{m,n}(t)=\min_B\psi_B(t).                       \tag{1.1}
```

Let `A,D` be arbitrary child signings and put

```math
L_\epsilon(C)=\log\mathbb E_{x,y}\cosh\left(
t\{H_A(x)+\epsilon H_D(y)+x^{\mathsf T}Cy\}\right).
```

For uniform independent switches `p,q`, set
`C_(p,q)=diag(p)B diag(q)`.  For fixed `x,y`, the random variable
`x^T C_(p,q)y` has the same law as `p^T Bq` and is symmetric.  Hence

```math
\begin{aligned}
\mathbb E_{p,q}e^{L_\epsilon(C_{p,q})}
&=e^{\psi_B(t)}\mathbb E_{x,y}
   \cosh\{t(H_A(x)+\epsilon H_D(y))\}\\
&=e^{\psi_B(t)}
  \{Z_A(t)Z_D(t)+\epsilon S_A(t)S_D(t)\}.           \tag{1.2}
\end{aligned}
```

Choose the relative orientation which makes the second term nonpositive,
then use `min log <= log average`.  Optimizing `B` proves

```math
\boxed{
\min_{\epsilon,C}L_\epsilon(C)
\le \phi_A(t)+\phi_D(t)+\Psi_{m,n}(t)
 +\log(1-|u_A(t)u_D(t)|).}                          \tag{1.3}
```

Dropping the last, nonpositive term gives the orientation-blind form

```math
\boxed{
\min_{\epsilon,C}L_\epsilon(C)
\le \phi_A(t)+\phi_D(t)+\Psi_{m,n}(t).}             \tag{1.4}
```

Now take `N=m+n`, `t=beta/sqrt(N)`,
`s_m=beta/sqrt(m)`, and `s_n=beta/sqrt(n)`.  Let `A,D` be exact
own-scale minimizers defining `P_m(beta),P_n(beta)`, and put

```math
\Delta_A=\phi_A(s_m)-\phi_A(t),
\qquad
\Delta_D=\phi_D(s_n)-\phi_D(t).
```

Parent minimization and (1.3) give the direct cross-order implication

```math
\boxed{
E_{m,n}(\beta)
\le \Psi_{m,n}(\beta/\sqrt N)-\Delta_A-\Delta_D
 +\log(1-|u_Au_D|).}                                \tag{1.5}
```

In particular,

```math
\boxed{
\Psi_{m,n}(\beta/\sqrt N)-\Delta_A-\Delta_D
\le C_\beta N^{1-\delta}
\quad\Longrightarrow\quad
E_{m,n}(\beta)\le C_\beta N^{1-\delta}.}           \tag{1.6}
```

Thus the pure bipartite quantity is not an unrelated proxy: (1.6) is its
immediate quantitative arrow to the desired defect.

## 2. Finite-order lower bound on the optimized orbit certificate

We now specialize to the balanced split `N=2r` and use the same exact
own-scale minimizer `A` on both sides.  Put

```math
s={\beta\over\sqrt r},
\qquad
t={\beta\over\sqrt{2r}},
\qquad
\Delta_A=\phi_A(s)-\phi_A(t).
```

Let `S_r=epsilon_1+...+epsilon_r` for independent fair signs and write
`mu_r=E|S_r|`.  For fixed `q`, exact integration over `p` gives

```math
\mathbb E_p e^{t p^{\mathsf T}Bq}
=\prod_{i=1}^r\cosh(t(Bq)_i).
```

Since `log cosh z>=|z|-log2`, Jensen's inequality and averaging uniform
`q` show

```math
\begin{aligned}
\psi_B(t)
&=\log\mathbb E_q\exp\left\{
  \sum_i\log\cosh(t(Bq)_i)\right\}\\
&\ge t\,\mathbb E_q\sum_i|(Bq)_i|-r\log2
=tr\mu_r-r\log2.                                    \tag{2.1}
\end{aligned}
```

The last equality holds because every row dot product has the law of
`S_r`.  Therefore

```math
\boxed{
\Psi_{r,r}(t)
\ge {\beta r\mu_r\over\sqrt{2r}}-r\log2.}           \tag{2.2}
```

For a rectangular `m` by `n` bridge, the same argument and its transpose
give the all-order bound

```math
\Psi_{m,n}(t)
\ge\max\{tm\mathbb E|S_n|-m\log2,
          tn\mathbb E|S_m|-n\log2\}.                \tag{2.2a}
```

Let `M(A)=max_x|H_A(x)|`.  Since the logarithmic derivative of
`phi_A` lies in `[0,M(A)]`,

```math
\Delta_A\le(s-t)M(A).                               \tag{2.3}
```

Write `U_r=M_r/r^(3/2)`.  Pressure minimality, evaluated at a ground-cap
minimizer, gives

```math
P_r(\beta)=\phi_A(s)\le sM_r=\beta U_r r.           \tag{2.4}
```

On the other hand the two maximizing configurations `x,-x` give

```math
\phi_A(s)\ge sM(A)-r\log2.
```

Consequently

```math
M(A)\le\left(U_r+{\log2\over\beta}\right)r^{3/2},
```

and (2.3) becomes

```math
\boxed{
\Delta_A
\le(1-2^{-1/2})(\beta U_r+\log2)r.}                 \tag{2.5}
```

Combining (2.2) and (2.5), the best orientation-blind full-orbit-average
certificate

```math
\mathcal C_r^{\rm orb}(\beta)
:=\Psi_{r,r}(\beta/\sqrt{2r})-2\Delta_A
```

satisfies the finite-order floor

```math
\boxed{
{\mathcal C_r^{\rm orb}(\beta)\over2r}
\ge {\beta\mu_r\over2\sqrt{2r}}
-{1\over2}\log2
-(1-2^{-1/2})(\beta U_r+\log2).}                    \tag{2.6}
```

Every term is explicit.  This is a lower bound on the *certificate value*,
not on `E_(r,r)`.

## 3. The asymptotic threshold

The central limit theorem gives

```math
{\mu_r\over\sqrt r}\longrightarrow\sqrt{2/\pi},
```

and the rigorous all-order cap construction gives
`limsup U_r<=1/2`.  Hence

```math
\liminf_{r\to\infty}{\mathcal C_r^{\rm orb}(\beta)\over2r}
\ge a_*\beta-b_*,                                   \tag{3.1}
```

where

```math
a_*={1\over2\sqrt\pi}-{1\over2}(1-2^{-1/2})
   =0.135648182367151...,
```

and

```math
b_*=(3/2-2^{-1/2})\log2
   =0.549591699105644....                            \tag{3.2}
```

Thus the critical value for this method is

```math
\boxed{
\beta_{\rm orb}={b_*\over a_*}
=4.051596486697427....}                              \tag{3.3}
```

For every fixed `beta>beta_orb`, there are `c_beta>0` and `r_beta` such
that for all `r>=r_beta`,

```math
\boxed{
\mathcal C_r^{\rm orb}(\beta)\ge c_\beta(2r).}       \tag{3.4}
```

Therefore pure bipartite pressure plus the two radial child payments,
when obtained only through full switching-orbit averaging, cannot certify
an `o(N)` cross-order defect in the low-temperature range needed for the
ground-state limit.

## 4. Exact orientation term and the unresolved escape

The exact oriented orbit certificate from (1.5) is

```math
\mathcal C_r^{\rm orb,or}(\beta)
=\mathcal C_r^{\rm orb}(\beta)+\log(1-|u_A(t)|^2).
                                                                    \tag{4.1}
```

Thus (3.4) gives a quantitative necessary condition.  If, for fixed
`beta>beta_orb`, this oriented certificate is `o(r)`, then

```math
-\log(1-|u_A(t)|^2)\ge c_\beta(2r)-o(r),             \tag{4.2}
```

or equivalently

```math
\boxed{
1-|u_A(t)|^2\le\exp\{-c_\beta(2r)+o(r)\}.}          \tag{4.3}
```

The rate in (4.3) is not unconstrained.  Since `E H_A=0`, Jensen gives

```math
Z_A(t)+S_A(t)=\mathbb E e^{tH_A}\ge1,
\qquad
Z_A(t)-S_A(t)=\mathbb E e^{-tH_A}\ge1.
```

Consequently

```math
1-|u_A(t)|\ge e^{-\phi_A(t)}.
```

For arbitrary children this yields the exact lower bound

```math
\boxed{
1-|u_A(t)u_D(t)|
\ge\exp\{-\min(\phi_A(t),\phi_D(t))\}.}             \tag{4.4}
```

There is a sharper way to propagate (4.4) in the balanced same-child case.
Since `Delta_A=P_r(beta)-phi_A(t)`, equations (4.1) and (4.4) give

```math
\mathcal C_r^{\rm orb,or}
\ge\Psi_{r,r}(t)-2P_r(\beta)+\phi_A(t).             \tag{4.5}
```

Put `q=t/s=2^(-1/2)`.  Represent `Z_A` as a uniform exponential sum on the
`2^r` oriented projective states `(tau,[x])`.  For nonnegative numbers
`a_omega` and `0<q<1`,

```math
\sum_\omega a_\omega^q
\ge\left(\sum_\omega a_\omega\right)^q.
```

Applied with `a_omega=exp(s tau H_A(x))`, this finite-support interpolation
gives

```math
\phi_A(t)\ge qP_r(\beta)-(1-q)r\log2.               \tag{4.6}
```

Combining (2.2), (2.4), (4.5), and (4.6) gives the fully quantified lower
bound

```math
\boxed{
{\mathcal C_r^{\rm orb,or}(\beta)\over2r}
\ge {\beta\mu_r\over2\sqrt{2r}}
-{2-2^{-1/2}\over2}(\beta U_r+\log2).}             \tag{4.7}
```

At the presently known endpoint `U_r<=1/2+o(1)`, the leading `beta`
coefficient in (4.7) tends to

```math
{1\over2\sqrt\pi}-{2-2^{-1/2}\over4}
=-0.041128512929485... .                            \tag{4.8}
```

Therefore existing one-sided information does **not** extend the positive
linear floor to the exact oriented certificate.  The conservative
conclusion is:

- orientation-blind full-orbit averaging is rigorously blocked above
  `beta_orb`;
- an oriented proof must establish the exponential imbalance (4.3), with
  its possible rate bounded by (4.4)--(4.5); and
- a proof which selects one exceptionally aligned switch rather than using
  the orbit average lies outside the obstruction altogether.

## 5. Scope

The result applies to actual own-scale pressure-minimizing children, all
finite balanced orders in (2.6), and every sufficiently large balanced
order in (3.4).  It establishes no lower bound on the globally optimized
parent pressure and does not exclude a genuinely correlated rare bridge.
It does exclude a broad and natural replacement of that rare minimum:

```math
\text{rare bridge}
\longmapsto
\text{full row/column-switch average}
\longmapsto
\text{pure bipartite pressure plus radial payments}.
```

At fixed low temperature the orientation-blind replacement has a positive
linear defect.  Exact orientation can evade the proved floor only through
the quantified one-sided mechanism in Section 4; present bounds neither
force nor exclude it.
