# A two-word rank-one obstruction to canonical row-lifetime closure

Status: **rigorous scalable generic-channel no-go**. This note constructs an
exact rank-one binary-channel likelihood whose canonical inverse row product
has linear reverse error, even though the best reverse row-product error is
sublinear. Every conditional row remains in a common bounded Renyi-two
neighborhood, and the row-refresh lifetime integral is already linear on a
fixed bounded time interval.

The construction is **not** a sequence of contracted-temperature minimizing
children. Its precise role is to prove that bounded conditional row
complexity and the exact rank-one latent representation do not, by
themselves, close the actual-child row-lifetime SML.

## 1. Exact two-word rank-one channel

Let

```math
m+n=N,\qquad {m\over N}\longrightarrow\theta\in(0,1),\qquad
u={\beta\over\sqrt N},\qquad \rho=\tanh u,
\tag{RN.1}
```

where `beta,lambda>0` are fixed. On the `d=mn` bridge bits let `U_d` be
fair, and let the latent word be

```math
Q=\sigma\mathbf 1_m\mathbf 1_n^{\mathsf T},\qquad
\Pr\{\sigma=1\}=\Pr\{\sigma=-1\}={1\over2}.
\tag{RN.2}
```

Thus the support consists of two rank-one sign matrices. Pass every entry
through the binary channel of amplitude `u`. If

```math
S_i(B)=\sum_{j=1}^nB_{ij},\qquad S(B)=\sum_iS_i(B),
\tag{RN.3}
```

then the output likelihood relative to `U_d` is exactly

```math
p(B)
={1\over2}\sum_{\sigma=\pm1}\prod_{ij}(1+\rho\sigma B_{ij})
={\cosh(uS(B))\over(\cosh u)^d}.                  \tag{RN.4}
```

Indeed, `1+sigma rho b=exp(u sigma b)/cosh u`. Its erased-row likelihood is

```math
p_i(b)={\cosh(uS_i(b))\over(\cosh u)^n}.          \tag{RN.5}
```

For a sum `S_k` of `k` fair signs, put

```math
Z_{k,N}=E_{U_k}\cosh(uS_k)^{-\lambda}.             \tag{RN.6}
```

The full inverse escort and canonical inverse row product are therefore

```math
{dq\over dU_d}={\cosh(uS)^{-\lambda}\over Z_{d,N}},
\qquad
r=\bigotimes_{i=1}^m r_i,\qquad
{dr_i\over dU_n}={\cosh(uS_i)^{-\lambda}\over Z_{n,N}}.
\tag{RN.7}
```

The factors `(cosh u)^k` cancel on normalization.

## 2. Uniform conditional row complexity

Fix the other rows and write `C=sum_(k ne i)S_k`. The conditional law of
row `i` under `q` is

```math
{dq_i^C\over dU_n}(b)
={\cosh\{u(C+S_i(b))\}^{-\lambda}
  \over E_{U_n}\cosh\{u(C+S_i)\}^{-\lambda}}.     \tag{RN.8}
```

Flipping one bit changes the argument by `2u`, while `log cosh` is
one-Lipschitz. Hence both `log(dq_i^C/dU_n)` and
`log(dr_i/dU_n)` have one-bit oscillation at most `2lambda u`, uniformly
in `C`.

We use the following elementary comparison. If `f,g` are normalized
positive densities on the fair `n`-cube and both logarithms have one-bit
oscillation at most `2a`, then

```math
D_2(fU_n\Vert gU_n),\ D_2(gU_n\Vert fU_n)\le5na^2. \tag{RN.9}
```

To verify it, put `V=na^2/2`. Bounded differences gives

```math
\log E e^{s(\log f-E\log f)}\le s^2V,
\qquad -V\le E\log f\le0,
\tag{RN.10}
```

and similarly for `g`. Holder with exponents `3/2,3` gives

```math
\log E{f^2\over g}
\le {2\over3}\log Ef^3+{1\over3}\log Eg^{-3}
\le {2\over3}(9V)+{1\over3}(12V)=10V.
\tag{RN.11}
```

Swap `f,g` for the reverse direction. With `a=lambda u` this proves

```math
\boxed{
\sup_C\max\{D_2(q_i^C\Vert r_i),D_2(r_i\Vert q_i^C)\}
\le5\lambda^2u^2n=O_{\beta,\lambda}(1).}          \tag{RN.12}
```

Thus neither escaping row complexity nor rare exceptional conditioning
causes the obstruction.

## 3. Linear canonical error and sublinear best-product error

Let `f(x)=log cosh x`. Direct substitution of (RN.7) gives

```math
\boxed{
\mathcal J:=D(r\Vert q)
=mD(r_i\Vert U_n)
+\lambda E_r f\!\left(u\sum_iS_i\right)
+\log Z_{d,N}.}                                   \tag{RN.13}
```

We estimate every term. Under `U_n`, the triangular-array CLT gives

```math
uS_n\Longrightarrow G,\qquad
G\sim N(0,\sigma^2),\qquad
\sigma^2=\beta^2(1-\theta).                       \tag{RN.14}
```

Both `exp(-lambda f)` and `f exp(-lambda f)` are bounded and continuous.
Consequently

```math
Z_{n,N}\longrightarrow z_0:=Ee^{-\lambda f(G)}\in(0,1)
\tag{RN.15}
```

and

```math
D(r_i\Vert U_n)\longrightarrow
d_0:=-\lambda {E[f(G)e^{-\lambda f(G)}]\over z_0}-\log z_0>0.
\tag{RN.16}
```

The limit is the KL divergence of the nonconstant tilted law with density
`z_0^(-1)e^(-lambda f)` relative to the Gaussian law, hence is strictly
positive.
Equation (RN.15) also bounds `Z_(n,N)` away from zero. By symmetry,

```math
E_{r_i}S_i=0,\qquad
E_{r_i}S_i^2
={E_{U_n}[S_i^2\cosh(uS_i)^{-\lambda}]\over Z_{n,N}}
\le {n\over Z_{n,N}}=O(n).                        \tag{RN.17}
```

Using `f(x)<=|x|` and row independence under `r`,

```math
E_r f\!\left(u\sum_iS_i\right)
\le u\{mE_{r_i}S_i^2\}^{1/2}=O(\sqrt N).          \tag{RN.18}
```

Also `Z_(d,N)<=1`. Chebyshev gives
`U_d{|S_d|<=sqrt(2d)}>=1/2`, so

```math
Z_{d,N}\ge {1\over2}\cosh(u\sqrt{2d})^{-\lambda},
\qquad -O(\sqrt N)\le\log Z_{d,N}\le0.            \tag{RN.19}
```

Equations (RN.13)--(RN.19) imply

```math
\boxed{\mathcal J=\theta d_0N+o(N)=\Theta(N).}     \tag{RN.20}
```

In contrast, `U_d` is itself a row product. Therefore

```math
\begin{aligned}
\mathcal I^{\leftarrow}
&:=\inf_{a=\otimes_i a_i}D(a\Vert q)\le D(U_d\Vert q)\\
&=\lambda E_{U_d}f(uS_d)+\log Z_{d,N}
\le\lambda uE|S_d|
\le\lambda u\sqrt d=O(\sqrt N).
\end{aligned}                                      \tag{RN.21}
```

Hence

```math
\boxed{\mathcal J-\mathcal I^{\leftarrow}
=\theta d_0N+o(N).}                                \tag{RN.22}
```

The obstruction is coherent product retuning, not irreducible extensive
reverse dependence. In this example the explicit product `U_d` exposes
that fact, but selecting an analogous product for a general actual-child
law is the unresolved retuning problem.

## 4. The bounded-time lifetime is already linear

Let

```math
h={dq\over dr},\qquad h_t=P_t^rh,\qquad
\mathfrak K_r(g)=\sum_iE_{r_{-i}}[(E_ig)(E_ig^{-1})-1].
\tag{RN.23}
```

The reverse row-noise identity gives

```math
\mathcal J=\int_0^\infty\mathfrak K_r(h_t)\,dt.    \tag{RN.24}
```

Its tail has a dimension-explicit bound. At time `T`, each row is retained
with probability `p=e^(-T)`. The term in which no row is retained has
weight `(1-p)^m` and conditional density `E_rh=1`; all other terms are
nonnegative. Thus

```math
h_T\ge(1-e^{-T})^m,\qquad
\int_T^\infty\mathfrak K_r(h_t)dt
=-E_r\log h_T\le-m\log(1-e^{-T}).                 \tag{RN.25}
```

Put `eta_0=theta d_0/2`. By (RN.20), `J>=eta_0N` for all large `N`. Choose
the fixed time

```math
T_0=-\log(1-e^{-\eta_0/2}).                        \tag{RN.26}
```

Then `-log(1-e^(-T_0))=eta_0/2`, and `m<=N` gives

```math
\boxed{
\int_0^{T_0}\mathfrak K_r(h_t)dt
\ge {\eta_0\over2}N.}                             \tag{RN.27}
```

Thus bounded row `D_2` and two-word rank-one support do not even make the
fixed bounded-time part of the row-lifetime integral sublinear.

For reference, the time-zero integrand has the exact form

```math
\mathfrak K_r(h)
=\sum_iE_{r_{-i}}\chi^2
\{r_i\Vert q(R_i\mid B_{-i})\}.                   \tag{RN.28}
```

This follows by writing `E_i h=q_{-i}/r_{-i}` in the pair formula. Equation
(RN.12) bounds every summand by a parameter-dependent constant. No lower
bound on (RN.28) is needed for the fixed-window obstruction (RN.27).

## 5. Consequence for the actual-child SML

The example proves that exact rank-one latent support plus uniformly
bounded conditional row Renyi-two complexity does not imply `J=o(N)`,
does not make the early row-noise lifetime sublinear, and does not separate
`J` from `J-I^leftarrow`.

Running the row-noise identity relative to the unknown optimal product
`p^*` would represent `I^leftarrow` exactly, but would restore the full
product oracle. Relative to the explicit canonical product, the lifetime
necessarily sees both genuine cross-row dependence and coherent factor
retuning.

Therefore a positive actual-child theorem must use a genuine minimizer
identity to suppress the common-sign/shared-latent retuning mechanism, or
must construct a different explicit product from coarser child data and
prove that its excess over `I^leftarrow` is sublinear. Even that would
resolve only the product phase; target reach remains a separate recurrence
obligation.
