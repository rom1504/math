# Independent verification of the common-sign row-lifetime no-go

Source audited:
[`generic_rank_one_common_sign_row_lifetime_no_go.md`](generic_rank_one_common_sign_row_lifetime_no_go.md).

Verdict: **pass**.  The conditional Rényi-two estimate, triangular-array
limits, canonical-error asymptotic, sublinear best-product upper bound, and
bounded-time lifetime witness are all valid with the stated normalizations.
The construction is a generic two-word rank-one channel, not an
actual-minimizer construction; the source states this limitation correctly.

## 1. Channel and canonical factors

For `Q=sigma 1_m 1_n^T`,

```math
E_\sigma\prod_{ij}(1+\rho\sigma B_{ij})
={\cosh(uS)\over(\cosh u)^{mn}}
```

follows from `1+rho*sigma*b=exp(u*sigma*b)/cosh(u)`.  Averaging all rows
except row `i` multiplies the numerator by `(cosh u)^((m-1)n)`, so the exact
row marginal is

```math
p_i(B_i)={\cosh(uS_i)\over(\cosh u)^n}.
```

Thus both normalized inverse escorts in (RN.7) are correct.  In particular,
the `r_i` used later is exactly the canonical inverse factor, not a chosen
surrogate.

## 2. Conditional Rényi-two calculation

Let `f,g` be normalized positive densities on the fair `n`-cube and suppose
each of `log f,log g` changes by at most `2a` when one bit is changed.  With

```math
V={na^2\over2},
```

bounded differences gives

```math
\log E e^{s(\log f-E\log f)}\le s^2V.             \tag{V.RN.1}
```

Normalization and Jensen imply

```math
-V\le E\log f\le0,                                \tag{V.RN.2}
```

because (V.RN.1) at `s=1` gives
`0=log E f<=E log f+V`.  Hölder with exponents `3/2,3` then gives

```math
\begin{aligned}
\log E{f^2\over g}
&\le {2\over3}\log Ef^3+{1\over3}\log Eg^{-3}\\
&\le {2\over3}(9V)+{1\over3}(12V)=10V=5na^2.
                                                               \tag{V.RN.3}
\end{aligned}
```

The same argument with `f,g` exchanged proves the reverse direction.  For
the conditional escort and canonical row factor, `log cosh` being
one-Lipschitz makes the one-bit oscillation exactly at most `2lambda*u`,
uniformly in the external row sum `C`.  Hence `a=lambda*u` in (V.RN.3), and

```math
\sup_C\max\{D_2(q_i^C\Vert r_i),D_2(r_i\Vert q_i^C)\}
\le5\lambda^2u^2n=O_{\beta,\lambda}(1).           \tag{V.RN.4}
```

No averaging over `C` or exceptional-conditioning assumption is hidden in
this step.

## 3. CLT and uniform-integrability checks

Since `u=beta/sqrt(N)` and `n/N->1-theta`, the fair row field satisfies

```math
uS_n\Longrightarrow G,
\qquad \operatorname {Var}G=\beta^2(1-\theta).    \tag{V.RN.5}
```

Both functions

```math
x\longmapsto\cosh(x)^{-\lambda},
\qquad
x\longmapsto\log\cosh(x)\cosh(x)^{-\lambda}
```

are bounded and continuous for fixed `lambda>0`.  Therefore ordinary weak
convergence, without a further tail argument, proves (RN.15) and the
numerator limit in (RN.16).  The limiting density
`z_0^(-1)cosh(G)^(-lambda)` is nonconstant because
`beta,lambda>0` and `theta in (0,1)`.  Its KL divergence from the Gaussian
base is consequently

```math
d_0=-\lambda {E[f(G)e^{-\lambda f(G)}]\over z_0}
       -\log z_0>0.                                \tag{V.RN.6}
```

For the only growing observable used later,

```math
E_{r_i}S_i^2
={E_U[S_i^2\cosh(uS_i)^{-\lambda}]\over Z_{n,N}}
\le {n\over Z_{n,N}}=O(n),                        \tag{V.RN.7}
```

because (RN.15) bounds `Z_(n,N)` away from zero.  Symmetry gives
`E_(r_i)S_i=0`.  These are sufficient second-moment bounds; no unstated
uniform integrability of `S_i^2` is required.

## 4. Canonical and optimal-product scales

The exact density ratio is

```math
\log{dr\over dq}
=-\lambda\sum_i f(uS_i)-m\log Z_{n,N}
 +\lambda f\!\left(u\sum_iS_i\right)+\log Z_{d,N}.
```

Taking `r` expectation proves (RN.13).  Equations (V.RN.6)--(V.RN.7), row
independence, and `f(x)<=|x|` give

```math
mD(r_i\Vert U_n)=\theta d_0N+o(N),
\qquad
E_rf\!\left(u\sum_iS_i\right)=O(\sqrt N).         \tag{V.RN.8}
```

The second statement follows from

```math
uE_r\left|\sum_iS_i\right|
\le u\sqrt{mE_{r_i}S_i^2}=O(\sqrt N).
```

Also, `Z_(d,N)<=1`.  Chebyshev gives
`P(|S_d|<=sqrt(2d))>=1/2`, and hence

```math
-\log2-\lambda\log\cosh(u\sqrt{2d})
\le\log Z_{d,N}\le0.                              \tag{V.RN.9}
```

The left side is `-O(sqrt(N))` on comparable splits.  Substitution in the
exact ratio therefore proves

```math
\mathcal J=D(r\Vert q)=\theta d_0N+o(N).           \tag{V.RN.10}
```

For the best reverse product, the fair law is an admissible row product and

```math
\begin{aligned}
\mathcal I^{\leftarrow}
&\le D(U_d\Vert q)\\
&=\lambda E_Uf(uS_d)+\log Z_{d,N}\\
&\le\lambda uE|S_d|
\le\lambda u\sqrt d=O(\sqrt N).                  \tag{V.RN.11}
\end{aligned}
```

Dropping `log Z_(d,N)` is valid because it is nonpositive.  Thus
`J-I^leftarrow=theta*d_0*N+o(N)` as claimed.

## 5. Bounded-time witness

For `h=dq/dr`, the row-refresh semigroup contains the all-refreshed term
with weight `(1-e^(-T))^m`; positivity of every other retained-set term and
`E_rh=1` imply

```math
h_T\ge(1-e^{-T})^m.                                \tag{V.RN.12}
```

Applying the reverse de Bruijn identity from time `T` onward gives

```math
\int_T^\infty\mathfrak K_r(h_t)dt
=-E_r\log h_T
\le-m\log(1-e^{-T}).                              \tag{V.RN.13}
```

Let `eta_0=theta*d_0/2`.  Equation (V.RN.10) gives
`J>=eta_0*N` for all large `N`.  With

```math
T_0=-\log(1-e^{-\eta_0/2}),
```

the right side of (V.RN.13) is at most
`m*eta_0/2<=N*eta_0/2`.  Consequently

```math
\int_0^{T_0}\mathfrak K_r(h_t)dt
\ge {\eta_0N\over2},                              \tag{V.RN.14}
```

exactly as in (RN.27).  The witness time depends on fixed
`beta,lambda,theta` through `d_0`, but not on `N`.

Finally, the claimed time-zero formula is also exact.  For fixed `B_(-i)`,

```math
h={q_{-i}\over r_{-i}}{q_i(\,·\mid B_{-i})\over r_i}.
```

Thus `(E_i h)(E_i h^(-1))-1` is precisely
`chi^2(r_i||q_i(.|B_(-i)))`.  The direction agrees with the second bound in
(V.RN.4).

## 6. Scope judgment

The source proves the intended generic no-go: exact two-word rank-one
support and uniformly bounded conditional row Rényi-two complexity do not
force the canonical row-lifetime integral, or even a fixed bounded-time
part of it, to be sublinear.  It simultaneously exhibits an explicit fair
product with only `O(sqrt(N))` reverse error, so the linear canonical cost is
coherent factor retuning rather than a lower bound on irreducible reverse
dependence.  No conclusion about actual contracted-temperature minimizing
children follows without an additional optimizer-specific theorem, and the
source does not claim one.
