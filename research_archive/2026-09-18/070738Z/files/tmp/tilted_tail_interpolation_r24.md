### Wave 24, Route 2: selector collapse and the full parent interpolation

**Status.** The variational identities, normalizations, chain rules, and
finite calculations below are **Verified** by
`tmp/tilted_tail_interpolation_r24.py`. They sharpen the exact overlap target
and give a single full-interpolation formulation. They do not prove its
remaining endpoint tail or transport estimate.

#### 1. Hard overlap: the selector optimization collapses exactly

Fix $A,m,\beta$ and tolerance $a$. Let $U=U_m$, let $\nu=\nu_\beta$, and
write

```math
G=\{(S,d):Q(A[S])-c_A(S,d)\le a\},
\qquad
\Omega_S=\nu(G_S),
\qquad
\overline\Omega=(U\otimes\nu)(G)=\mathbb E_U\Omega_S.
```

If $\overline\Omega>0$, define

```math
\pi_*(S)={U(S)\Omega_S\over\overline\Omega}.
\tag{R24.1}
```

It assigns zero mass to every selector with $\Omega_S=0$. For any selector
law $\pi$ supported on $\{\Omega_S>0\}$,

```math
\begin{aligned}
\mathbb E_\pi[-\log\Omega_S]+D_{\rm KL}(\pi\Vert U)
&=\sum_S\pi(S)\log{\pi(S)\over U(S)\Omega_S}\\
&=D_{\rm KL}(\pi\Vert\pi_*)-\log\overline\Omega.
\end{aligned}
\tag{R24.2}
```

Hence

```math
\boxed{
\inf_\pi\left\{
\mathbb E_\pi[-\log\Omega_S]+D_{\rm KL}(\pi\Vert U)
\right\}
=-\log\overline\Omega.
}
\tag{R24.3}
```

If $\overline\Omega=0$, both sides are $+\infty$ under the usual extended
KL conventions, so (R24.3) still holds.

There is an equivalent one-line chain-rule proof. Put
$R=U\otimes\nu$ on selector--cut pairs and $P=R(\cdot\mid G)$. Then

```math
D_{\rm KL}(P\Vert R)=-\log R(G)=-\log\overline\Omega.
```

The selector marginal of $P$ is (R24.1), and for every selector in its
support, $P(d\mid S)=\nu(d\mid G_S)$. Therefore

```math
\boxed{
-\log\overline\Omega
=D_{\rm KL}(\pi_*\Vert U)
+\mathbb E_{S\sim\pi_*}[-\log\Omega_S].
}
\tag{R24.4}
```

The conditional channel has distortion at most $a$. Consequently the exact
hard sufficient lemma replacing (10.727) is simply

```math
\boxed{
\Pr_{S\sim U_m,\,D\sim\nu_\beta}
\{Q(A[S])-c_A(S,D)\le a\}
\ge e^{-O(n^{1/2-2c})},
\quad a=O(n^{3/2-c}).
}
\tag{R24.5}
```

This is equivalent to the optimized selector statement for the fixed
$A,\beta,a$; the selector need not be guessed separately.

#### 2. Combination with the soft row-square interface

The same identity materially relaxes the overlap exponent when combined with
the output-weighted Hoeffding/reference theorem of
`tmp/soft_row_entropy_r24.md`. Let $\nu$ now be any law supported on exact
parent grounds and keep the same event $G$. Write
$Z=(U\otimes\nu)(G)$ and condition $P=(U\otimes\nu)(\cdot\mid G)$. The full
product-reference chain rule is

```math
\boxed{
-\log Z
=I_P(S;D)+D_{\rm KL}(P_S\Vert U)
+D_{\rm KL}(P_D\Vert\nu).
}
\tag{R24.6}
```

In particular the information plus selector cost is at most $-\log Z$, and

```math
\mathbb E_P R_2(D)
={\mathbb E_{U\otimes\nu}[R_2(D)\mathbf1_G]\over Z}.
\tag{R24.7}
```

Thus the newer downstream bound proves a power-saving restriction edge from
the rigorously sufficient joint conditions

```math
\boxed{
\begin{aligned}
a&=O(n^{3/2-c}),\\
-\log Z&=O(n^{3/4-c}),\\
\mathbb E_{U\otimes\nu}[R_2(D)\mid G]
&=O(n^{9/4-c}).
\end{aligned}}
\tag{R24.8}
```

Indeed, with $\lambda\asymp n^{-3/4}$, the information contribution is
$n^{3/4}\,O(n^{3/4-c})=O(n^{3/2-c})$, and the row term is
$\lambda O(n^{9/4-c})=O(n^{3/2-c})$. This is the exact factor check behind
the relaxed overlap exponent.

The row condition is not automatic. A ground only satisfies
$R_2\le(n-1)q_n=O(n^{5/2})$, losing $n^{1/4+c}$, and conditioning can bias
toward less regular grounds. For $A_9$ with the uniform ground prior, the
unconditional mean is $104.32$, while conditioning on exact child-ground
compatibility raises it to $107.428571\ldots$. Thus (R24.7) is analyzable as
one joint moment, but scalar overlap alone does not bound it.

#### 3. Soft channels: exact selector optimization

For the two-temperature construction, set

```math
\ell_S(y)=F_S(y)+(\beta-\gamma)c_S(y),
\qquad
J_S=\log\mathbb E_{\mu_{\gamma,S}}e^{\ell_S}
-\mathbb E_{\mu_{\gamma,S}}\ell_S.
```

Here $J_S=D_{\rm KL}(\mu_{\gamma,S}\Vert\nu_{\beta,S})$ exactly. The same
Gibbs variational calculation gives

```math
\boxed{
\inf_\pi\{\mathbb E_\pi J_S+D_{\rm KL}(\pi\Vert U)\}
=-\log\mathbb E_{S\sim U}e^{-J_S},
}
\tag{R24.9}
```

attained by
$\pi_J(S)=U(S)e^{-J_S}/\mathbb E_Ue^{-J_S}$. At finite temperatures every
$J_S$ is finite. The extended-value statement remains valid if some costs
are infinite, provided the normalizing expectation is interpreted in
$[0,1]$.

Using the full interpolation (10.751) for every selector, (R24.9) becomes
the exact optimized soft target

```math
\boxed{
\mathbb E_{S\sim U}
\exp\left\{-\int_0^1(1-t)
\operatorname{Var}_{\mu_{S,t}}(\ell_S)\,dt\right\}
\ge e^{-O(n^{1/2-2c})}.
}
\tag{R24.10}
```

Unlike a bound at $t=0$, this retains the entire rare-tail transition.

#### 4. One normalized interpolation with a uniform-parent endpoint

There is also a single joint interpolation, rather than one path per
selector. Put

```math
Z_0=\mathbb E_{S\sim U}Z_S(\gamma),
\qquad
\pi_0(S)={U(S)Z_S(\gamma)\over Z_0},
```

and on the disjoint union of pairs $(S,y)$ define

```math
\begin{aligned}
Z_t&=\mathbb E_{S\sim U}
\sum_y
e^{[(1-t)\gamma+t\beta]c_S(y)}K_{\beta,S}(y)^t,\\
M_t(S,y)&={U(S)
e^{[(1-t)\gamma+t\beta]c_S(y)}K_{\beta,S}(y)^t
\over Z_t}.
\end{aligned}
\tag{R24.11}
```

Every normalization is explicit:

```math
M_0(S,y)=\pi_0(S)\mu_{\gamma,S}(y),
\qquad
Z_1=Z_A(\beta),
\qquad
M_1(S,y)=U(S)\nu_{\beta,S}(y).
\tag{R24.12}
```

In particular, the endpoint selector is exactly uniform and independent of
a parent Gibbs cut. Since $\partial_t\log$ of the unnormalized weight is
$\ell_S(y)$,

```math
\boxed{
\begin{aligned}
D_{\rm KL}(M_0\Vert M_1)
&=D_{\rm KL}(\pi_0\Vert U)
+\mathbb E_{S\sim\pi_0}J_S\\
&=\log{Z_A(\beta)\over Z_0}-\mathbb E_{M_0}\ell_S\\
&=\int_0^1(1-t)\operatorname{Var}_{M_t}(\ell_S)\,dt.
\end{aligned}}
\tag{R24.13}
```

The first line is KL chain rule, and the last follows by twice integrating
$\partial_t^2\log Z_t=\operatorname{Var}_{M_t}(\ell_S)$. The child-Gibbs
distortion bound holds for every $S$, hence also under $\pi_0$. Therefore

```math
\boxed{
\int_0^1(1-t)\operatorname{Var}_{M_t}(\ell_S)\,dt
=O(n^{1/2-2c})
}
\tag{R24.14}
```

is a precise single-path sufficient lemma for the original parent-Gibbs
interface.

#### 5. Exact endpoint decomposition and the missing transport estimate

Lift $M_t$ to full cuts using the unchanged parent conditional law
$P_{S,y}=\nu_\beta(\cdot\mid D[S]=y)$. At $t=1$ the lifted endpoint is the
product law

```math
\widehat M_1(S,d)=U(S)\nu_\beta(d).
```

Set $C=Z_A(\beta)/Z_0$. Direct substitution, including all normalizers,
gives

```math
f(S,d):={d\widehat M_0\over d\widehat M_1}(S,d)
=C e^{-\ell_S(d[S])}.
\tag{R24.15}
```

Define the endpoint functions

```math
a(d)=\mathbb E_{S\sim U}e^{-\ell_S(d[S])},
\qquad
\rho_d(S)={U(S)e^{-\ell_S(d[S])}\over a(d)}.
\tag{R24.16}
```

Then $\mathbb E_{\nu_\beta}Ca=1$, the cut marginal of
$\widehat M_0$ is $Ca\,\nu_\beta$, and its conditional selector law given
$d$ is $\rho_d$. Entropy chain rule at the uniform-parent endpoint gives the
exact decomposition

```math
\boxed{
D_{\rm KL}(M_0\Vert M_1)
=\operatorname{Ent}_{\nu_\beta}(Ca)
+\mathbb E_{d\sim Ca\nu_\beta}
D_{\rm KL}(\rho_d\Vert U).
}
\tag{R24.17}
```

This states precisely what endpoint transport must control. Applying the
Johnson modified log-Sobolev inequality to $\rho_d$ for each fixed $d$
shows that the second term is at most

```math
\mathcal D_{\rm sel}
=\tau_{\rm mls}{m(n-m)\over2n}\,C
\mathbb E_{d\sim\nu_\beta,\,S\sim S'}
\left[
(e^{-\ell_{S'}}-e^{-\ell_S})(\ell_S-\ell_{S'})
\right].
\tag{R24.18}
```

Here the adjacent selectors use the same full cut $d$, and the integrand is
nonnegative. Thus the following is an exact endpoint sufficient lemma:

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}(Ca)
+\mathcal D_{\rm sel}
=O(n^{1/2-2c}).
}
\tag{R24.19}
```

Equation (R24.19) implies (R24.14) through (R24.17). It separates the two
missing pieces: transport of the cut marginal under the parent Gibbs law,
and exponential—not merely quadratic—selector-gradient control.

The endpoint estimates (10.746)--(10.750) do not prove (R24.19).
Equation (10.750) controls an unweighted square of
$\ell_{S'}-\ell_S$ (at best $O(n^{1/2+2c})$ before the Johnson factor),
whereas (R24.18) contains the rare-tail weight $e^{-\ell}$. Equation
(10.746) controls changes of the parent conditional law under the parent
marginal; it does not bound $\operatorname{Ent}_{\nu_\beta}(Ca)$. A new
exact-minimizer transport theorem would have to control these two endpoint
functionals at the displayed target exponent. Abstract product measures
already rule out a dimension-free conversion from an average one-coordinate
conditional KL to total relative entropy, so such a theorem cannot be
generic.

#### 6. Mean truncation is exactly circular

Under the endpoint product law $U\otimes\nu_\beta$, let
$\Delta=Q(A[S])-c_A(S,D)$. Since
$\mathbb E_Uc_A(S,d)=p_2\langle A,d\rangle$,

```math
\boxed{
\mathbb E\Delta
=\mathbb E_UQ(A[S])-p_2\mathbb E_{\nu_\beta}\langle A,D\rangle.
}
\tag{R24.20}
```

The thermal deficit is harmless:

```math
q_n-\mathbb E_{\nu_\beta}\langle A,D\rangle
\le {n\log2\over\beta}
=O(n^{3/2-c}).
\tag{R24.21}
```

But after adding and subtracting $p_2q_n$, the other part of (R24.20) is
exactly the unknown uniform optimized-restriction excess. Consequently a
Markov/truncated-mean proof of (R24.5) assumes the desired restriction bound;
it is not an independent overlap argument.

Temperature interpolation of the hard event is also exact. If
$p(\beta)=(U\otimes\nu_\beta)(G)$, then

```math
\boxed{
{d\over d\beta}\log p(\beta)
=\mathbb E_{U\otimes\nu_\beta}
[\langle A,D\rangle\mid G]
-\mathbb E_{\nu_\beta}\langle A,D\rangle.
}
\tag{R24.22}
```

At $\beta=0$ the covariance is nonnegative: conditional on a child cut, the
uniform outside energy has mean zero, and $G$ is an upper level set of the
child energy. There is no generic sign at positive temperature.

#### 7. Finite audits and a generic monotonicity counterexample

For $A_9$ at matched temperatures, the optimized soft selector cost in
(R24.9), the single joint path in (R24.13), and their common
$\beta\to\infty$ limit are:

| $\beta$ | optimized soft cost | joint-path cost |
|---:|---:|---:|
| 0.10 | 0.01435953 | 0.01437331 |
| 0.50 | 0.59565431 | 0.60143860 |
| 1.00 | 1.17865599 | 1.17983260 |
| 2.00 | 1.38318529 | 1.38318741 |
| 4.00 | 1.39074373 | 1.39074373 |

For the uniform parent-ground prior, the hard global probabilities and exact
collapsed costs are

| tolerance $a$ | $\overline\Omega$ | $-\log\overline\Omega$ |
|---:|---:|---:|
| 0 | 0.24888889 | 1.39074871 |
| 4 | 0.58666667 | 0.53329848 |
| 8 | 0.86222222 | 0.14824224 |
| 12 | 0.96888889 | 0.03160534 |

Thus the soft and hard selector collapses agree in the zero-tolerance,
zero-temperature limit, as their exact formulas predict.

Monotonicity in parent temperature is false for generic signings. On five
vertices, in lexicographic edge order

```text
(01,02,03,04,12,13,14,23,24,34),
```

take signs

```text
(1,1,1,1,1,1,1,1,-1,-1).
```

For $S=\{0,1,4\}$, tolerance $a=4$, and $\beta=1/2$,

```math
p(\beta)=0.6791552231529845\ldots,
\qquad
\operatorname{Cov}_{\nu_\beta}
(\langle A,D\rangle,\mathbf1_G)
=-0.0867533751995673\ldots.
```

Therefore (R24.22) is negative. This signing has $Q(A)=12$ while $q_5=8$,
so it is **not** an exact minimizer. Exhaustive exact minimizers through
order six and the audited $A_9$ tolerances show no such failure. Hence the
generic monotonicity claim is falsified, while an exact-minimizer-specific
version remains open. Monotonicity alone would still be insufficient: it
does not supply the quantitative increase from the typically tiny
$\beta=0$ near-ground fraction to (R24.5) or (R24.8).

#### Disposition

The main rigorous advance is that selector adaptation is no longer an
independent obstacle. The hard common-prior objective is exactly one global
endpoint tail, and the soft objective has both an optimized local full-path
form (R24.10) and a single normalized interpolation (R24.13).

For the original parent-Gibbs route, the missing theorem is either the hard
endpoint probability (R24.5) or the endpoint transport estimate (R24.19).
For the stronger soft-row downstream interface, it is the relaxed joint
ground-prior condition (R24.8). The latter permits overlap cost
$O(n^{3/4-c})$ but requires a conditioned row-square moment; overlap alone
does not provide it. Mean truncation is circular, generic temperature
monotonicity is false, and the existing quadratic/conditional-KL endpoint
bounds do not control the exponential density in (R24.18).
