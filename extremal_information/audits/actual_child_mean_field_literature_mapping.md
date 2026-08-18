# Primary-source audit: mean-field decompositions versus the actual bridge escort

Status: **retrieval-grounded mapping audit**.  The primary arXiv sources and
their TeX were read, including theorem statements and proofs where needed.
The conclusion is negative but sharp: the actual directed product projection
is exactly a mean-field variational gap, yet none of the three requested
theorems gives an `o(N)` or even a competitive `O(N)` bound for the bridge
Hamiltonian.  Mixture conclusions are not silently converted into a
single-product reverse-KL conclusion.

Primary sources:

- Eldan--Gross, [*Decomposition of mean-field Gibbs distributions into
  product measures*](https://arxiv.org/abs/1708.05859), especially the Main
  Structural Theorem and its definitions of discrete gradient complexity.
- Eldan, [*Taming correlations through entropy-efficient measure
  decompositions with applications to mean-field approximation*](https://arxiv.org/abs/1811.11530),
  especially the Main Decomposition Theorem, its mean-field theorem, and its
  Schatten/rank corollaries.
- Jain--Koehler--Risteski, [*Mean-field approximation, convex hierarchies,
  and the optimality of correlation rounding*](https://arxiv.org/abs/1808.07226),
  especially Theorem 1 and the general `k`-MRF theorem.

Because the Eldan--Gross paper explicitly imports its decomposition and
partition estimates from Eldan's earlier paper, the corresponding primary
source was also checked:

- Eldan, [*Gaussian-width gradient complexity, reverse log-Sobolev
  inequalities and nonlinear large deviations*](https://arxiv.org/abs/1612.04346),
  especially its mean-field Corollary 1.

## 1. Exact bridge mapping

Fix contracted-temperature children `A,D`, an orientation `epsilon`, and a
comparable split `m+n=N`.  Put

```math
d=mn,
\qquad t={\beta\over\sqrt N},
\qquad f(B)=-\lambda L_\epsilon(B),
\qquad
L_\epsilon(B)=\log\overline Z_N(A,\epsilon D,B;t).       \tag{MF.1}
```

The bridge escort is the Gibbs law

```math
q(B)={e^{f(B)}\over\mathbb E_{U_d}e^f},U_d(B).          \tag{MF.2}
```

Let

```math
\mathcal F(f)=\log\mathbb E_{U_d}e^f
```

and let `F_bit(f)` be the same Gibbs variational problem restricted to
products of the `d` individual bridge bits.  The Gibbs identity gives

```math
\boxed{
\mathcal F(f)-\mathcal F_{\rm bit}(f)
=\inf_{p\ {\rm bit\ product}}D(p\Vert q).}              \tag{MF.3}
```

Thus the usual naive mean-field deficit is **exactly** a reverse KL
projection.  If `P_row` denotes products of the `m` whole bridge rows, then

```math
\mathcal I_\lambda^{\leftarrow}
:=\inf_{p\in\mathcal P_{\rm row}}D(p\Vert q)
\le \mathcal F(f)-\mathcal F_{\rm bit}(f),              \tag{MF.4}
```

because every bit product is a row product.  A published bitwise mean-field
bound would therefore upper-bound the row projection, but not conversely.

The task-local exact theorem already proves

```math
0\le\mathcal I_\lambda^{\leftarrow}
\le D(U_d\Vert q)
\le {\lambda^2t^2d\over2}=O_{\beta,\lambda}(N).         \tag{MF.5}
```

Any imported result must beat (MF.5), ideally by proving `o(N)`.

## 2. Gradient geometry, with the discrete-gradient correction

Extend the bridge entries to `b in [-1,1]^d`.  Under the parent Gibbs law at
`b`, let

```math
Q=\tau xy^T\in\{\pm xy^T:x\in\{\pm1\}^m,
                           y\in\{\pm1\}^n\}.
```

Direct differentiation gives

```math
\nabla f(b)=-\lambda t\,\mathbb E_bQ
\in\lambda t\operatorname{conv}\{\pm xy^T\},           \tag{MF.6}
```

and

```math
\nabla^2f(b)=-\lambda t^2\operatorname{Cov}_b(\operatorname{vec}Q).
                                                                  \tag{MF.7}
```

In particular,

```math
\operatorname{GW}(\{\nabla f(b)})
\le\lambda t\,\mathbb E\max_{x,y}|x^TGy|
\le\lambda t\sqrt{mn}(\sqrt m+\sqrt n),               \tag{MF.8}
```

where `G` is an `m by n` standard Gaussian matrix and
`E||G||op<=sqrt(m)+sqrt(n)` was used.

The Eldan papers use the **discrete** half-difference gradient, not the
continuous gradient in (MF.6).  This distinction can be repaired, but must
not be omitted.  Since

```math
|\partial_{ee}^2f|\le\lambda t^2,
```

the discrete secant in coordinate `e` differs from the continuous derivative
at the same vertex by at most `lambda t^2`.  Therefore

```math
\boxed{
\operatorname{Comp}(f)
\le\lambda t\sqrt{mn}(\sqrt m+\sqrt n)
 +\sqrt{2/\pi}\,\lambda t^2mn
=O_{\beta,\lambda}(N).}                               \tag{MF.9}
```

Also

```math
\operatorname{Lip}(f)\le\lambda t=O(N^{-1/2}).          \tag{MF.10}
```

The discrete-gradient `l_1` Lipschitz constant used for the Eldan--Gross
fixed-point conclusion is only bounded here by

```math
L_2\le\max\{1,C\lambda t^2d\}=O_{\beta,\lambda}(N).     \tag{MF.11}
```

The low Gaussian width is real; it is not, by itself, at the scale needed
for the pressure problem.

## 3. Eldan--Gross 1708.05859

### Precise imported theorem

For a Boolean Hamiltonian on an ambient cube of dimension `d`, their Main
Structural Theorem sets

```math
D=\operatorname{Comp}(f),
\quad L_1=\max(1,\operatorname{Lip}f),
\quad
L_2=\max\left(1,
 {\|\nabla f(x)-\nabla f(y)\|_1\over\|x-y\|_1}
 \right).
```

It gives an exact mixture into small linear tilts of the Gibbs law; all but
`3(D/d)^(1/4)` of the mixture is close in normalized Hamming transportation
to product measures, and its product means lie within

```math
5000L_1L_2^{3/4}D^{1/4}d^{3/4}                         \tag{MF.12}
```

in `l_1` of the mean-field fixed-point equation.

### Actual scaling

Equations (MF.9)--(MF.11) give

```math
(D/d)^{1/4}=O(N^{-1/4}),                               \tag{MF.13}
```

so the theorem does recognize a low-complexity mixture.  But its Hamming
transport error is `O(N^(7/4))`.  Since changing one bridge bit changes
`f` by at most `2lambda t=O(N^(-1/2))`, converting this transport error into
an energy error costs

```math
O(N^{-1/2}N^{7/4})=O(N^{5/4}),                         \tag{MF.14}
```

already larger than the target `N` scale.  The fixed-point radius (MF.12)
is worse: with `L_2=O(N)`, it is `O(N^(5/2))`, exceeding the ambient `l_1`
diameter and hence vacuous.

Most importantly, a mixture of nearby products is not a single product with
small `D(p||q)`.  The paper's statement alone does not bound (MF.3): a
continuous mixture can have no component of macroscopic latent weight, and
Wasserstein proximity does not imply reverse KL.  One needs the additional
entropy and same-marginal productization argument from the imported 2016
paper.  Applying that argument yields the next bound, not `o(N)`.

## 4. The exact Gaussian-width mean-field corollary

Eldan's 2016 Corollary 1 states, with the paper's exact constant, that some
bit-product law `p` satisfies

```math
\mathcal F(f)
\le \mathbb E_pf-D(p\Vert U_d)
 +64\operatorname{Lip}(f)^{2/3}
     \operatorname{Comp}(f)^{1/3}d^{2/3}.              \tag{MF.15}
```

By (MF.3), (MF.9), and (MF.10), this rigorously gives

```math
\mathcal I_\lambda^{\leftarrow}
\le O_{\beta,\lambda}(N^{4/3}).                        \tag{MF.16}
```

This is a legal single-product reverse-KL conclusion, but it is strictly
worse than the elementary `O(N)` estimate (MF.5).  Even a hypothetical
dimension-free `Comp(f)=O(1)` would make (MF.15) only `O(N)`; this corollary
would need `Comp(f)=o(1)` to give `o(N)` at the present Lipschitz and ambient
dimension scales.

**Verdict on Eldan--Gross:** genuine structural mixture information, but no
competitive control of the reverse row-product projection.  Quantitatively:
`O(N^(4/3))` through the exact imported mean-field corollary, versus the
already proved `O(N)`.

## 5. Eldan 1811.11530

### Precise imported theorems

The Main Decomposition Theorem applies to an arbitrary measure `mu` and any
positive definite matrix `L`.  It gives an exact mixture
`mu=int mu_theta dm(theta)` satisfying

```math
H(mu)-\mathbb E H(mu_\theta)
\le\log\det(\operatorname{Cov}(mu)L+I),                 \tag{MF.17}
```

```math
\mathbb E\operatorname{Cov}(mu_\theta)\preceq L^{-1},
\qquad
\mathbb E[\operatorname{Cov}(mu_\theta)L
          \operatorname{Cov}(mu_\theta)]
\preceq\operatorname{Cov}(mu).                         \tag{MF.18}
```

For a **quadratic** Ising/Potts Hamiltonian with fixed interaction matrix
`J`, its mean-field theorem then produces a product measure with free-energy deficit at
most

```math
3\log\det(\operatorname{Cov}(mu)|J|+I).                 \tag{MF.19}
```

Its corollaries include

```math
10{p+1\over p}(D_0^2d\|J\|_{S_p})^{p/(p+1)}            \tag{MF.20}
```

for single-spin support diameter `D_0`, and

```math
3\operatorname{rank}(J)
 \log(D_0^2d\|J\|_{op}+1).                             \tag{MF.21}
```

### Why it does not apply to the bridge escort

The bridge Hamiltonian `f=-lambda L` is not quadratic in `B`.  Formula
(MF.7) is a *state-dependent* Hessian, not a single interaction matrix.
The proof of (MF.19) uses quadraticity at exactly the decisive step:
replacing a component by the product of its marginals changes the energy by
`Tr(J Cov)`.  There is no such identity for `-lambda log Z`.

The arbitrary-measure decomposition (MF.17)--(MF.18) leaves its components
correlated.  It does not assert that they are products and therefore does not
control `inf_product D(p||q)`.  Calling it a product mixture would be an
incorrect strengthening.

For scale only, (MF.7) gives

```math
\|\nabla^2f(b)\|_{S_1}
\le\lambda t^2d=O(N).                                  \tag{MF.22}
```

Even if one illegally froze one Hessian and treated it as a quadratic
interaction, the `p=1` Schatten estimate (MF.20) would be `O(N^(3/2))` in
ambient dimension `d=Theta(N^2)`.  This is not an application, merely a
warning that the published scale is not close to `o(N)`.

**Verdict on Eldan 2018:** no bound on the actual `I^leftarrow`.  The general
theorem supplies an entropy/covariance mixture, while its product conclusion
requires a fixed quadratic interaction absent here.

## 6. Jain--Koehler--Risteski 1808.07226

### Precise imported theorems

For a quadratic Ising Hamiltonian on `d` Boolean sites with interaction
matrix `J`, Theorem 1 proves

```math
0\le\mathcal F-\mathcal F_{\rm bit}
\le3d^{2/3}\|J\|_F^{2/3}.                              \tag{MF.23}
```

The proof conditions at most `ell` coordinates, applies correlation rounding,
and uses the exact quadratic covariance error paired with `J`; optimizing
`2 epsilon d ||J||F+epsilon^(-2)` gives (MF.23).

For a fixed-order `k`-MRF over alphabet size `q`, their general theorem is

```math
\mathcal F-\mathcal F_{\rm bit}
\le3\left({k\log q\over\sqrt{k!}}
 d^{k/2}\|J\|_F\right)^{2/3}.                          \tag{MF.24}
```

Here `||J||F^2` is the sum of squared sup norms of the declared `k`-body
potentials.

### Mapping failure

Again, `f=-lambda L` is neither quadratic nor a fixed-order MRF.  Its Fourier
expansion contains interactions of unbounded order.  Correlation rounding of
pairs does not control expectation of this nonlinear function, and the proof
of (MF.23) cannot be repeated by substituting the gradient hull (MF.6).

As a scale diagnostic only, a frozen Hessian has Frobenius norm at most its
trace norm `O(N)` by (MF.22); putting that into (MF.23) would give `O(N^2)`,
and is not a valid application anyway.  Declaring the entire function to be
one interaction of order `d` makes (MF.24) vastly worse and carries no useful
asymptotics.

**Verdict on JKR:** no legal control of `I^leftarrow`; the crucial fixed-order
energy/covariance pairing is missing.

## 7. What the sources do and do not prove

| source | legal conclusion for actual escort | scale for `I^leftarrow` | obstruction |
|---|---|---:|---|
| Eldan--Gross 2017 | low-complexity tilt/product mixture in Hamming transport | no direct bound | mixture is not one reverse-KL product; target transport scale is too fine |
| Eldan 2016 mean-field corollary imported there | one bit-product reverse-KL competitor | `O(N^(4/3))` | generic Gaussian-width exponent loses `N^(1/3)` |
| Eldan 2018 arbitrary decomposition | entropy-efficient correlated components | no bound | components are not products |
| Eldan 2018 quadratic mean field | not applicable | n/a | `-lambda L` has a varying Hessian, not fixed `J` |
| JKR quadratic / fixed-order MRF | not applicable | n/a | unbounded-order nonlinear Hamiltonian |
| task-local bounded-odds theorem | uniform competitor and optimal row-product shadow | `O(N)` | does not distinguish `o(N)` from `Theta(N)` |

No requested theorem proves:

```math
\mathcal I_\lambda^{\leftarrow}=o(N),
\qquad
\mathcal I_\lambda^{\leftarrow}=\Omega(N),
```

or a corresponding `o(N)` row-product pressure gain.  The only rigorous
general scale remains

```math
0\le\mathcal I_\lambda^{\leftarrow}=O(N).              \tag{MF.25}
```

## 8. Narrowed literature-informed target

The missing theorem is not another generic product-mixture decomposition.
It must exploit the special *concave log-partition of rank-one features*

```math
f(B)=-\lambda\log\sum_z w_z
 e^{t\langle B,Q_z\rangle},
\qquad Q_z=\pm xy^T,                                   \tag{MF.26}
```

together, if needed, with the actual child-minimizer flip inequalities.
The precise open alternative is:

```math
\inf_{p\in\mathcal P_{\rm row}}D(p\Vert q)=o(N)
\quad\text{or}\quad
\inf_{p\in\mathcal P_{\rm row}}D(p\Vert q)=\Theta(N). \tag{MF.27}
```

The cited theories show why this is a bona fide nonlinear mean-field-gap
problem, but they do not decide it.  Any future use of them must add a theorem
that converts the moving covariance geometry (MF.7), not merely its Gaussian
width, into an `o(N)` row-product variational gap.
