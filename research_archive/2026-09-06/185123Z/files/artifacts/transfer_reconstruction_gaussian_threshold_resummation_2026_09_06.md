# A valid attenuated spectral resummation through Gaussian thresholds

Date: 2026-09-06. The unattenuated fixed-variance Gaussian comparison fails
on a finite actual signing. A threshold-channel construction nevertheless
proves a dimension-uniform attenuated version with `O(sqrt N)` error.
At the current very small flip density this theorem is rigorously weaker
than the clipped correlated-cluster credit; no further cap improvement is
claimed from it.

## 1. A finite failure of the unattenuated covariance functional

Write

```math
\mathcal S(J)=\sup_{R\succ0,\ \operatorname{diag}R=1}
 \frac12\{\operatorname{tr}(JR)+\log\det R\}.
```

Take the negative triangle, `A_ij=-1` for `i!=j`, and inverse temperature
`b=log2` (thus `beta=sqrt3 log2` in the normalization below). Its globally
optimal product means are zero. Indeed `lambda_max(A)=1` and
`h((1-m)/2)<=log2-m^2/2`, so every nonzero vector has strictly smaller
product value because `b<1`.

The exact partition is `Z=2e^{-3b}+6e^b=49/4`, hence its gap over the
product value is `log(49/32)`. The positive definite correlation matrix
with off-diagonal entries `-1/3` has determinant `16/27`, and gives

```math
\mathcal S(bA)\ge\log2+\frac12\log(16/27).
```

Its excess over the actual gap is at least

```math
\frac12\log\left(\frac{65536}{64827}\right)>0.           \tag{1}
```

This disproves a universal finite lower comparison even at a global
product optimum. Here `||bA||op=2log2>1`, so it does not refute a strict
weighted-norm hypothesis or a genuinely small-variance asymptotic theorem.
The reconnaissance script tested 240 cases of orders three through eight;
its other optimizers are only numerical multistart estimates. No failure
was found among those sampled cases with weighted norm below one. That
observation is not a theorem, and small-order tiny-variance tests are
particularly weak because the global spin-reversal contribution is `O(1)`.

## 2. Finite threshold-channel theorem

Let `A` be a symmetric hollow sign matrix of order `N`, and let
`b=beta/sqrt N`, `beta>0`. For arbitrary product means `m_i in (-1,1)`,
put

```math
P_A(m)=\sum_i h\!\left(\frac{1-m_i}{2}\right)+\frac b2m^TAm,
\quad t_i=\Phi^{-1}\!\left(\frac{1-m_i}{2}\right),
\quad a_i=2\phi(t_i),\quad v_i=1-m_i^2,\quad w_i=v_i-a_i^2,
```

where `Phi,phi` are the standard Gaussian distribution and density.
Define

```math
J_{\rm eff}=\frac\beta{\sqrt N}\operatorname{diag}(a)A\operatorname{diag}(a),
\qquad \kappa=\|J_{\rm eff}\|_{op}<1,
\qquad w_* =\max_i w_i.
```

Then

```math
\log Z_A(b)-P_A(m)
\ge\mathcal S(J_{\rm eff})
 -\frac{\beta w_*\kappa}{2(1-\kappa)}\sqrt N.             \tag{2}
```

No stationarity assumption is needed. In particular, (2) may be applied
at a globally optimal product mean or at the prescribed uniform cluster
around a maximizing spin. If `beta` is bounded and `kappa` stays below
one by a fixed margin, its error is `o(N)` uniformly over all dimensions,
signings, means, and thresholds. No uniform Gaussian approximation of
the spins is asserted.

### Proof: an actual spin law and its entropy

For a positive definite correlation matrix `R`, sample `G~N(0,R)` and
put `X_i=sign(G_i-t_i)`. Each spin has mean `m_i`. Under the same
coordinatewise map, `N(0,I)` becomes the independent product law of these
means. Relative-entropy data processing gives

```math
D(\mathcal L(X)\Vert\mu_m)
\le D(N(0,R)\Vert N(0,I))=-\frac12\log\det R.            \tag{3}
```

Data processing here is just conditional Jensen for the Gaussian density
ratio on each threshold cell. Since the means agree,
`H(X)=sum_i h((1-m_i)/2)-D(L(X)||mu_m)`.

### Proof: Hermite covariance and its uniform remainder

Expand the centered threshold function in orthonormal Gaussian Hermites.
Its first coefficient is `a_i=E[G sign(G-t_i)]=2phi(t_i)`, by one
integration by parts, and the sum of all squared coefficients is `v_i`.
For two standard Gaussians with correlation `r`, the Hermite identity
`E[H_k(G_1)H_l(G_2)]=1_{k=l} r^k` follows by comparing coefficients in
`E exp(sG_1-s^2/2) exp(tG_2-t^2/2)=exp(rst)`.
Thus Parseval and Cauchy--Schwarz give

```math
\operatorname{Cov}(X_i,X_j)=a_i a_jR_{ij}+E_{ij},\qquad
|E_{ij}|\le\sqrt{w_iw_j}\,R_{ij}^2.                     \tag{4}
```

In particular `w_i>=0`. The estimate is uniform even as a mean approaches
the boundary of the spin interval.

The Gibbs variational inequality, (3), (4), and `|A_ij|=1` off the
diagonal therefore yield, for every such `R`,

```math
\log Z_A(b)-P_A(m)
\ge\frac12\{\operatorname{tr}(J_{\rm eff}R)+\log\det R\}
 -\frac\beta{2\sqrt N}\sum_{i\ne j}\sqrt{w_iw_j}R_{ij}^2. \tag{5}
```

The means being fixed is what removes all local-field contributions.

### Proof: controlling the optimizing Gaussian correlation matrix

The maximum defining `S(J_eff)` exists in the interior of the correlation
cone: its closure is compact, while `log det R` tends to negative infinity
at a singular boundary. At its maximizer,

```math
R^{-1}=\Lambda-J_{\rm eff},\qquad \Lambda\text{ diagonal}.
```

Since `J_eff` is hollow and `R_ii=1`, the Schur complement for the `i`th
diagonal entry gives `Lambda_ii=1+u_i^T K_{-i}^{-1}u_i>=1`, where
`K=R^{-1}`. Hence

```math
R\preceq(1-\kappa)^{-1}I,\qquad
\sum_{i\ne j}R_{ij}^2=\operatorname{tr}R^2-N
\le\frac{N\kappa}{1-\kappa}.                            \tag{6}
```

Substitute (6) into (5). This proves the precise error in (2); the
argument uses `tr R=N`, not a dimension-dependent entrywise estimate.

## 3. Explicit pressure-to-cap mapping and attenuation

For uniform means of magnitude `q=1-2delta`, let
`a=2phi(Phi^{-1}(delta))` and `v=4delta(1-delta)`.
For any hollow `J` with norm `kappa<1`,

```math
\frac{\operatorname{tr}J^2}{4(1+\kappa)}
\le\mathcal S(J)\le
\frac{\operatorname{tr}J^2}{4(1-\kappa)}.                 \tag{7}
```

For the lower bound, use `R=I+J/(1+kappa)` and the scalar inequality
`log(1+x)>=x-x^2/[2(1-|x|)]`. For the upper bound, the Gaussian diagonal
dual permits `Lambda=I`, and the logarithmic series bounds
`-logdet(I-J)/2` using `tr J=0`. These arguments prove (7) without
approximating the optimizer.

When `||A||op<=C sqrt N`, one has `kappa<=beta C a^2`. Thus (2), (7)
give the following asymptotic credit per spin, if `beta C a^2<1`:

```math
g_{\rm cop}=\frac{\beta^2a^4}{4(1+\beta C a^2)}.          \tag{8}
```

For the retained weave, `beta=2t/sqrt p+o(1)` and
`C<=1/sqrt p+o(1)`. An annealed pressure `A_*m^2+o(m^2)` therefore
implies the valid cap extraction

```math
c\le\frac{A_*-p h(\delta)-p g_{\rm cop}}
              {2t\sqrt p(1-2\delta)^2}.                 \tag{9}
```

All parameters and the mean vector choice precede the order limit. The
`O(sqrt N)` error in (2) disappears on the `m^2` pressure scale. As usual,
choose the appropriate objective sign and gauge around an absolute
maximizer; the argument applies to both signs separately.

At the current `p=31/32,t=4,delta=2^{-24}`, the threshold is approximately
`-5.29470408485`, `a^2≈4.25515912054e-13`, and
`a^2/v≈1.78474319837e-6`. Consequently (8) is approximately
`2.990472942e-24`, versus the already proved clipped-law credit
`beta^2v^2/128≈2.933853527e-14`. These decimal evaluations are diagnostics.

The inferiority can also be certified without a Gaussian quantile oracle.
Write `z=-Phi^{-1}(delta)`. Elementary tail bounds give `1<z<6` at this
delta. Integration by parts gives the lower Mills bound
`Phi(-z)>=phi(z)z/(1+z^2)`. Therefore

```math
\frac{a^2}{v}\le\frac{1369\delta}{36(1-\delta)}<51\delta.
```

For completeness, `Phi(-6)<e^{-18}/6<delta` follows from `e>5/2`;
`Phi(-1)>=int_1^2 phi(u)du>=phi(2)>1/27>delta` follows from
`e<3` and `sqrt(2pi)<3`. These are very loose bounds but sufficient.
Here `beta C a^2<=beta C v<1/2`. Even the upper bound in (7) is therefore
at most the following fraction of the clipped-law credit:

```math
64(51\delta)^2<10^{-8}.                                 \tag{10}
```

Thus the fully optimized **attenuated** Gaussian functional certified by
(2), not merely its quadratic approximation, cannot improve the current
clipped-law coefficient. This does not bound the actual entropy of the
thresholded law: the Gaussian data-processing entropy estimate may lose
information. A useful stronger theorem would have to recover that lost
binary entropy, or use a different correlation channel.

## 4. Primary-source scope check

The diagonal covariance functional is closely related to adaptive TAP and
expectation-consistent approximations. Opper and Winther explicitly separate
these approximations from variational bounds and explain that the approximate
free energy need not bound the true one; their discussion in Sections 3--4
does not supply the desired universal Ising lower inequality. This agrees
with the finite counterexample (1).
[Opper--Winther, JMLR 6 (2005), Sections 3--4](https://www.jmlr.org/papers/volume6/opper05a/opper05a.pdf).

Gaussian-copula variational proposals are established methodology. The
Han--Liao--Dunson--Carin treatment uses continuous marginal transformations;
it is not a theorem about the threshold channel's discrete entropy or the
uniform flat-interaction error in (2). Those steps were proved directly here.
[Variational Gaussian Copula Inference, AISTATS 2016](https://proceedings.mlr.press/v51/han16.pdf).

## 5. Reproduction and preservation

`computations/transfer_reconstruction_gaussian_threshold_exact_checks_2026_09_06.py`
certifies the negative-triangle witness and the current-parameter attenuation
comparison in rational arithmetic, and separately labels finite Gaussian
integral regressions as numerical.

The reconnaissance script is
`computations/transfer_reconstruction_fixed_variance_gaussian_probe_2026_09_06.py`;
its unverified multistart output and seed are preserved in
`computations/results/transfer_reconstruction_fixed_variance_gaussian_probe_2026_09_06.json`.
