# Independent audit of the flat-threshold entropy theorem

Date: 2026-09-06, approximately 19:18 UTC. Status: **PASS** for
`transfer_seed_threshold_entropy_resampling_2026_09_06.md`, including
its uniform finite theorem, rational factor 22, and pressure consequence.
The operator mechanism was reconstructed independently before the source
artifact arrived. The source's own constants and count bound were then
checked in full; no correction was required.

## 1. Gaussian densities, smoothing direction, and the L2 domain

Write `B=A/sqrt(N)`, `tau=s/r`, and let `L_tau` be the density of
`N(0,I+tau B)` relative to standard Gaussian measure. The hypothesis
`|s|C/r<1` places every covariance eigenvalue strictly between zero and
two. This is the precise sufficient condition for `L_tau` to belong
to Gaussian `L2`, not just for its covariance to be positive definite.

Under the Gaussian Markov kernel with coefficient `alpha`, a Gaussian
measure of covariance `R` becomes one of covariance
`alpha^2 R+(1-alpha^2)I`. The same kernel acts on densities because it
is self-adjoint relative to its invariant standard Gaussian measure.
With `alpha=sqrt(r)`, this proves the required density identity
`L_s=T_sqrt(r)^(tensor N)L_tau`. There is no reversal of the smoothing
parameter and no appeal to a pointwise comparison of likelihoods.

## 2. Threshold projection and tensor PSD comparison

Under one standard Gaussian coordinate, let
`u=(sign(G-t)-m)/sqrt(v)`. The threshold conditional-expectation
projection is exactly `P=|1><1|+|u><u|`. The discrete likelihood ratio
is the Gaussian likelihood conditioned on the threshold cells. Hence

```math
1+\chi^2(p_s\Vert p_0)
=\langle L_\tau,(T_{\sqrt r}PT_{\sqrt r})^{\otimes N}L_\tau\rangle.
```

Put `w=T_sqrt(r)u`. It is centered and has squared norm `kappa`.
On the centered subspace, Cauchy--Schwarz gives
`|w><w| <= kappa I`. Therefore

```math
0\le T_{\sqrt r}PT_{\sqrt r}
\le \Pi+\kappa(I-\Pi),
```

where `Pi` is projection onto constants. This is a Hilbert-space
positive-semidefinite order, not pointwise domination and not domination
of two-coordinate probability densities. Tensoring is valid: the
difference of the N-fold tensor products is a sum of positive tensor
products, each containing one positive difference factor. All operators
are bounded and `L_tau` is in `L2`, so the quadratic-form use is legitimate.

The scalar is the normalized covariance of two threshold bits at
Gaussian correlation `r`. Differentiating the bivariate Gaussian density
in its correlation equals differentiating it once in each coordinate;
integrating over the two threshold tails leaves the boundary density.
Thus its derivative is
`phi(t)^2 exp(t^2 r/(1+r))/sqrt(1-r^2)`. Accounting for the factor four
in the covariance of sign bits and dividing by `v` gives precisely

```math
\kappa=\eta\int_0^r
       \frac{\exp(t^2z/(1+z))}{\sqrt{1-z^2}}\,dz.
```

Equivalently its full Hermite expansion is a sum of nonnegative
coefficients times `r^k`, with coefficient sum one and `k>=1`.
In particular `0<kappa<1`. No higher Hermite terms are discarded.

## 3. Resampling, determinant normalization, and flatness

The upper comparison operator is `kappa I+(1-kappa)Pi` in each
coordinate. Its product is an average of orthogonal conditional-
expectation projections onto independently retained coordinate sets.
Consequently the quadratic form equals the average squared norm of
the marginal likelihood on `S~Ber(kappa)^(tensor N)`.

Integrating the ambient likelihood over the omitted coordinates yields
the Gaussian marginal with covariance `I+tau B_S`; it does not yield a
Gaussian conditional covariance or a Schur complement. Direct integration
gives

```math
\|L_{I+H}\|_2^2
=\det(I+H)^{-1}\det(2(I+H)^{-1}-I)^{-1/2}
=\det(I-H^2)^{-1/2}.
```

This confirms both the square-root power and the admissible domain
`||H||op<1`. Compression does not increase the operator norm. Summing
`-log(1-x)<=x/(1-(|s|C/r)^2)` over the squared eigenvalues therefore
gives the stated coefficient
`lambda=s^2/[2(r^2-s^2 C^2)]`.

The decisive signing-specific identity is
`tr B_S^2=|S|(|S|-1)/N`. It uses the hollow diagonal and entrywise
magnitude `1/sqrt(N)`. It is not implied by an operator bound alone.
The resulting finite determinant/count inequality in the source is exact
up to the explicitly stated logarithmic eigenvalue bound.

## 4. Bernoulli count bound

The elementary binomial type bound, followed by summing at most `N+1`
types, yields the variational expression
`sup_x {lambda x^2-D(Ber(x)||Ber(kappa))}`. Endpoint derivatives force
every maximizer into `(0,1)`. Uniqueness or concavity is not needed.
Its stationary equation implies successively

```math
x\le b=\frac{\kappa e^{2\lambda}}{1-\kappa},\qquad
x\le\frac{\kappa e^{2\lambda b}}{1-\kappa}.
```

Discarding the nonnegative divergence now gives the source's exact
`lambda kappa^2 F_B` upper bound. All quantities are finite for fixed
parameters. The entropy inequality
`D(p_s||p_0)<=log(1+chi^2(p_s||p_0))` is Jensen under `p_s`, with
the correct measure and the required additive one.

As a separate algebra sanity check, solving the stationary equation
for `x` with its full denominator gives the slightly stronger bound
`x<=kappa exp(2lambda kappa exp(2lambda))`. This is not needed for the
audited factor 22; the source's more conservative `F_B` already suffices.

## 5. Rational constants at the fixed parameters

For `delta=2^-24`, `C=33/32`, `s=1/32`, `r=1/25`, all displayed
inequalities in the source checker were independently followed.

- For the lower threshold endpoint, use the lower Mills bound,
  `sqrt(2pi)<3`, and `e^(25/2)<(11/4)^13`. This yields
  `Phi(-5)>5/[78(11/4)^13]>delta`.
- For the upper endpoint, use the upper Mills bound,
  `sqrt(58pi)>13`, `e>8/3`, and `sqrt(e)>8/5`. This yields
  `Phi(-sqrt(29))<5/[104(8/3)^14]<delta`.
- Therefore `5<z=-t<sqrt(29)`. The function `z+1/z` increases for
  `z>1`; both Mills inequalities give
  `25/4 < x=a^2/v^2 < 225/[29(1-delta)^2] < 8`, hence `eta<32delta`.
- The positive Taylor series through degree 12 with its geometric
  remainder proves `exp(29/25)<16/5`. The squared rational inequality
  `(55/29)^2*(625/624)<(19/10)^2` then proves
  `kappa/(eta r)<19/10`.
- `lambda=320000/367951<1` and `kappa<10^-6`. Since
  `e^(2lambda)<e^2<(11/4)^2<8`, the source's bounds `b<9*10^-6`
  and `F_B<1001/1000` follow using `exp(y)<=1/(1-y)`.

The resulting entropy coefficient is exactly bounded by

```math
\frac{23682154496}{1149846875}<22.
```

No decimal tuning or numerical integration is needed for this result.

## 6. Deterministic pressure consequence

Switching an absolute maximizing spin to all ones, with its energy
sign, preserves hollow signing entries and operator norm. The threshold
covariance has diagonal one, so every marginal still has mean `m`.
Consequently the entropy identity is exactly
`H(X)=N h(delta)-D(p_s||p_0)`.

For a centered threshold function, Parseval bounds the absolute sum of
its squared Hermite coefficients by `v`. Thus, for either sign of a
pair correlation, the terms after the linear one have absolute value
at most `v rho_ij^2`. Summing over the `N(N-1)/2` undirected edges and
multiplying by physical inverse temperature `beta/sqrt(N)` gives at
most `beta v s^2 sqrt(N)/2`, as claimed.

The linear correlation contribution is
`beta a^2 s (N-1)/2`, not `beta a^2 s N/2` exactly. The source retains
the missing constant as `-beta a^2 s/2`. The remaining entropy error
is `-log(N+1)`. These are the three error terms in its pressure formula,
and all are subextensive at fixed parameters. The selected phase is
bounded by the sum of both phase partition functions, so no factor
two or trace correction is missing.

For `8<=beta<=33/4`, the source's quadratic lower bound increases on
`x in [25/4,8]`. It gives
`g_new/v^2>=18725/32768`, while the prior clipped credit is at most
`1089/2048` times `v^2`. Hence the strict gain is exactly at least

```math
g_{new}-g_{clip}\ge\frac{1301}{32768}v^2>0.
```

The retained-weave limits `beta=8/sqrt(31/32)` and normalized operator
norm `1/sqrt(31/32)` lie strictly inside the stated rational ranges.
Thus sufficiently large orders satisfy the hypotheses. Integration
into the all-order upper proof is owned by the root; this audit does
not assert a new decimal endpoint or settle convergence.

## Replay and preservation

Ran `computations/transfer_seed_threshold_entropy_resampling_verify_2026_09_06.py`:
PASS for all rational comparisons, 10 finite exact-formula threshold
instances checked numerically, and 56 general flat determinant instances.
The finite floating diagnostics are regression checks, not the proof of
the dimension-uniform theorem. Additional rational cross-checks were
evaluated without creating any files. This audit is the only new file
from this task; no temporary research output was created.
