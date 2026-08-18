# Independent audit of the diffuse common-phase obstruction

Status: **passed**, with the scope qualification in Section 6.  This note
checks
[`generic_diffuse_common_phase_row_lifetime_no_go.md`](generic_diffuse_common_phase_row_lifetime_no_go.md)
without using its derivation as a black box.  The construction is a generic
rank-one channel, not a sequence of quadratic Gibbs sectors of actual
minimizing children.

## 1. Factorwise conditional spread

Write `X=sigma xi`, where `sigma` is fair and the coordinates of `xi` are
independent with crossover probability `epsilon_L`.  Conditional on an
exterior value `X_(U^c)=x`, only the posterior weight of `sigma` changes.
Conditional on either value of `sigma`, the coordinates in `U` remain
independent and every atom has probability at most
`(1-epsilon_L)^|U|`.  A mixture of the two conditional laws has the same
upper bound.  Therefore

```math
\left\|\mathcal L(X_U\mid X_{U^c})\right\|_\infty
\le (1-\epsilon_L)^{|U|}.                         \tag{VDC.1}
```

The identical argument applies to `Y`.  This is stronger than marginal
spread: it holds after revealing every complementary coordinate.  Thus the
example genuinely survives the factorwise conclusion of Theorem 37.64;
there is no hidden frozen positive-density block.

## 2. Exact erased-row law

The global sign in `X_iY` is fair and independent of the right BSC word.
For a bridge row `z`, put `rho=tanh u`, `b=1-2epsilon_R`, and
`v=arctanh(b rho)`.  Direct product expansion gives

```math
\begin{aligned}
{E_\eta\cosh(u\langle z,\eta\rangle)\over(\cosh u)^n}
&={1\over2}\prod_j(1+b\rho z_j)
  +{1\over2}\prod_j(1-b\rho z_j)\\
&=(1-b^2\rho^2)^{n/2}\cosh\!\left(v\sum_jz_j\right).
                                                               \tag{VDC.2}
\end{aligned}
```

The prefactor cancels from the inverse escort.  Hence the canonical row
factor is exactly

```math
{d\nu_N\over dU_n}(z)
={\cosh(vS(z))^{-\lambda}\over
  E_{U_n}\cosh(vS)^{-\lambda}}.                   \tag{VDC.3}
```

Since `v sqrt(n) -> b beta sqrt(1-theta)`, the triangular-array CLT applies.
The numerator in the entropy formula involves
`log cosh(g) exp(-lambda log cosh(g))`, which is bounded for positive
`lambda`; the normalizer integrand is bounded by one.  Consequently the
claimed one-row KL limit `d_0(b)` follows without an unproved uniform
integrability step.  It is strictly positive whenever `b beta lambda` is
nonzero and is continuous in `b` near one.

## 3. Quantitative noisy-to-hard comparison

The all-positive latent BSC word has probability at least

```math
w_0=(1-\epsilon_L)^m(1-\epsilon_R)^n.
```

Its contribution to the mixture immediately gives
`log p-log p_0 >= log w_0=-A_NN` pointwise.

For the other direction, let `Delta_Q=Q-Q_0`.  A one-bit flip in (VDC.3)
changes the log density by at most `2lambda v`.  This pointwise ratio bound
survives summing over any unrevealed coordinates, so under arbitrary
sequential conditioning every row-bit mean has magnitude at most

```math
\tanh(\lambda v)\le\lambda u.                      \tag{VDC.4}
```

Rows are independent under the canonical product.  Iterating the elementary
one-bit MGF bound therefore yields

```math
E_r e^{u\langle B,\Delta\rangle}
\le \exp\left\{{u^2\over2}\|\Delta\|_2^2
               +\lambda u^2\|\Delta\|_1\right\}. \tag{VDC.5}
```

For `Delta_Q in {0,+-2}^(m n)` with Hamming support `d`, this exponent is
at most `2(1+lambda)beta^2 d/N`.  The extra factor two used to remove the
absolute value contributes only `log 2`.  Finally, if `F,G` are the two BSC
flip counts,

```math
d(Q,Q_0)=nF+mG-2FG\le nF+mG.                       \tag{VDC.6}
```

The two binomial MGFs are exactly the two logarithms in (DC.12).  This
proves (DC.13), uniformly for both `P=U` and `P=r`, and its error coefficient
tends to zero as the two fixed crossover probabilities tend to zero.

## 4. Linear coherent product gap

For any row product `P`, the inverse-escort identity is

```math
D(P\Vert q)=D(P\Vert U)+\lambda E_P\log p+\log Z_q. \tag{VDC.7}
```

Taking `P=r` and comparing the infimum over row products with the admissible
fair product cancels the common normalizer and gives

```math
J-I^\leftarrow
\ge m d_N+\lambda(E_r-E_U)\log p.                  \tag{VDC.8}
```

For the hard two-word likelihood, constants cancel, `log cosh` is
nonnegative, and

```math
E_U\log\cosh\!\left(u\sum_{ij}B_{ij}\right)
\le u\sqrt{mn}.
```

Combining this with the two one-sided bounds from Section 3 gives precisely

```math
J-I^\leftarrow
\ge m d_N-\lambda\{(A_N+R_N)N+\log2+u\sqrt{mn}\}. \tag{VDC.9}
```

First choose the crossover probabilities positive but sufficiently small.
Continuity gives `theta d_0(b)>lambda limsup(A_N+R_N)`; then let `N` grow.
The square-root term is `o(N)`, so (VDC.9) proves a fixed positive linear
gap.  This is a continuity proof at the physical scale
`u=beta/sqrt(N)`, not an interchange with an assumed thermodynamic limit.

## 5. The claimed distinction from subgroup examples is real

The BSC prior is not uniform on a switching-group orbit.  Hence there is no
symmetry reason for averaged posterior retuning to vanish.  The finite
two-projective-word calculation in (DC.26)--(DC.27) is enough to disprove
such invariance.  Accordingly, this example does not rule out a future
coarse posterior observable; it rules out conditional factor spread alone
as a closure input.

## 6. Verdict and exact scope

The construction rigorously proves

```math
\boxed{
\begin{gathered}
\text{uniform exponential conditional min-entropy on every factor subset}
\\[-2pt]
\centernot\Longrightarrow
J-I^\leftarrow=o(N)
\end{gathered}}
```

for generic diffuse rank-one channels.  The canonical row law, the linear
gap, and all constants in the noisy-to-hard comparison pass the audit.
The construction is **not** shown realizable by actual contracted-temperature
minimizing children.  Its proper consequence is therefore a sharper new
SML: any actual-child closure must use an optimizer-specific constraint
beyond even factorwise conditional min-entropy, or expose the diffuse phase
through a genuinely lower-information posterior statistic.

## 7. Polynomial orbit quotient

The refinement DC.29--DC.30 is also exact. Row and column permutations
preserve the BSC prior, fair bridge law, likelihood, inverse escort, and
Bayes kernel. Bridge inversion preserves the inverse escort because the
latent prior is antipodally symmetric. A rank-one word is described, up to
simultaneous factor complementation, by the two factor Hamming weights;
row/column permutations are transitive at fixed weights. Passing to

```math
K=\min\{|\xi^-|,m-|\xi^-|\},
\qquad L=\min\{|\eta^-|,n-|\eta^-|\}               \tag{VDC.10}
```

merges the two antipodal orbits, and bridge inversion gives them equal
conditional weight. Hence both `mu` and `bar mu` are uniform conditional
on `(K,L)`. Applying the KL chain rule proves

```math
D(\bar\mu\Vert\mu)
=D((K,L)_\#\bar\mu\Vert(K,L)_\#\mu).               \tag{VDC.11}
```

The quotient has at most
`(floor(m/2)+1)(floor(n/2)+1)=O(mn)` states. This refinement does not
weaken the no-go: it shows that the diffuse coherent phase is invisible to
atom spread but remains exposed by a genuinely low-information statistic.
