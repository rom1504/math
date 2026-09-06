# Arbitrary-arity Schur factorization and the exact seed interface

Date: 2026-09-06. Independent seed-transfer track. The multivariate
inequality below is proved. It extends the scalar recursive certificate,
but does not identify the original minimax liminf with its limsup.

Write `Q(A)=max_x |x^T A x|/2` for hollow symmetric sign matrices, and
`a_n=M_n/n^(3/2)`.

## 1. A precise sufficient seed-transfer criterion

Suppose there are hollow sign seeds `A_j` of orders `k_j -> infinity`
such that

```math
\frac{Q(A_j)}{k_j^{3/2}}\longrightarrow\liminf_n a_n.
```

It is sufficient to have errors `epsilon_j -> 0` with the following
property: for each FIXED `j`, there is an unbounded sequence of actual
hollow sign matrices `C_(j,l)` of orders `N_(j,l)`, with

```math
\frac{N_{j,l+1}}{N_{j,l}}\longrightarrow1,
\qquad
\limsup_l\frac{Q(C_{j,l})}{N_{j,l}^{3/2}}
\le\frac{Q(A_j)}{k_j^{3/2}}+\epsilon_j.
```

Indeed principal restriction fills all sufficiently large orders with
an additional multiplicative factor tending to one, first at fixed `j`.
Thus `limsup_n a_n <= Q(A_j)/k_j^(3/2)+epsilon_j`; then let `j` tend
to infinity. No uniform construction-depth bound in `j` is necessary.
Only existence for one near-optimal seed sequence is needed. A theorem
for every near-optimal signing would be stronger.

The present recursive-weave theorem gives such a dense sequence with a
constant independent of `A_j`. It does not give the displayed
seed-dependent inequality. Selecting a better numerical ternary envelope
alone does not supply the missing dependence.

## 2. Scalar facts used in the multivariate extension

Let

```math
T_t(\nu)=\sup_L\{\mathbb E g_t(\operatorname{Var}(X\mid L))-I(X;L)\},
\qquad X\sim\nu.
```

This definition makes sense for any finite real law, not necessarily a
symmetric one. The Gaussian function `g=g_t` is decreasing and convex,
`g(0)=0`, and `h(s)=g(exp(s))` is concave. Explicitly,

```math
g'(v)=-t(1-\rho),\qquad
\frac{d}{ds}g(e^s)=-\frac{\rho}{2(1+\rho)},
\qquad 2tv=\frac{\rho}{1-\rho^2}.
```

In particular `-t <= g'(v) <= 0`, and the logarithmic derivative lies
between `-1/4` and zero. These identities include the continuous value
at `v=0`.

For a finite mixture `nu=sum_i p_i nu_i`, a label which first reveals the
mixture component and then uses a channel for that component proves

```math
T_t(\nu)\ge\sum_i p_i T_t(\nu_i)
                 -\sum_i p_i D(\nu_i\Vert\nu).                 (1)
```

Equivalently `T_t+H` is concave on each finite alphabet. This is a
one-sided inequality with the exact information payment, not a claim
that `T_t` itself is concave.

## 3. Cholesky rewards dominate every rotated diagonal

Let `Sigma` be a positive definite `d`-by-`d` covariance matrix. In the
chosen coordinate order let `delta_i` be its successive linear-regression
residual variances, so that `delta_i=det(Sigma_[i])/det(Sigma_[i-1])`.
Let `lambda_i` be its eigenvalues. For every real orthogonal `O`,

```math
\sum_i g(\delta_i)\ge\sum_i g(\lambda_i)
                    \ge\sum_i g((O\Sigma O^T)_{ii}).          (2)
```

For completeness, write `Sigma=R^T R` with upper-triangular Cholesky
factor `R`; its diagonal squares are `delta_i`. For any subset `J` of
`r` coordinates, the principal `J`-minor of `R` is the product of its
selected diagonal entries. Its absolute determinant is bounded by the
operator norm of the `r`th exterior power of `R`, namely the product of
its `r` largest singular values. Consequently the product of the `r`
largest `delta_i` is at most the product of the `r` largest `lambda_i`,
with equality for `r=d`. Thus `log(delta)` is majorized by `log(lambda)`.
Concavity of `h` gives the first inequality in (2). For the second,
diagonal entries of `O Sigma O^T` are convex combinations of the
eigenvalues with a doubly stochastic matrix of squared orthogonal
entries; ordinary convexity of `g` applies.

Positive semidefinite matrices follow by replacing `Sigma` with
`Sigma+epsilon I` and taking `epsilon` to zero. The actual conditional
variance after nonlinear regression is at most the linear-regression
variance. One can see the regularized comparison directly: for every
linear predictor coefficient vector `b`, adding `epsilon I` increases
its residual quadratic form by `epsilon(1+||b||^2)`. Hence no inverse
or singular-pivot convention is hidden in the limiting argument.

## 4. General orthogonal sequential-label inequality

Let `X=(X_1,...,X_d)` have any finite joint law and put `Y=OX`. Denote
total correlation by `TC(X)=sum_i H(X_i)-H(X)`. Then

```math
\boxed{\quad
\sum_i T_t(\mathcal L(Y_i))
\le\sum_i T_t(\mathcal L(X_i))+\operatorname{TC}(X).
\quad}                                                       (3)
```

Choose channels `L_i|Y_i`, independent conditionally on `Y`, and write
`L=(L_1,...,L_d)`. Label parent coordinate `X_i` by `(L,X_1,...,X_(i-1))`.
The exact information identity and bound are

```math
\sum_i I(X_i;L,X_{<i})=I(X;L)+\operatorname{TC}(X)
\le\sum_i I(Y_i;L_i)+\operatorname{TC}(X).                     (4)
```

The deficit in the second inequality is `TC(L)`, because the channels
are conditionally independent and `X <-> Y` is bijective.

At fixed `L`, Jensen, decreasingness of `g`, and conditional regression
bound the total parent reward below by the Cholesky sum in (2).
Equation (2) bounds this below by
`sum_i g(Var(Y_i|L))`. Refining `L_i` to `L` cannot decrease the
averaged reward, again by total variance, decreasingness, and Jensen.
Subtract (4), then take the supremum over the child channels. This
proves (3). The argument does not require Gaussian posteriors, equal
conditional variances, or equal child laws.

## 5. Every finite-arity orthogonal type operator has this supersolution

Fix a symmetric finite parent law `nu`, and a real orthogonal mixer `O`
of arity `d`. Let `pi` range over input laws whose average absolute
coordinate marginal equals `|nu|`. Let `mu_i` be the symmetrized law
of `(OX)_i`. Define

```math
(B_O f)(\nu)=\sup_\pi\left\{
\frac1d\sum_i f(\mu_i)-\frac1dD(\pi\Vert\nu^{\otimes d})\right\}.
```

Then

```math
\boxed{B_O T_t\le T_t.}                                    (5)
```

First average `pi` with its global reversal. Each symmetrized output law
is unchanged, while relative entropy decreases. Thus one may assume
`pi` is centrally symmetric. Its actual coordinate marginals `nu_i`
and its output marginals are then symmetric, and `d^(-1)sum_i nu_i=nu`.
Use (3), then (1), and the exact identity

```math
D(\pi\Vert\nu^{\otimes d})
=\operatorname{TC}(X)+\sum_iD(\nu_i\Vert\nu).
```

This yields (5). No permutation symmetry of the particular mixer `O`
is assumed. In particular, the admissible input laws need not have
equal coordinate marginals. These are the two potential gaps that a
naive extrapolation from the binary proof would leave open.

## 6. Where an actual seed can enter, and where it currently disappears

The exact recursive Hadamard realization replaces the binary mixer by
`O otimes I_s`, between independent flat orthogonal child bases. A
physical entry of the result has magnitude `|O_ij|/sqrt(s)`. Therefore
it is a full sign Hadamard after the intended normalization exactly
when every `|O_ij|=1/sqrt(d)`: the mixer must itself be a normalized
real Hadamard matrix. An arbitrary hollow minimizing signing does not
satisfy this condition, and adjoining a sign diagonal does not force
orthogonality.

The theorem is therefore stronger in arity but not yet seed-transfer:
its exact flat sign realization admits Hadamard mixers, while the
original seed datum is the Boolean support value `Q(A)`. Neither
orthogonality nor a suitable entrywise-flat orthogonal dilation follows
from this datum.

In the outer weave, independent edge-sign averaging also removes the
chosen outer sign matrix from the PSD-kernel bound. Even without that
averaging, a bound by the product of local Gaussian-kernel tensor norms
is insensitive to the edge reflection signs. A seed-sensitive theorem
must retain information beyond those local norms, or use a different
exact sign realization. The matching Boolean eigenvectors of the full
rank-one weave remain an additional obstruction before restriction.

## Status and next discriminating target

The exact useful extension is (3)--(5). The convergence criterion in
Section 1 remains open. The next issue is a quantitative nonorthogonal
or seed-conditioned factorization whose loss is `o(k^(3/2))` for a
near-optimal SEQUENCE of actual sign seeds, together with an exact
sign realization. Bounding such a loss only by operator norm, singular
values, or an absolute PSD majorant is not automatically a Boolean-cap
comparison.
