# A quantitative obstruction to finite additive mean representations

Date: 2026-09-07. Status: elementary theorem, pending independent audit.
This concerns an exact hypothesis of a classical mean-field limit theorem;
it does not disprove original convergence or a different interpolation.

Let A_n be actual hollow symmetric signings and `B_n=A_n/sqrt(n)`,
with `||B_n||op<=L` for a fixed L>0. Write
`E_n(x)=x^T B_n x/(2n)=H_(A_n)(x)/n^(3/2)`.

## 1. The theorem

For fixed k,M, let

```math
m_{n,l}(x)=b_{n,l}+\frac1n\sum_{i=1}^n a_{n,l,i}x_i,
\quad |a_{n,l,i}|\le M,\quad |b_{n,l}|\le M,
\quad 1\le l\le k.
```

For ANY fixed continuous `g:[-2M,2M]^k -> R`,

```math
\liminf_n\ \sup_{x\in\{\pm1\}^n}
 |E_n(x)-g(m_n(x))|\ \ge\ \frac1{4\pi L}.             (1)
```

Thus the deterministic critical-scale signing energy cannot be uniformly
approximated to o(1) by a fixed continuous function of finitely many
bounded additive spin statistics. Site-dependent coefficients, changing
with n, are explicitly allowed. A bounded single-spin function on
`{+1,-1}` is affine, so this includes arbitrary bounded additive
single-spin statistics, not just magnetization.

There is also a finite-dimensional quantitative form. If g is K-Lipschitz
in Euclidean distance, with k now permitted to depend on n, then

```math
\sup_x|E_n(x)-g(m_n(x))|
\ge \frac{n-1}{4\pi L n}
 -\frac{KM}{2}(1+\sqrt{3/2})\sqrt{\frac{k}{n}}.       (2)
```

In particular uniformly bounded K,M and k=o(n) still cannot achieve
vanishing uniform error. Fixed k continuous g is handled by uniform
continuity, without a Lipschitz assumption.

## 2. Two spin laws with indistinguishable additive statistics

Under uniform independent spins X, `E E_n(X)=0` and

`Var(m_(n,l)(X)) <= M^2/n`.

For a second law use a centered Gaussian vector with covariance

`R_n=I+B_n/(2L)`

and take its coordinate signs F. This is a genuine covariance matrix:
its spectrum lies in [1/2,3/2] and its diagonal is one. Gaussian
hyperplane rounding gives

`C_n=E F F^T=(2/pi) arcsin[R_n]`,

where the function is entrywise. Its operator norm is at most 3/2.
To check this last assertion explicitly, the odd Taylor coefficients
of arcsin are nonnegative and sum to pi/2. For a positive semidefinite
correlation matrix R, the Schur map `X -> R circ X` is positive and
unital, hence is a contraction on self-adjoint operator norm. Therefore
every positive Schur power of R has norm at most ||R||; summing the
arcsin series proves `||C_n||op<=||R_n||op<=3/2`.

It follows that each mean coordinate under F is still exactly b_(n,l)
and

`Var(m_(n,l)(F)) <= (3/2) M^2/n`.

Thus under both spin laws all fixed additive statistics concentrate
at the SAME vector b_n. This statement is uniform in the coefficient
arrays, and requires neither independence of the rounded spins nor
an asymptotic central limit theorem.

## 3. Their normalized energies remain separated

Gaussian rounding and flatness give the exact identity

```math
\mathbb E E_n(F)
=\frac{n-1}{\pi\sqrt n}\arcsin\frac1{2L\sqrt n}
\ge\frac{n-1}{2\pi L n}.                            (3)
```

Indeed every off-diagonal product
`B_ij arcsin(B_ij/(2L))` is the same positive number, and diagonal
terms of B vanish. The argument of arcsin lies in [-1,1] because
R is a covariance matrix. This calculation tests an actual Boolean
law on the actual signing, not a weighted surrogate.

Let `D_n=sup_x |E_n(x)-g(m_n(x))|`. Comparing expectations under the
two laws gives

`2 D_n >= E E_n(F)-|E g(m_n(F))-E g(m_n(X))|`.

For fixed k and continuous g, the absolute difference tends to zero
by the two variance bounds and uniform continuity on the fixed compact
cube. This holds even if b_n does not converge, proving (1). If g is
K-Lipschitz, the same difference is bounded by
`KM(1+sqrt(3/2))sqrt(k/n)`, which proves (2).

## 4. Exact primary-theorem implication and exclusions

[Guerra--Toninelli, The infinite volume limit in generalized mean field
disordered models](https://arxiv.org/pdf/cond-mat/0208579), equations
(4)--(6) and Theorem 1, requires the normalized deterministic mean to
be a fixed C1 function of finitely many bounded order parameters whose
unnormalized values are additive in system size, with a uniform O(1/n)
error. Equations (7)--(9) impose analogous additive-overlap and convex
covariance conditions. I read these exact hypotheses directly in the
primary paper.

Our theorem rules out the stated deterministic-mean representation
for bounded additive single-spin statistics already for EVERY
bounded-operator flat-signing sequence. Gaussian smoothing changes
the covariance but not this mean. Multiplying the mean by a fixed
positive a scales the obstruction in (1) by a. Including the global
orientation does not evade it: restrict to its positive branch.

This does not rule out a representation using order parameters whose
number grows proportionally to n, unbounded coefficients or derivatives,
nonadditive matrix-dependent statistics, or a substantially different
interpolation theorem. It also does not obstruct Ghirlanda--Guerra
identities enforced by small generic perturbations; such identities
are distinct from this finite-dimensional mean hypothesis. No claim
is made that the optimized sequence lacks a limit.
