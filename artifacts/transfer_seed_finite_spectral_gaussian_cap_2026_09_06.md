# A finite spectral Gaussian witness requiring only trace and diagonal control

Date: 2026-09-06. Seed-transfer track. This is an unconditional finite
inequality for actual hollow signings. The random-compression moment
hypotheses in Section 2 are a separate obligation, not assumed proved.
They have since been proved for sparse restrictions of every bounded-cap
parent in `transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`.
The conditional statement and its original diagnostic record below are
retained to make the dependency explicit.

## 1. Finite inequality, with no operator-norm assumption

Let `A` be a hollow symmetric sign matrix of order `n` and `L=A/sqrt(n)`.
Let `R` be any positive semidefinite matrix with `R_ii>=a>0`. Set

```math
v=\operatorname{tr}R/n,\qquad
\delta^2=\frac1n\sum_i(R_{ii}-v)^2,\qquad
s=\operatorname{tr}R^2/n.
```

Then

```math
\boxed{\quad
\frac{Q(A)}{n^{3/2}}
\ge\frac1\pi\left[
 \frac{\operatorname{tr}(LR)}{nv}
 -\frac{\delta\sqrt s}{a^2}
 -\frac{s}{a^2\sqrt n}\right].\quad}                       (1)
```

In particular no bound on `||L||op`, no entrywise approximation of a
covariance matrix, and no delocalization of individual eigenvectors is
required. The row/column squared norms of `L` are exactly `(n-1)/n`;
this sign-flat identity is what controls diagonal normalization.

Proof. Put `D=diag(R)` and `C=D^(-1/2) R D^(-1/2)`. Then `C` is a
correlation matrix. Gaussian hyperplane rounding produces Boolean spins
with correlation `(2/pi)arcsin(C_ij)`, so

```math
\frac{Q(A)}{n^{3/2}}
\ge\frac1{\pi n}\sum_{i,j}L_{ij}\arcsin(C_{ij}).             (2)
```

The positive Taylor coefficients imply
`|arcsin z-z|<=|z|^3` for `|z|<=1` (their sum after the linear term is
`pi/2-1<1`). Since `|L_ij|<=1/sqrt(n)` and `|C_ij|<=1`,

```math
\left|\frac1n\sum_{ij}L_{ij}(\arcsin C_{ij}-C_{ij})\right|
\le\frac{\operatorname{tr}C^2}{n\sqrt n}
\le\frac{s}{a^2\sqrt n}.                                  (3)
```

For the linear normalization error put
`E=D^(-1/2)-v^(-1/2)I`. The derivative bound for `t^(-1/2)` on
`[a,infinity)` gives `||E||F<=sqrt(n) delta/(2a^(3/2))`.
Decompose

```math
C-R/v=E R D^{-1/2}+v^{-1/2}R E.
```

The flat row/column norms give `||LE||F,||EL||F<=||E||F`.
Cyclically moving factors inside the trace and using Frobenius
Cauchy--Schwarz consequently yields

```math
|\operatorname{tr}L(C-R/v)|
\le2a^{-1/2}\|E\|_F\|R\|_F
\le n\delta\sqrt s/a^2.                                   (4)
```

Equations (2)--(4) prove (1).

## 2. A rational cubic polynomial makes semicircle-type moments sufficient

Fix

```math
P(x)=x^3+\frac85x^2-\frac25x-\frac35,\qquad
R_n=\frac1{10}I+P(L_n)^2.
```

Assume along an actual signing sequence that

```math
\frac1n\operatorname{tr}R_n\to\frac{361}{50},\qquad
\frac1n\operatorname{tr}(L_nR_n)\to\frac{288}{25},           (5)
```

and that

```math
\frac1n\operatorname{tr}R_n^2=O(1),\qquad
\frac1n\sum_i\left((R_n)_{ii}-\frac1n\operatorname{tr}R_n\right)^2
\longrightarrow0.                                        (6)
```

The finite inequality, with `a=1/10`, implies

```math
\liminf\frac{Q(A_n)}{n^{3/2}}
\ge\frac{576}{361\pi}
>\frac{2016}{3971}>\frac{507}{1000}>\frac12.                 (7)
```

The rational comparison uses only `pi<22/7`.

For clarity, the values in (5) are exactly the centered semicircle
moments: with orthonormal polynomials
`1,x,x^2-1,x^3-2x`, the coefficients of `P` are `1,8/5,8/5,1`.
Thus `E P(S)^2=178/25` and `E S P(S)^2=288/25` for a variance-one
semicircle variable `S`. Adding `I/10` gives `361/50`.

Sufficient elementary moment conditions are global trace convergence
to the semicircle moments through degree 7, bounded global moments
through degree 12, and mean-square concentration of the diagonal of
the fixed degree-6 polynomial `P(L)^2`. These conditions concern more
than the empirical eigenvalue distribution alone. In particular an
ESD theorem without diagonal control is not enough to invoke (7).

For a random principal compression, conditions (5)--(6) may instead
hold in probability; the identical argument then gives (7) in
probability. Establishing them uniformly for all vanishing retentions
of symmetric Hadamards remains separate work.

The exact polynomial coefficients, semicircle moments, and rational
comparison are reproduced in
`computations/transfer_seed_spectral_gaussian_verify_2026_09_06.py`.
Its optional finite compression experiment reports actual Gaussian
rounding expectations and moment diagnostics, not asymptotic proofs.

The default optional diagnostic was run with parent Walsh order `2^18`
and RNG seed `20260906`. At retained orders `64,128,256,512`, its
floating exact-Gaussian expected normalized energies were respectively
`0.4668181470859914`, `0.4944349520971711`,
`0.5107060789611981`, `0.5041652982289402`.
These are reproducible finite lower-witness diagnostics, not exact cap
certificates or asymptotic claims. The source also checks each displayed
normalization and arcsine-error inequality on those instances.
