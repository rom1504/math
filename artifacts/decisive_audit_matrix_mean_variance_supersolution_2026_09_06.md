# Independent audit: noncommuting matrix precision alignment

Date: 2026-09-06. Status: PASS for the finite-dimensional vector
supersolution below. The director proposed the block Schur argument.
This note does not silently import a vector stopped-tree theorem, and
does not assert a new sign-matrix construction.

## 1. The matrix envelope is exact

For a finite `R^k`-valued source `X`, define

```math
E_t^{(k)}(X)=\sup_L\{\operatorname{Tr}g_t(\mathbb E\operatorname{Cov}(X\mid L))-I(X;L)\}.
```

For `0<Lambda<=tI`, set

```math
C_t(Lambda)=\tfrac14\log\det[Lambda(2tI-Lambda)/t^2],
\quad J_Lambda(X)=\inf_L\{I(X;L)+\mathbb E[(X-\mathbb E[X\mid L])^T
Lambda(X-\mathbb E[X\mid L])]\}.
```

Then exactly

```math
E_t^{(k)}(X)=\sup_{0<Lambda\le tI}\{C_t(Lambda)-J_Lambda(X)\}.       \tag{1}
```

To check the spectral representation, diagonalize a positive semidefinite
matrix `V`. Replacing `Lambda` by its diagonal pinching in this basis
preserves `Tr(Lambda V)` and increases BOTH `log det Lambda` and
`log det(2tI-Lambda)` by Hadamard's determinant inequality. The pinched
matrix is still in `(0,tI]`. Optimizing its diagonal entries separately
is precisely the scalar supporting-line formula for `sum_i g_t(V_ii)`.
This also covers zero covariance eigenvalues, whose optimum precision
is `t`. Interchanging the two suprema over `L` and `Lambda` gives (1).

## 2. Noncommuting square completion

Let `A,B` be jointly finite `R^k`-valued sources and put
`U=(A+B)/sqrt(2)`, `V=(A-B)/sqrt(2)`. Choose arbitrary child precisions
`Lambda_1,Lambda_2` in `(0,tI]`. They need NOT commute. Put

```math
P=(Lambda_1+Lambda_2)/2,\quad D=(Lambda_1-Lambda_2)/2,
\quad S=P-DP^{-1}D.
```

Both `P,S` are positive definite and at most `tI`. Moreover

```math
S=Lambda_1 P^{-1}Lambda_2
 =2(Lambda_1^{-1}+Lambda_2^{-1})^{-1}.                    \tag{2}
```

The nonsymmetric-looking product in the middle is symmetric because
it equals the displayed Schur complement.

For arbitrary errors `e_A,e_B`, rotated errors satisfy the exact identity

```math
e_U^T Lambda_1 e_U+e_V^T Lambda_2 e_V
=(e_A+P^{-1}D e_B)^T P(e_A+P^{-1}D e_B)+e_B^T S e_B.       \tag{3}
```

Use conditionally independent child channels `M|U,N|V`, transform their
conditional-mean reconstructions back to `widehat A,widehat B`, and use
parent estimators `widehat A-P^{-1}D(B-widehat B)` with label `(B,M,N)`
and `widehat B` with label `(M,N)`. The conditional mean minimizes every
fixed positive quadratic loss. Therefore (3) and the unchanged chain-rule
information identity give

```math
J_P(A)+J_S(B)\le J_{Lambda_1}(U)+J_{Lambda_2}(V)+I(A;B).   \tag{4}
```

Arbitrarily accurate child channels suffice, so label attainment is not
an extra hypothesis.

## 3. Why the block pivots are log-majorized

The positive definite block matrix

```math
Q=\begin{pmatrix}P&D\\D&P\end{pmatrix}
=\begin{pmatrix}I&0\\DP^{-1}&I\end{pmatrix}
 \begin{pmatrix}P&0\\0&S\end{pmatrix}
 \begin{pmatrix}I&P^{-1}D\\0&I\end{pmatrix}               \tag{5}
```

is orthogonally similar to `diag(Lambda_1,Lambda_2)`.
Conjugate (5) by a block-diagonal orthogonal matrix which separately
diagonalizes `P` and `S`. Its middle factor becomes diagonal; its left
factor remains unit block-lower-triangular, hence unit scalar-lower-
triangular. Thus the scalar squared Cholesky pivots of the conjugated
matrix are exactly `eigen(P)` followed by `eigen(S)`.

Here is the needed scalar fact without assuming block commutation.
For a positive covariance matrix with scalar Cholesky pivots `delta_i`,
the product of any selected `r` pivots is at most the determinant of
the same principal `r`-coordinate submatrix: remove the unselected
preceding regressors. That principal determinant is at most the product
of the largest `r` eigenvalues. The total pivot product equals the full
determinant. Hence `log(delta)` is majorized by the log eigenvalue vector.

Since `c_t(exp(s))` is concave, this proves

```math
C_t(P)+C_t(S)\ge C_t(Lambda_1)+C_t(Lambda_2).              \tag{6}
```

Combining (1), (4), and (6), and optimizing the two child precisions,
proves the vector entropy inequality

```math
E_t^{(k)}(U)+E_t^{(k)}(V)-I(A;B)
\le E_t^{(k)}(A)+E_t^{(k)}(B).                            \tag{7}
```

Equal signed vector marginals therefore give the corresponding binary
Bellman supersolution. No scalarization of the conditional covariance
and no common eigenbasis of the child precisions is needed.

## 4. Precisely what is and is not imported

For a fixed dimension `k`, an independently established vector analogue
of the conditional-copy lower envelope and Gaussian stopped-tree boundary
theorem would give `E_t^(k)<=H_t^(k)`; (7), together with the trivial-label
bound, then gives equality. The ingredients appear compatible with a
fixed-dimensional vector extension, but that extension is NOT a logical
consequence of the block algebra alone and is not fully audited here.

The dimension is fixed in such a boundary passage. This note gives no
uniformity when the seed dimension grows with the eventual matrix order.
It also supplies no sign-flat realization carrying arbitrary finite
seed cap into the vector source. Both are separate original-problem
obligations.

## 5. Orthogonal invariance and exact independent tensor additivity

For possibly different vector dimensions `k,l` and jointly finite `X,Y`,

```math
E_t^{(k)}(X)+E_t^{(l)}(Y)
\le E_t^{(k+l)}(X,Y)
\le E_t^{(k)}(X)+E_t^{(l)}(Y)+I(X;Y).                    \tag{8}
```

For the first inequality, choose block-diagonal precision and independent
local channels `M|X,N|Y`. Their joint information is the sum of their
local informations minus `I(M;N)`; conditional reconstruction with the
refined label `(M,N)` only improves each block's quadratic error. Hence
the joint `J` is at most the sum of the local `J` values, and optimization
proves the lower bound. This argument does NOT require `X,Y` independent.

For the upper inequality, take an arbitrary joint precision
`Lambda=[[P,D],[D^T,Q]]` and a single arbitrary channel `L|(X,Y)`.
Use its reconstructed `X,Y` means in the block square completion, with
`S=Q-D^T P^{-1}D`. Parent labels `(Y,L)` and `L` give information sum
`I(X,Y;L)+I(X;Y)`. Their weighted reconstruction cost is at most the
original joint cost, and the same block-pivot proof gives
`C_t(P)+C_t(S)>=C_t(Lambda)`. Optimizing yields the upper bound in (8).

Therefore independent vector sources give EXACT tensor additivity.
Orthogonal invariance is immediate from the covariance definition (or
by conjugating the precision domain in (1)). Iterating (8) gives

```math
sum_i E_t(X_i) <= E_t^{(d)}(X)
              <= sum_i E_t(X_i)+TC(X),
```

so for every real orthogonal `O`,

```math
sum_i E_t((OX)_i)-TC(X) <= sum_i E_t(X_i).                \tag{9}
```

This is a finite exact, orthogonally invariant entropic tensor functional.
It is not a claim that an arbitrary signing's Boolean cap is determined
by that functional.

## 6. Noncommuting finite diagnostics

`computations/decisive_audit_matrix_precision_checks_2026_09_06.py`
uses the exact rational noncommuting pair
`Lambda_1=diag(1,1/3)`, `Lambda_2=[[2/3,1/6],[1/6,2/3]]` and checks
seven rational square identities plus the matrix harmonic formula.
It also passes 500 random positive-definite cases in dimensions
`1,2,3,4,8`, checking every partial log-majorization sum and the curvature
inequality. Maximum normalized harmonic error is `2.18e-15`, maximum
log-majorization numerical error is `3.20e-14`. These checks complement
the noncommutative proof; they are not optimizer or boundary calculations.
