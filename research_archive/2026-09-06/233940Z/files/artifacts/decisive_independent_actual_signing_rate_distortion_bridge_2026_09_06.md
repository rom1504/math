# A finite actual-signing rate-distortion constraint, and the remaining gap

Status: elementary extension of the director's audited precision argument;
pending independent audit. It is NOT an unrestricted cap lower bound.

Let `J_lambda(X)=inf_L[I(X;L)+lambda E Var(X|L)]` as in the director's
precision-Schur theorem. For ANY invertible real T, not only orthogonal T,
put U=TX and factor

```math
T^T\operatorname{diag}(\lambda_i)T
=R^T\operatorname{diag}(p_i)R
```

with R unit upper triangular and p_i>0. Then

```math
\sum_i J_{p_i}(X_i)
\le\sum_i J_{\lambda_i}(U_i)+\operatorname{TC}(X).       \tag{1}
```

The proof is the same completed-square/channel proof, with reconstructed
estimate `T^(-1) E[U_i|M_i]`. Parent coordinate i receives label
`(X_(i+1),...,X_n,M_1,...,M_n)`. Weighted errors sum exactly before replacing
estimators by conditional means, and the information chain is invariant under
the invertible map T. Orthogonality was not used in this part; it is needed
for the later supporting-constant majorization in the director's E inequality.

## Exact specialization to every hollow signing

Let A be any hollow symmetric signing, `B=A/sqrt(n-1)`, and take independent
Rademacher X_i. For a real u with I+uB invertible put T=I+uB. Every coordinate
U_i has EXACTLY the same law

```math
W_{n,u}=\varepsilon+{u\over\sqrt{n-1}}
\sum_{j=1}^{n-1}\varepsilon_j,
```

where all displayed signs are independent. Hollowing makes the first term
independent of the sum, and coefficient signs disappear by Rademacher symmetry.
Thus if p_i are precision Schur pivots of `lambda(I+uB)^2`, (1) gives

```math
{1\over n}\sum_i J_{p_i}(\varepsilon)
\le J_\lambda(W_{n,u}).                                  \tag{2}
```

This is an actual-signing statement with a universal one-dimensional right
side, not a statement only about a random-sign or orthogonal ensemble.
For fixed u,lambda the right side tends to
`J_lambda(epsilon+uG)`. Indeed W_2 convergence follows from the CLT plus
second moments, and channel transport gives W_2 continuity of J_lambda on
uniformly second-moment-bounded source families:

```math
J_\lambda(Y)\le J_\lambda(X)
+2\lambda\sqrt{\operatorname{Var}X}\,W_2(X,Y)
+\lambda W_2(X,Y)^2,
```

and conversely with X,Y interchanged. Translation may be matched first.

## A Gibbs specialization and a diagnostic barrier

For any globally reversal-invariant law of cube X, its one-coordinate
marginals are Rademacher and `TC(X)=n log2-H(X)`. Choosing constant child
channels in (1) gives

```math
H(X)\le n\log2+\mathbb E X^TQX-\sum_i J_{p_i}(\varepsilon),
\qquad Q=T^T\operatorname{diag}(\lambda_i)T.             \tag{3}
```

If Q=lambda(I-uB)>0 and X is the Gibbs law proportional to
`exp(beta q_A(X)/sqrt(n-1))`, choose beta=2lambda u. The energy term cancels
against the Gibbs entropy identity, yielding

```math
{1\over n}\log\sum_xe^{\beta q_A(x)/\sqrt{n-1}}
\le\log2+\lambda-{1\over n}\sum_iJ_{p_i}(\varepsilon).    \tag{4}
```

This is an UPPER pressure bound. It is not the universal LOWER comparison
needed to identify the original minimum with the recursive construction.
Its positive-definiteness requirement reintroduces spectral/diagonal-majorant
constraints. Dropping pivot losses, assuming orthogonality of B, or reversing
(4) would be invalid.

## Audit judgment

The new exact scalar Bellman identity H=E supplies a genuine analytic closure
inside the recursive construction. Equations (1)--(2) provide a concrete first
map from arbitrary signings into its scalar rate-distortion language. They do
not yet involve the Boolean cap. The direct Gibbs attempt (3)--(4) has the
wrong comparison direction for a universal lower theorem and retains the
known spectral barrier. No construction-optimality assertion follows.
