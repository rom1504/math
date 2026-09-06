# Precision Schur factorization closes temperature alignment

2026-09-06. Director theorem. Two independent researchers reconstructed
the binary argument immediately; detailed audits are preserved separately.
This proves an analytic statement previously left open in the recursive
upper construction. It is not a proof that this ensemble is optimal among
all signings, and is not a convergence theorem for M_n.

## 1. Definitions and exact scalar duality

Fix t>0 and let g_t(v) be minus one half the entropic quadratic
self-transport cost of a centered Gaussian of variance v. Direct Gaussian
optimization gives, with rho in [0,1),

```math
g_t(v)=-tv(1-\rho)+\tfrac14\log(1-\rho^2),
\qquad 2tv=\frac{\rho}{1-\rho^2}.
```

Since g'_t(v)=-t(1-rho), its exact supporting-line representation is

```math
g_t(v)=\sup_{0<\lambda\le t}\{c_t(\lambda)-\lambda v\},
\qquad c_t(\lambda)=\tfrac14\log\frac{\lambda(2t-\lambda)}{t^2}.
```

For a finite-second-moment real random variable X define

```math
J_\lambda(X)=\inf_L\{I(X;L)+\lambda\mathbb E\operatorname{Var}(X\mid L)\},
\qquad
E_t(X)=\sup_L\{g_t(\mathbb E\operatorname{Var}(X\mid L))-I(X;L)\}.
```

Swapping two suprema, not a minimax interchange, gives

```math
E_t(X)=\sup_{0<\lambda\le t}\{c_t(\lambda)-J_\lambda(X)\}.
```

All channels may be approximated in objective value; no optimizing lambda
or channel is assumed to exist. Conditional means are square integrable.
The constant channel shows finiteness of the needed infima.

## 2. Two-input precision inequality

Let (A,B) have finite second moments and finite mutual information, and put
U=(A+B)/sqrt(2), V=(A-B)/sqrt(2). For lambda_1,lambda_2 in (0,t], set

```math
a=(\lambda_1+\lambda_2)/2,\quad
h=2\lambda_1\lambda_2/(\lambda_1+\lambda_2),\quad
d=(\lambda_1-\lambda_2)/(\lambda_1+\lambda_2).
```

Then

```math
J_a(A)+J_h(B)\le J_{\lambda_1}(U)+J_{\lambda_2}(V)+I(A;B). \tag{1}
```

Choose child channels M|U and N|V, independently conditional on (U,V).
Use child conditional means uhat=E[U|M], vhat=E[V|N], and rotate them
back to ahat,bhat. The exact identity is

```math
\lambda_1(U-\widehat U)^2+\lambda_2(V-\widehat V)^2
=a[A-\widehat A+d(B-\widehat B)]^2+h(B-\widehat B)^2.
```

For A the parent label (B,M,N) permits the estimator
ahat-d(B-bhat); for B the label (M,N) permits bhat. Conditional means
can only lower these weighted squared errors. Their information cost is
exactly

```math
I(A;B,M,N)+I(B;M,N)
=I(U;M)+I(V;N)+I(A;B)-I(M;N).
```

Dropping the last nonnegative mutual information and taking approximate
child infima proves (1). In particular no Gaussian conditional laws,
equal child temperatures, or covariance regularity is assumed.

The precisions a,h stay in (0,t], have product lambda_1 lambda_2, and
their logarithms are less spread than those of lambda_1,lambda_2.
Moreover

```math
\frac{d^2}{ds^2}c_t(e^s)
=-\frac{t e^s}{2(2t-e^s)^2}<0.
```

Hence c_t(a)+c_t(h)>=c_t(lambda_1)+c_t(lambda_2). Combining with (1)
and then taking the two independent child suprema yields

```math
E_t(U)+E_t(V)-I(A;B)\le E_t(A)+E_t(B).                 \tag{2}
```

The parent precisions need not coincide: they are two admissible
supporting lines evaluated at the same parent law when A,B are identically
distributed. Requiring them to coincide was an unnecessary proof obligation.

## 3. Arbitrary orthogonal mixing

More generally, if A=(A_1,...,A_d), U=OA for orthogonal O, and the
total correlation TC(A)=D(law(A)||product_i law(A_i)) is finite, then

```math
\sum_i E_t(U_i)-\operatorname{TC}(A)\le\sum_i E_t(A_i). \tag{3}
```

Choose positive child precisions lambda_i<=t, and factor the positive
definite precision Q=O^T diag(lambda_i) O as R^T diag(p_i) R, where R
is unit upper triangular. Its Schur pivots satisfy 0<p_i<=t. Complete
the square row by row, using parent label (A_{i+1},...,A_d,M_1,...,M_d)
for coordinate i. The error identity and conditional-mean optimality
give parent cost at most child cost. The information chain gives exactly

```math
\sum_i I(A_i;A_{>i},M)
=\operatorname{TC}(A)+\sum_i I(U_i;M_i)-\operatorname{TC}(M).
```

Here the M_i are conditionally independent local child channels.
Thus sum_i J_{p_i}(A_i)<=sum_i J_{lambda_i}(U_i)+TC(A).

For completeness, log(p) is majorized by log(lambda). Write Q=S^T S
with upper triangular S having diagonal sqrt(p_i). For any index set I,
Cauchy--Binet gives det(Q_{II})>=product_{i in I} p_i by selecting the
row set I. The variational eigenvalue bound gives
det(Q_{II})<=product_{j<=|I|} lambda_j^downarrow. Full products agree
by determinant. These are exactly the required partial-product inequalities.
Concavity of c_t(exp(s)) gives sum_i c_t(p_i)>=sum_i c_t(lambda_i).
Taking separate child suprema proves (3). The proof is finite-dimensional;
d is not sent to infinity inside an uncontrolled estimate.

## 4. Bellman consequence and its boundary

For the repository's binary Bellman operator with symmetric parent law nu,
safe simultaneous reversal and input interchange produce equal signed
parent marginals without changing child absolute laws, while only decreasing
the entropy charge. Equation (2) therefore proves

```math
\mathcal B E_t(\nu)\le E_t(\nu).
```

Together with the separately proved lower bound E_t<=H and the existing
Gaussian-boundary identification H=lim_r B^r g_t, monotonicity gives
H<=E_t (since g_t<=E_t), hence **H=E_t**. Those boundary/lower-envelope
results must be reconstructed before using this consequence in a cap bound;
they are not proved merely by the Schur inequality above.

This identifies the exact deep Gaussian-boundary Bellman certificate.
It does not assert that the row-permanent/Finner bound is tight for the
actual ensemble pressure or cap, nor identify the unrestricted M_n with it.
The remaining original-problem gap is an actual-signing lower comparison or
a separate limit theorem, not temperature alignment.
