# Independent proof audit: the mean-variance envelope is a supersolution

Date: 2026-09-06. Status: PASS. The director supplied the precision-side
Schur factorization; the following is an independent reconstruction with
the original definition of `E_t` checked against
`continued_convergence_temperature_alignment_next_target_2026_09_06.md`.
This closes the old analytic supersolution question. It does not identify
the original signing minimax value or prove its convergence.

## 1. Exact envelope representation

Fix `t>0` and put, for `0<lambda<=t`,

```math
c_t(lambda)=\frac14\log\frac{lambda(2t-lambda)}{t^2},\qquad
J_lambda(X)=\inf_L\{I(X;L)+lambda\,\mathbb E\operatorname{Var}(X\mid L)\}.
```

The Gaussian potential has the exact supporting-line representation

```math
g_t(v)=\sup_{0<lambda\le t}\{c_t(lambda)-lambda v\}.
```

Indeed the stationary parameter is `lambda=t(1-rho)`, where
`2tv=rho/(1-rho²)`; substitution gives the archived formula for `g_t`.
At `v=0`, `lambda=t` is allowed and gives zero. Consequently

```math
E_t(X)=\sup_L\{g_t(\mathbb E\operatorname{Var}(X\mid L))-I(X;L)\}
      =\sup_{0<lambda\le t}\{c_t(lambda)-J_lambda(X)\}.       \tag{1}
```

This interchanges two suprema, not a supremum with an infimum.

## 2. The weighted reconstruction calculation

Let `(A,B)` have any finite joint law, and define
`U=(A+B)/sqrt(2)`, `V=(A-B)/sqrt(2)`. Choose arbitrary positive precisions
`lambda_1,lambda_2<=t` and child channels `M|U`, `N|V`, sampled
independently conditional on `(U,V)`. Write

```math
\widehat U=\mathbb E[U\mid M],\quad
\widehat V=\mathbb E[V\mid N],\quad
\widehat A=(\widehat U+\widehat V)/\sqrt2,\quad
\widehat B=(\widehat U-\widehat V)/\sqrt2.
```

Set `a=(lambda_1+lambda_2)/2`,
`h=2lambda_1lambda_2/(lambda_1+lambda_2)`, and
`d=(lambda_1-lambda_2)/(lambda_1+lambda_2)`. Direct expansion gives

```math
lambda_1(U-\widehat U)^2+lambda_2(V-\widehat V)^2
=a[A-\widehat A+d(B-\widehat B)]^2+h(B-\widehat B)^2.       \tag{2}
```

Give `A` the label `(B,M,N)` and estimator
`widehat A-d(B-widehat B)`, and give `B` the label `(M,N)` and estimator
`widehat B`. Conditional expectations minimize squared error, so

```math
a\,\mathbb E\operatorname{Var}(A\mid B,M,N)
+h\,\mathbb E\operatorname{Var}(B\mid M,N)
\le lambda_1\,\mathbb E\operatorname{Var}(U\mid M)
   +lambda_2\,\mathbb E\operatorname{Var}(V\mid N).         \tag{3}
```

No independence of posterior residuals, linear optimal predictor, or
Gaussian conditional law is used.

The exact information identity is

```math
I(A;B,M,N)+I(B;M,N)
=I(A;B)+I(U;M)+I(V;N)-I(M;N).                             \tag{4}
```

The last term is nonnegative. Combining (3)--(4), and taking independently
arbitrarily accurate child channel infima, proves

```math
J_a(A)+J_h(B)\le J_{lambda_1}(U)+J_{lambda_2}(V)+I(A;B).   \tag{5}
```

The channels may be approximating channels; attainment is unnecessary.

## 3. Precision curvature has the correct direction

Both `a,h` lie in `(0,t]`, and `ah=lambda_1lambda_2`. The pair
`(log a,log h)` is majorized by `(log lambda_1,log lambda_2)`: `a` and
`h` both lie between the two child precisions, and their products agree.
Furthermore

```math
\frac{d^2}{ds^2}c_t(e^s)
=-\frac{t e^s}{2(2t-e^s)^2}<0.
```

Thus

```math
c_t(a)+c_t(h)\ge c_t(lambda_1)+c_t(lambda_2).              \tag{6}
```

Using (1), (5), and (6), then taking the two precision suprema, yields

```math
E_t(U)+E_t(V)-I(A;B)\le E_t(A)+E_t(B).                    \tag{7}
```

For actual equal signed parent marginals `nu`, divide by 2 and optimize
the admissible Bellman coupling. This proves `B E_t(nu)<=E_t(nu)`.
If the operator is initially written with averaged absolute marginals,
the archived safe reversal/input-swap reduction must first be applied;
(7) itself does not silently assume that reduction for an arbitrary
different operator.

## 4. Scope of the resulting closure

The older stopped-tree theorem gives `E_t<=H_t`, where `H_t` is the
Gaussian-boundary Bellman fixed point. Since `G_t<=E_t` and the present
supersolution bounds every lower iterate, `H_t<=E_t`. Therefore the
already established Bellman limit is exactly `H_t=E_t`.

This last identification imports the stopped-tree boundary theorem, not
just the algebra above. It does not imply `T_t=E_t`: an envelope averaging
the potential over individual posterior variances can still be strictly
larger. Nor does it bridge the construction's entropy/Finner certificate
to the minimax over all finite signings. Those distinctions are essential.

## 5. Independent finite diagnostics

`computations/decisive_audit_mean_variance_supersolution_checks_2026_09_06.py`
passes 400 exact rational square factorizations and 500 finite-channel
instances with arbitrary parent marginals, degenerate supports and child
precision ratios down to `exp(-10)`. The largest information-identity
error is `3.56e-15`; the smallest reconstruction slack is zero up to
`4e-18` roundoff. These are diagnostics of the elementary steps, not
numerical optimization of `J` or a substitute for the proof.
