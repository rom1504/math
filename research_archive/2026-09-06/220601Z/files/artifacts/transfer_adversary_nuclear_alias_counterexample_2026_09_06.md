# Why small cuts and two-root local limits do not remove the alias condition

2026-09-06. Proved counterexample in a finite-color Boolean polynomial
model. It verifies the necessity of the explicit same-degree alias
condition in the dual-split theorem. It is NOT a counterexample to that
theorem with its hypotheses, nor an extremal-signing construction.

Let n tend through powers of two and let H be the NORMALIZED symmetric
Sylvester Hadamard matrix, so H^2=I and |H_ij|=n^(-1/2). Fix an odd P>1.
Take P independent colors of Boolean seeds and put

```
eta_a=product_(c=1)^P S_(a,c),   X=H eta.
```

The eta_a are independent signs, while X has exact original seed degree
P. Every proper fixed-root cut of its kernel is O(n^-1/2), and all its
global coefficient cuts are bounded uniformly. Its covariance is I.

Use k=3 and its exact squarefree main

```
U_i=(X_i^3-(3-2/n)X_i)/sqrt(6),       Y=H U.
```

The linear correction removes exactly the repeated eta labels in the
cubic. Thus U and Y have exact original degree 3P. The two-factor
transport theorem gives every proper fixed-root cut of Y the bound
O(n^-1/2), and its global cuts remain uniformly bounded.

The covariance of U is EXACTLY

```
Cov(U)=v_n I,       v_n=(n-1)(n-2)/n^2.
```

To check this, its diagonal is `6*binom(n,3)/n^3`. For distinct rows
i,j, the n signs `n H_ia H_ja` are balanced. The third elementary
symmetric sum of those signs is zero, since their generating polynomial
is `(1-t^2)^(n/2)`. This gives the claimed off-diagonal zeros.
Consequently `Cov(Y)=v_n I` and

```
Cov(Y,U)=v_n H.                                         (1)
```

Let `C_i=tanh(X_i+Y_i)`. Every fixed collection (X_i,X_j,Y_i) has the
usual uniform Gaussian local comparison: all its nontrivial proper
contractions vanish; different original degrees have zero covariance.
The comparison keeps the exact covariance of X_i,X_j. In particular
the natural two-root mixed prediction is a DIAGONAL matrix D, with

```
d_i=E[tanh'''(Y_i+N)]/sqrt(6),    ||d||infinity<=2/sqrt(6).
```

Uniformly over i,j,
`E[C_i U_j]-d_i 1_{i=j} -> 0`. This is a genuine local conclusion,
not nuclear convergence.

Indeed dual-test the difference with the norm-one matrix H. Exact
symmetry and (1) give

```
sum_(i,j) H_ij E[C_i U_j]/n = average_i E[C_i Y_i]
  -> E[tanh(G_1+G_2)G_2]
   = E[sech^2(G_1+G_2)] > 0,                            (2)
```

where G_1,G_2 are independent standard Gaussians. The local comparison
justifies the linear-growth test by fixed-degree moment bounds. Meanwhile
`|tr(H D)|/n <= ||d||infinity/sqrt(n) -> 0`. Therefore

```
liminf_n ||Cov(C,U)-D||_*/n
 >= E[sech^2(G_1+G_2)] > 0.                             (3)
```

Here the literal remainder is R=Y. It has no degree-P part and has
small proper cuts at its higher degree 3P. All the other elementary
cut and covariance premises of the dual-split construction hold.
What fails is precisely its alias condition:

```
sum_j |E[R_i U_j]|^2 = v_n^2 -> 1,                      (4)
```

not zero. Y is a LATER noise formed by transporting the probe source
itself. Its old source contains X, whose degree is P. Hence it does
not meet the strict same-old-frame primitive-degree bound M<P used
by the Hall proof of the alias condition.

This example shows that neither bounded global cuts, small proper cuts,
nor uniform two-root bounded-test limits alone imply the desired
nuclear statement. The extra source/alias requirement is substantive.
It does not show that the actual same-source rich-frame extension fails.

Finite replay: `computations/transfer_adversary_nuclear_alias_counterexample_2026_09_06.py`
enumerates all 16, 256, and 65536 eta patterns at n=4,8,16. The exact
covariance formulas pass to floating-point tolerance. The H-dual cross
is approximately 0.214262, 0.283866, 0.390130, respectively. Subtracting
the universal diagonal bound 2/sqrt(6n) gives approximately -0.193987,
-0.004809, and +0.186006. These are finite numerical diagnostics, not
interval-certified transcendental values; the asymptotic strict
counterexample is proved analytically in (2)--(4).
