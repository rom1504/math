# Exact cubic response of the finite-anchor certificate

Date: 2026-09-05. Status: proved closed formula and rational interval
certificate; independently checked by the variational and literature agents.

This supplies a nonzero direction required by the variance-normalized
unmarked transport theorem. It is not a floating-point inference about
a Gaussian coefficient.

Let `(G_0,V,W)` be the centered Gaussian triple from the banked finite-anchor
fixed-point construction. Its covariance is

```
Var(G_0)=Var(V)=1,  Var(W)=p,
Cov(G_0,V)=rho, Cov(V,W)=w, Cov(G_0,W)=p.
```

Here `rho=2027/2500`, `alpha=361/500`; the exact intervals for `p,w` are
provided by `fresh_finite_anchor_fixed_point_certificate.py`. Define

```
F=sign(W) 1{|V|>alpha},  h3(g)=(g^3-3g)/sqrt(6),
sigma=sqrt(p-w^2), c=w/sigma, k=(p-rho*w)/sigma.
```

Conditioning on `V,Z`, where `W=w V+sigma Z` and `V,Z` are independent
standard Gaussians, gives

```
G_0 = rho V + k Z + independent Gaussian residual,
Var(residual)=1-rho^2-k^2 > 0.
```

The third Hermite polynomial consequently has conditional expectation

```
E[He_3(G_0)|V,Z]
 = rho^3 He_3(V) + 3 rho^2 k He_2(V) Z
   + 3 rho k^2 V He_2(Z) + k^3 He_3(Z).
```

Write `s(v)=2 Phi(cv)-1`. The signed Gaussian moments against
`sign(Z+cv)` of `1,Z,He_2(Z),He_3(Z)` are respectively

```
s(v), 2 phi(cv), -2cv phi(cv), 2(c^2 v^2-1) phi(cv).
```

The negative and positive `V` tails give equal contributions. Set

```
d=sqrt(1+c^2), a=alpha*d,
J0=phi(0) [1-Phi(a)]/d,
J2=phi(0) [a phi(a)+1-Phi(a)]/d^3,
A=6 rho^2 k-6 rho k^2 c+2 k^3 c^2,
D=-6 rho^2 k-2 k^3.
```

These are exactly the integrals of `phi(v)phi(cv)` and
`v^2 phi(v)phi(cv)` over `v>=alpha`. Integration by parts using
`(phi He_2)'=-phi He_3` yields the closed formula

```
E[F h3(G_0)] = (2/sqrt(6)) {
 rho^3 phi(alpha)(alpha^2-1)[2 Phi(c alpha)-1]
 +(2c rho^3+A)J2 + (D-2c rho^3)J0 }.
```

The executable certificate uses outward rational intervals throughout,
including the previously audited finite-series Gaussian CDF and density
bounds. No numerical quadrature or optimization enters the calculation.
It proves

```
-0.063914331766535141564927573349386587491619191945279997110215
 <= E[F h3(G_0)] <=
-0.063914331766535141564927573349386587491619191945279997110143.
```

It also certifies

```
Var(V|W) > 0.07385942754623718,
Var(G_0|V,W) > 0.30532226606700294.
```

Thus the cubic response is bounded away from zero and the Gaussian event
`H=0, |W|<=t/4` has probability bounded below by a positive constant times
`t` for all sufficiently small `t`. The latter follows directly from the
nondegenerate bivariate Gaussian density of `(V,W)`, for instance by
integrating over any fixed interval strictly outside `[-alpha,alpha]`.
Both properties persist uniformly in a sufficiently small covariance
neighborhood. This is the uniformity needed to choose finite smooth tree
approximations after fixing the operator bound and perturbation scale.

Reproduction:

```
.venv/bin/python -B computations/fresh_finite_anchor_fixed_point_certificate.py
.venv/bin/python -B computations/fresh_finite_anchor_cubic_coefficient_certificate.py
```

The second output is saved in
`computations/results/fresh_finite_anchor_cubic_coefficient_certificate.json`.
