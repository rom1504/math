# Optimizing the smooth first mask in the second-response theorem

Date: 2026-09-05. Status: independently audited analytic theorem and exact
rational numerical certificate. The director and two independent agents
checked both the two-site theorem and this extension. See the two custom
response audit files linked from `fresh_limit_custom_response_2026_09_05.md`.

## 1. Scale-invariant form of the response formula

The two-site theorem in `fresh_limit_custom_response_2026_09_05.md` applies to
every fixed bounded smooth even `h`, not only to indicator masks. If its final
outer response is a hard threshold with own-spin bias `α`, its limiting
oriented half-energy is

\[
 {2\phi(\alpha/\sqrt v)\over\sqrt v}
      \mathbb E[h(G)\mathbf1_{|Y|\le\alpha}],
 \qquad v=\mathbb Eh(G)^2,
 \quad \operatorname{Cov}(G,Y)=\mathbb Eh(G),
 \quad \operatorname{Var}Y=v.                         \tag{1}
\]

Scaling `h` and `α` by the same positive constant leaves (1) unchanged.
Normalize `Eh²=1` and write `ρ=Eh∈(0,1)`. Then
`Y=ρG+sqrt(1-ρ²)Z`, where `G,Z` are independent standard normals, and (1) is

\[
 d(h,\alpha)=2\phi(\alpha)\,\mathbb E[h(G)k(G)],\quad
 k(g)=\Phi\!\left({\alpha-\rho g\over\sqrt{1-\rho^2}}\right)
      -\Phi\!\left({-\alpha-\rho g\over\sqrt{1-\rho^2}}\right).
                                                               \tag{2}
\]

Here `k` is bounded, smooth, and even. Set

\[
 p=\mathbb Ek(G)=2\Phi(\alpha)-1,\qquad
 s^2=\mathbb E(k(G)-p)^2>0.
\]

## 2. Exact Hilbert-space optimizer at fixed `(ρ,α)`

Among all square-integrable even `h` with `Eh=ρ` and `Eh²=1`,

\[
 \mathbb E hk
 =\rho p+\mathbb E[(h-\rho)(k-p)]
 \le\rho p+\sqrt{1-\rho^2}\,s.                        \tag{3}
\]

Equality is attained by the explicit bounded smooth even function

\[
 \boxed{\quad h_*(g)=\rho+
        {\sqrt{1-\rho^2}\over s}(k(g)-p).\quad}        \tag{4}
\]

Its derivatives of every fixed order are bounded. Thus it belongs to the
actual class required by the two-site proof, not merely to an abstract
relaxation. Intermediate fields need not be Boolean; if a magnitude bound of
one is desired, rescale `h*` and the outer threshold together.

Consequently the two-site theorem gives the explicit universal lower bound

\[
 \liminf_n {M_n\over n^{3/2}}\ge
 D(\rho,\alpha):=2\phi(\alpha)
       [\rho p+\sqrt{1-\rho^2}\,s].                    \tag{5}
\]

This optimizes `h` at fixed `ρ,α`. It does not assert that a central hard
threshold is globally optimal among all possible outer even masks.

Numerically, optimizing the two remaining scalar parameters gives
`ρ≈0.9400063157674`, `α≈0.807232538311`, and
`D≈0.38578983745603`. The simple rational choice `ρ=47/50`, `α=81/100` gives
approximately `0.38578587690878`.

## 3. A positive series for exact certification

For the probabilists' Hermite polynomials,

\[
 \mathbb E[\mathbf1_{|Z|\le\alpha}H_{2j}(Z)]
                 =-2\phi(\alpha)H_{2j-1}(\alpha).
\]

This follows by integrating `(φH_(2j-1))'=-φH_(2j)`. Conditional Gaussian
expectation multiplies degree `r` by `ρ^r`. Parseval therefore gives

\[
 \boxed{\quad
 s^2=4\phi(\alpha)^2\sum_{j=1}^{\infty}
       {\rho^{4j}H_{2j-1}(\alpha)^2\over(2j)!}.
 \quad}                                              \tag{6}
\]

Every summand is nonnegative. If `s_K²` is the sum through `K`, then

\[
 s_K^2\le s^2\le s_K^2+
               \rho^{4(K+1)}p(1-p).                  \tag{7}
\]

The upper bound uses the full unweighted Hermite variance `p(1-p)` and the
largest remaining geometric multiplier. At rational `ρ,α`, the polynomial
factor in every finite summand is rational. Thus the exact interval machinery
already used for `φ,Φ` provides a short rigorous certificate, with no
two-dimensional numerical quadrature and no inferred sign of a truncation
error.

The exact-arithmetic script
`computations/fresh_limit_response_variational_certificate.py` uses `K=100`
and the rational choices `ρ=47/50`, `α=81/100`. Its stdout is saved in
`computations/results/fresh_limit_response_variational_certificate.json`.
It proves

\[
 \begin{split}
 0.385785876908778466066127787300791530025224664775532506435400
 &\le D(47/50,81/100)\\
 &\le
 0.385785876909691210517840971183287938897036818354614542045282.
 \end{split}
\]

In particular the value is rigorously greater than `0.38578`. This certifies
the scalar consequence, not independently the preceding two-site theorem.
