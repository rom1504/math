# Exact Gaussian tree fixed points: the critical criterion

Date: 2026-09-05. Status: independent proof. The parent independently found
the same Galton--Watson interpretation. This note includes the critical
case, not just the strict contraction case used by the numerical certificate.

## 1. Setup

Index independent standard Gaussians `Z_T` by finite rooted trees whose
external root has degree one and whose other vertices have odd degree.
If `T` has child-branch multiplicities `m_tau`, define

\[
h_T(Z)=\prod_\tau H_{m_\tau}(Z_\tau)/\sqrt{m_\tau!}.
\]

The total child multiplicity is even. Every finite even Hermite multiindex
defines exactly one such parent tree. Therefore the functions `h_T` form
the complete jointly-even Hermite basis, and

\[
\mathcal U H=\sum_T E[Hh_T]Z_T
\]

is an isometry from jointly-even Gaussian `L²` to first Gaussian chaos.

Let `g` be a real even Gaussian `L²` function of squared norm one, with
normalized Hermite expansion

\[
g(z)=\sum_{r\ {m even}} b_r\frac{H_r(z)}{\sqrt{r!}},
\qquad \sum_r b_r^2=1.
\]

Put `D=sum_r r b_r²`, allowing `D=+infinity`. When `g` has a Gaussian weak
derivative in `L²`, this is exactly `E g'(Z)²`.

## 2. The forced coefficient recursion

Suppose a unit first-chaos Gaussian `V=sum_T u_T Z_T` satisfies

\[
V=\mathcal U[g(V)].                                  \tag{1}
\]

For a parent tree `T` with `r=sum_tau m_tau` children, the Hermite addition
formula gives

\[
\boxed{\quad
u_T=b_r\sqrt{\frac{r!}{\prod_\tau m_\tau!}}
             \prod_\tau u_\tau^{m_\tau}.
\quad}                                               \tag{2}
\]

Only the degree-`r` Hermite component of `g(V)` contributes to `h_T`.
The identity holds for infinite first-chaos series by `L²` convergence and
the usual finite-coordinate Hermite addition formula. Every child tree is
strictly smaller than its parent, so (2) uniquely forces every coefficient.
In particular, the single-edge coefficient is `u_edge=b_0`.

Define all coefficients recursively by (2), without assuming their squared
sum is one. Put `p_T=u_T²`. Then

\[
p_T=b_r^2\frac{r!}{\prod_\tau m_\tau!}
                \prod_\tau p_\tau^{m_\tau}.           \tag{3}
\]

This is exactly the probability recursion for an unordered Galton--Watson
family tree whose offspring law is `Pr(N=r)=b_r²`. An elementary generating-
function argument suffices, so no branching-process theorem is needed.

## 3. Total finite-tree mass

Let `s_k` be the total mass of trees with height at most `k`, taking
`s_0=0`. The multinomial formula, valid also for countable child types by
nonnegative monotone summation, gives

\[
s_{k+1}=F(s_k),\qquad F(s)=\sum_r b_r^2 s^r.          \tag{4}
\]

Thus `s_k` increases to the smallest fixed point `q∈[0,1]` of `F`. Every
finite tree has finite height, so

\[
\sum_T u_T^2=q.                                      \tag{5}
\]

If `D<1`, convexity or the derivative bound gives `F(s)>s` for `0≤s<1`,
and hence `q=1`. If `D=1`, the same strict inequality holds: the offspring
law has only even values and cannot be deterministic one-child. Unless it
is constant zero (which has `D=0`), `F` is strictly convex. Its tangent at
one is the line `s`, so `F(s)>s` for every `s<1`. Again `q=1`.

If `D>1`, including infinite mean, then `F(s)<s` for some `s<1` sufficiently
close to one, because the left derivative at one exceeds one. Since
`F(0)≥0`, the smallest fixed point satisfies `q<1`.

Consequently

\[
\boxed{\quad
\sum_Tu_T^2=1\ \Longleftrightarrow\ D\le1.
\quad}                                               \tag{6}
\]

For `D<1`, (4) also gives the quantitative height tail
`1-s_k≤D^k`. If `g` is a finite Hermite polynomial, its offspring law has
bounded support, and only finitely many trees have each fixed height. Thus
this is a direct finite-coordinate approximation bound for the numerical
certificate, albeit potentially with very many coordinates.

## 4. Exact existence and uniqueness criterion

If `D≤1`, the recursively defined coefficient vector has squared norm one.
Hence `V=sum_T u_T Z_T` is a well-defined unit Gaussian on the same
probability space. The Hermite addition formula used in (2) shows that
every coefficient of `U g(V)` equals `u_T`; therefore (1) holds in `L²`.

Conversely, every unit first-chaos solution of (1) must obey the uniquely
forced recursion (2). If `D>1`, its total squared coefficient mass would
equal `q<1`, a contradiction. This proves:

\[
\boxed{\quad
\text{There is a unit first-chaos fixed point }V=\mathcal U g(V)
\text{ iff }D\le1;\text{ it is then unique.}
\quad}                                               \tag{7}
\]

No positive-mean condition on `g` is needed. The sign of `b_0` and the signs
of the other Hermite coefficients propagate through (2), while squared
coefficient mass is governed by (3).

For `D<1`, this agrees with the independent Banach-contraction proof on the
complete unit sphere of first Gaussian chaos. The present argument also
handles the critical equality `D=1`, where the contraction rate disappears.

## 5. Relation to the original Boolean lower certificate

For any final mask `H=1_{|V|≤alpha}`, the first-chaos variable `W=UH` is
jointly Gaussian with `V`, with

\[
\operatorname{Var}W=2\Phi(\alpha)-1,
\qquad E VW=E[g(Z)1_{|Z|\le\alpha}].
\]

Finite-coordinate approximation of this bounded mask, followed by fixed
smooth approximation and the hierarchical tree-energy theorem, gives the
valid Boolean lower certificate `E|W|1_{|V|>alpha}`. The Gaussian fixed point
is not an infinite iteration run on the original signing. Every desired
accuracy is realized by a finite response construction before the dimension
limit is taken.

The strict-margin degree-200 candidate used in the exact numerical
certificate already has `D<1`; this critical extension is not needed to
validate its bound above `0.426`.
