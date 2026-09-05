# Independent audit: rooted Gaussian response lower bound

Date: 2026-09-05. Auditor: fresh-limit variational agent.

This is a self-contained condensation of the successful audit, separated
from the earlier variational exploration. The underlying proof was developed
by the parent and algebra agent in `fresh_limit_rooted_response_2026_09_05.md`
and `fresh_limit_rooted_gaussian_2026_09_05.md`. No convergence theorem for
the extremal sequence is asserted.

## Statement and normalization

Let `A` be symmetric, hollow, and off-diagonally ±1. Put `m=n−1`,
`B=A/sqrt(m)`, `Q=B²`, and `q(B)=max_x |x^TBx|/2` over sign vectors.
For a fixed low-cap sequence `q(A)≤C n^(3/2)`, the key theorem is:

> For every fixed smooth even `h` with polynomial-growth derivatives,
> `sum_j Bij Sj h((BS)_j)` converges to `N(0,E h(Z)²)`, uniformly over
> root coordinates `i`, where `S` has independent Rademacher entries and
> `Z` is standard normal.

Only this rooted marginal is claimed, not a full AMP limit or a Gaussian
limit for the nonrooted field `B f(BS)`.

The spectral bootstrap is valid. Every row of `A` belongs to the cube, so
`||A²||_(∞→∞)≤β(A)≤4q(A)`, where `β(A)=max_{x,y}x^TAy` over two cubes.
Consequently

\[
Q_{ii}=1,\quad \|Q\|_{op}=O(\sqrt n),\quad
\operatorname{Tr}Q^2=O(n^{3/2}),\quad
(Q^2)_{ii}\le\|Q\|_{op}.
\]

## 1. Even-function transport without an operator-norm loss

Set `G=BS`, `w_j=S_j g(G_j)` for fixed smooth even `g`. For `j≠k`, write
`Gj=U+Bjk Sk`, `Gk=V+Bjk Sj`. Averaging the two endpoint spins exactly gives

\[
Ew_jw_k=\tfrac14 E[g(U+B_{jk})-g(U-B_{jk})]
                       [g(V+B_{jk})-g(V-B_{jk})].
\]

Centered Taylor expansion followed by two-dimensional smooth Lindeberg
replacement yields `Ew_jw_k=m^(-1)K_g(Qjk)+O_g(m^(-3/2))`, where
`K_g(q)=E g'(Z)g'(Z')` for standard normal correlation `q`.
The common-field variances are `1−1/m`; adding independent Gaussian
increments changes this expectation by only `O_g(1/m)`.
This argument does not invert the covariance, so near-perfect correlations
cause no degeneracy. Polynomial-growth derivative moments are bounded by
row subgaussianity.

Because `g'` is odd, `K_g(Q)=sum_{r odd} c_r Q^{∘r}`, with nonnegative
coefficients summing to `E g'(Z)²`. Correlation-matrix Schur multiplication
contracts operator norm, hence `||K_g(Q)||op≤E g'(Z)² ||Q||op`.
Sandwiching by a flat unit row of `B` gives

\[
\boxed{E(Bw)_i^2=E g(Z)^2+o(1).}
\]

The kernel contribution is `O_g(||Q||op/m)` and the entrywise remainder
contributes `O_g(m^(-3/2))||B_i||_1²=O_g(m^(-1/2))`. The estimate applies
to a bounded smooth function minus a fixed polynomial. Limits must be
taken in the order dimension first, polynomial degree second.

## 2. Explicit polynomial transfer

For a row with coefficients `a_j=±m^(-1/2)`, let `P_r` be the ordered
distinct-index degree-`r` sum. Its exact recurrence is

\[
P_{r+1}=G P_r-r(1-(r-1)/m)P_{r-1}.
\]

Thus `H_r=P_r+sum_{s<r,s≡r mod2}O_r(1/m)P_s` at fixed degree; for example
`H4=P4−(8/m)P2−2/m`. The transport estimate in §1 bounds all lower terms
inductively. Root multiplication by `S_j` introduces no repetition because
`Bjj=0`.

In the transported homogeneous degree-`r+1` polynomial, any fixed support
has at most `r+1` root choices and coefficient at most
`C_r m^(-(r+1)/2)`. Each variable influence is therefore `O_r(1/n)`.
One-input Taylor replacement, with fixed-degree hypercontractivity for its
third remainder, transfers any fixed finite chaos list and the linear
field to Gaussian inputs with error `O_r(n^(-1/2))` on smooth tests.

For Gaussian inputs, the corresponding rooted Wick polynomial is

\[
F_{i,r}=\sum_j B_{ij}Z_jH_r(b_j\cdot Z),\qquad b_j=B_j.
\]

It is pure chaos of degree `r+1`, since `e_j·b_j=0`. Its repeated-index
tensor terms have at most `r` distinct indices and squared total norm
`O_r(1/n)`, so the multilinear comparison is valid. Endpoint Gaussian
integration by parts gives the exact covariance

\[
E[Z_jH_r(G_j)Z_kH_r(G_k)]
=r!1_{j=k}+\frac{r r!}{m}1_{j\ne k}Q_{jk}^{r-1}.
\]

For even `r≥2`, the variance of `F_i,r` tends uniformly to `r!`.

## 3. Complete contraction audit

Write `a_j=Bij`, `T=sum_j a_j e_j⊗b_j^{⊗r}`, and symmetrize `T`.
For contraction order `1≤ell≤r`, set `q=r−ell+1≥1` and `P=Q^{∘q}`.
Every contraction between permuted copies has exactly one of four root
placements, up to a norm-preserving permutation of uncontracted coordinates:

| Placement | Squared-norm bound |
|---|---|
| Roots paired to each other | `sum_jk a_j²a_k² Qjk^(2q)≤Tr Q²/m²` |
| Neither root contracted | `sum_jk a_j²a_k² Qjk^(2ell)≤Tr Q²/m²` |
| One root meets an opposite branch | `sum_k a_k² v_k^T P v_k≤||Q||op/m` |
| Both roots meet opposite branches | `Tr(CPC^TP)≤||Q||op²/m²` |

Here `(v_k)_j=a_j Bkj Qjk^(ell−1)`, so `||v_k||²≤1/m`. In the last row,
`Cjk=a_j a_k Bjk² Qjk^(ell−2)`, with zero diagonal and `||C||F²≤1/m²`.
The Schur contraction `||P||op≤||Q||op` is valid because `q≥1`.
All bounds vanish. The finite number of permutation terms depends only on
the fixed degree, so the triangle inequality closes the symmetrized case.

I read the primary statements directly:
[Nualart--Peccati (2005), Theorem 1](https://arxiv.org/pdf/math/0503598)
applies after variance normalization, and
[Nualart--Ortiz-Latorre (2007), Theorem 7](https://arxiv.org/pdf/math/0703240)
gives joint Gaussian convergence from componentwise convergence and limiting
identity covariance. Distinct chaos degrees, including the degree-one field,
are orthogonal. The finite-dimensional inputs embed in the common space
`ell²`. These hypotheses match exactly; no extra mixed-contraction claim
is imported.

Fixed-degree Rademacher transfer followed by §1's `L²` tail control proves
the stated rooted Gaussian theorem. Bounded second moments suffice to pass
expectations of linearly growing Lipschitz functions. Uniformity follows by
applying the argument to any chosen sequence of root coordinates.

## 4. From the marginal theorem to the Boolean lower bound

For fixed `t,τ>0`, set `ψτ(z)=2Φ(z/τ)−1`,
`f(g)=[ψτ(g+t)+ψτ(g−t)]/2`, and
`h(g)=[ψτ(g+t)−ψτ(g−t)]/2`. Then `f` is odd and `h` is even.
Put `u=f(G)`, `v=S h(G)`, `a=E f'(Z)`, `b=E h(Z)`.
Independent Gaussian-dither rounding with conditional means
`μσ=σu+v=ψτ(σG+tS)` gives actual Boolean probes `Xσ`.

The baseline can be checked without another theorem. With
`e=[E H_B(X+)−E H_B(X−)]/2`, conditional independence gives
`e=E u^TBv`. Averaging endpoint `Sj` in each term gives

\[
e=\frac1m\sum_{i\ne j}E[f'(Z_1)h(Z_2)]+O_\tau(\sqrt n),
\quad\operatorname{Corr}(Z_1,Z_2)=Q_{ij}.
\]

The two functions are even, so their centered covariance is bounded by
`Cτ Qij²`. Therefore `e=ab n+Oτ(sqrt(n))`.

Convexity bounds the mean best-response field by
`E sum_i max(|Bu_i|,|Bv_i|)`. Removing the own spin gives
`Bu_i=aS_i+W_i+o_(L²)(1)`. The derivative average concentrates because it
is even and `Tr Q²/n²→0`. The analogous own-spin dependence of `Bv_i`
has leading term `S_i m^(-1)sum_j Sj h'(G_j^(i))`; endpoint extraction
gives off-diagonal covariance `O(1/m)`, hence vanishing `L²` norm.
The replacement pair `(W_i,V_i)` is jointly independent of `S_i`.
Averaging that spin gives
`E max(|Bu_i|,|Bv_i|)≥E max(a,|Bv_i|)+o(1)`.
No tightness of `W_i` is required.

The rooted marginal theorem and then `τ↓0` yield

\[
L(t)=a(t)[2\Phi(a(t)/\sqrt{b(t)})-1]
 +2\sqrt{b(t)}\phi(a(t)/\sqrt{b(t)}),
\quad a(t)=2\phi(t),\ b(t)=2\Phi(t)-1.
\]

The exact random partial-update inequality
`(1+p²)q(B)≥(1−p)²e+p(1−p)ell` consequently gives

\[
\liminf_n\frac{M_n}{n^{3/2}}
\ge\frac{(1-p)^2a(t)b(t)+p(1-p)L(t)}{1+p^2}.
\]

The conversion from `q(B)/n` multiplies by `sqrt((n−1)/n)→1`.
Low-cap restriction loses no minimizing subsequence.

## 5. Independent exact-arithmetic certificate audit

Audited and ran `computations/fresh_limit_rooted_lower_certificate.py`.
It uses only `Fraction` and integer arithmetic, with outward rounding to
the grid `10^(-60)`. Addition, multiplication, reciprocal, and integer-square-
root enclosures are correct, including the extra upper grid unit for square
roots. Machin's identity gives π from alternating reciprocal-arctangent
series; the first omitted term bounds each remainder. The exponential
Taylor remainder is conservatively bounded by `3u^(d+1)/(d+1)!` on `[0,1]`.
The Gaussian-CDF series is the integrated alternating exponential series,
with remainder `x u^(d+1)/[(d+1)!(2d+3)]`; its terms decrease in the asserted
domain. Interval dependence only widens these bounds.

All argument-domain assertions hold for `t=7/8`, `p=8/125` and the derived
`a/sqrt(b)`. The exact lower interval endpoint exceeds `849/2500=0.3396`.
The script therefore certifies the decimal consequence of the displayed
formula without floating-point assumptions.

**Audit conclusion:** no fatal gap found in the structural proof or the
rational certificate. The result improves the lower bound to `0.3396`;
it does not settle existence of the original limit.
