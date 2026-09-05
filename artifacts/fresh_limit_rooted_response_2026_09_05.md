# Rooted-response decoupling for flat sign matrices

Date: 2026-09-05. Author: fresh-limit algebra subagent.

Status: proposed proof, undergoing independent adversarial audit. The finite
identities below are exact. The asymptotic polynomial-transfer argument is
spelled out so that a mistaken universality claim can be isolated. No Gaussian
limit for the nonlinear residual is assumed or needed.

## 1. Target and normalization

Let `m=n-1`, let `A` be hollow symmetric with off-diagonal entries `±1`, and set
`B=A/sqrt(m)`, `Q=B²`, and `G=BS`, where `S` has independent Rademacher entries.
Write `q(B)=max_x |xᵀBx|/2`. Consider any sequence with `q(A)≤C n^(3/2)` for a
fixed finite `C`. All little-oh statements below are uniform in the coordinate
`i` and in such sequences, with constants allowed to depend on a fixed smooth
function or polynomial degree.

The proposed conclusion is a strict universal improvement over the one-probe
constant `c*=0.336493364431...`, not convergence of the original minima.

For `τ>0` and fixed `t>0`, put `ψτ(z)=2Φ(z/τ)-1` and

\[
 f(g)={\psi_\tau(g+t)+\psi_\tau(g-t)\over2},\qquad
 h(g)={\psi_\tau(g+t)-\psi_\tau(g-t)\over2}.
\]

Then `f` is odd, `h` is even, and the conditional means of the two oriented
Boolean probes are

\[
 \mu_i^\sigma=\psi_\tau(\sigma G_i+tS_i)
                 =\sigma f(G_i)+S_i h(G_i).
\]

Define `u_i=f(G_i)`, `v_i=S_i h(G_i)`, and

\[
 a=\mathbb E f'(Z)=2\phi_{1+\tau^2}(t),\qquad
 b=\mathbb E h(Z)=2\Phi(t/\sqrt{1+\tau^2})-1,
\]

where `Z` is standard normal. Let `X^σ` be actual independent-coordinate
Gaussian-dither rounding conditional on `S`, with mean `μ^σ`. Convexity gives

\[
 \ell:={1\over2}\sum_\sigma\mathbb E\|BX^\sigma\|_1
 \ge\mathbb E\sum_i\max\{|(Bu)_i|,|(Bv)_i|\}.       \tag{1}
\]

The existing one-probe energy theorem gives

\[
 e:={\mathbb EH_B(X^+)-\mathbb EH_B(X^-)\over2}
       \ge ab\,n-o(n).                                \tag{2}
\]

This note supplies the proposed new lower bound for (1).

## 2. Elementary spectral bootstrap

Let `β(A)=max_{x,y∈{±1}^n}|xᵀAy|`. Multilinearity extends this maximum to the
two cubes. If `Az=λz`, use a random sign vector of mean `z/||z||∞`, and choose
the other signs to maximize the resulting linear functional. Then

\[
 \beta(A)\ge |\lambda|\,{\|z\|_1\over\|z\|_\infty}
                    \ge\lambda^2.
\]

The last inequality follows from the eigen-equation and `|aij|≤1`:
`|λ| ||z||∞≤||z||1`. Polarization on the cube gives `β(A)≤4q(A)`.
Consequently

\[
 \|B\|_{\rm op}^2=O(\sqrt n),\qquad
 \operatorname{Tr}Q^2\le\|Q\|_{\rm op}\operatorname{Tr}Q
                         =O(n^{3/2}),                 \tag{3}
\]

and `(Q²)ii≤||Q||op=O(sqrt(n))`. Here `Qii=1` and `TrQ=n` exactly.

## 3. Transport of every fixed even row function

Let `g` be fixed and even, smooth with derivatives through a sufficient fixed
order of at most polynomial growth. This includes a bounded smooth function
minus a fixed polynomial. Set `w_j=S_j g(G_j)`.

For `j≠k`, write `β=Bjk`,

\[
 G_j=U+\beta S_k,\qquad G_k=V+\beta S_j,
\]

where `(U,V)` is independent of `(Sj,Sk)`. Averaging the two endpoint spins
gives the exact identity

\[
 \mathbb Ew_jw_k={1\over4}\mathbb E
 [g(U+\beta)-g(U-\beta)]
 [g(V+\beta)-g(V-\beta)].                              \tag{4}
\]

Taylor expansion and two-dimensional Lindeberg replacement of the remaining
flat Rademacher sums therefore yield

\[
 \mathbb Ew_jw_k={1\over m}K_g(Q_{jk})+O_g(m^{-3/2}),
 \qquad K_g(q)=\mathbb Eg'(Z)g'(Z'),                    \tag{5}
\]

where the normal pair has unit variances and correlation `q`. The actual
common-field variances are `1-1/m`; replacing them by one costs `O_g(m^-1)`
inside the expectation, hence only `O_g(m^-2)` in (5). The derivative-product
test has polynomial growth, whose moments remain uniformly bounded because
all row sums are subgaussian. Thus the ordinary Lindeberg third-order
remainder is uniformly `O_g(m^-1/2)` before multiplication by `1/m`.

Because `g'` is odd, its Hermite kernel is

\[
 K_g(Q)=\sum_{r\ \mathrm{odd}}c_r Q^{\circ r},\qquad
 c_r\ge0,\quad\sum_rc_r=\mathbb Eg'(Z)^2.
\]

Schur multiplication by a positive semidefinite matrix of diagonal one is an
operator-norm contraction. Therefore

\[
 \|K_g(Q)\|_{\rm op}
       \le\mathbb Eg'(Z)^2\,\|Q\|_{\rm op}.            \tag{6}
\]

The diagonal of `Cov(w)` tends uniformly to `Eg(Z)²`. Sandwich (5) between
row `i` of `B` and its transpose. The leading off-diagonal contribution is
`O_g(||Q||op/m)=o(1)`. The entrywise remainder contributes at most
`O_g(m^-3/2)||B_i||1²=O_g(m^-1/2)`. Hence

\[
 \boxed{\quad\mathbb E(Bw)_i^2=\mathbb Eg(Z)^2+o(1).\quad} \tag{7}
\]

This is the crucial tail-transport estimate: a fixed polynomial approximation
to `h` may be transported by `B` without paying its growing operator norm.
Polarization also gives the corresponding covariance formula for two fixed
even functions.

## 4. Finite rooted Hermite chaoses and their local derivative

Let `r≥2` be fixed and even. The Rademacher rooted polynomial is

\[
 F_{i,r}=\sum_jB_{ij}S_jH_r(G_j),                       \tag{8}
\]

where `H_r` is the probabilists' Hermite polynomial. Up to an `o(1)` L² error,
this is the homogeneous multilinear degree-`r+1` polynomial obtained by
retaining only distinct indices in the `r` branch factors. Here is an explicit
finite recurrence proving that assertion, without an appeal to a polynomial
limit theorem. Set

\[
 W_{j,r}=r!\sum_{\substack{T\subset[n]\setminus\{j\}\\|T|=r}}
             \prod_{k\in T}B_{jk}S_k.
\]

Counting whether a multiplied index is already in `T` gives exactly

\[
 G_jW_{j,r}=W_{j,r+1}
       +r\left(1-{r-1\over m}\right)W_{j,r-1}.        \tag{8a}
\]

Compare (8a) with the Hermite recurrence. For every fixed `r`,
`H_r(G_j)=W_j,r+Σ_(s<r) c_(r,s)(m) W_j,s`, with matching parity and all
`c_(r,s)(m)=O_r(1/m)`. Formula (7), or induction in this triangular identity,
bounds the L² norm of every transported lower-degree term uniformly. For
example `H_2=W_2` and `H_4=W_4-(8/m)W_2-2/m` exactly.

For a fixed set of `r+1` distinct variables, its coefficient in the transported
homogeneous polynomial is the sum of at most `r+1` root choices. Its absolute
value is at most `C_r m^(-(r+1)/2)`. Thus each variable influence is `O_r(1/n)`
and the sum of influences is bounded. Joint Lindeberg replacement for this
finite family of multilinear polynomials and `G_i` has error `O_r(n^-1/2)`
against a fixed smooth test with bounded third derivatives: expand in one
input, use the matching first two moments, and bound the third moment of its
influence polynomial by fixed-degree hypercontractivity.

After replacing the inputs by independent standard normals `Z_j`, the same
multilinear polynomial differs in L² by `O_r(n^-1/2)` from

\[
 \mathcal F_{i,r}=\sum_jB_{ij}Z_jH_r((BZ)_j).           \tag{9}
\]

Indeed (9) is a Wick polynomial of degree `r+1`, since `Bjj=0`. In the product
Hermite basis, a multi-index `α` of total degree `r+1` can receive a root
contribution only from an index `j` with `αj=1`, so there are at most `r+1`
contributions. Each coefficient is bounded by `C_r m^(-(r+1)/2)`. The omitted
multi-indices have at most `r` distinct coordinates, numbering `O_r(n^r)`;
their Hermite factorial weights are bounded in terms of `r`. Their squared
norm is therefore `O_r(1/n)`.

For completeness, the exact Gaussian covariance of the rooted summands is

\[
 \mathbb E[Z_jH_r(G_j)Z_kH_r(G_k)]
 =r!\,\mathbf1_{j=k}
  +{r\,r!\over m}\,\mathbf1_{j\ne k}\,Q_{jk}^{r-1}.
                                                               \tag{10}
\]

It follows that `E Fcal_i,r²=r!+o(1)` uniformly. Distinct `r` are exactly
orthogonal because their total Gaussian chaos degrees differ.

Let `D_i` be differentiation in the unit Gaussian direction `B_i`. Directly,

\[
 D_i\mathcal F_{i,r}
 ={1\over m}\sum_{j\ne i}H_r(G_j)
   +r\sum_jB_{ij}Q_{ij}Z_jH_{r-1}(G_j).               \tag{11}
\]

The first term has squared L² norm at most
`r! Tr(Q²)/m²=o(1)`, because `r` is even and at least two. For the second,
the covariance matrix of `Zj H_(r-1)(Gj)` has operator norm bounded by a
constant depending only on `r`. At `r=2`, its extra matrix is a constant
multiple of `(J-I)/m`; at larger `r`, use the same Schur contraction as in
(6). The squared coefficient norm is

\[
 \sum_j(B_{ij}Q_{ij})^2
    \le{(Q^2)_{ii}\over m}=o(1).
\]

Therefore

\[
 \mathbb E|D_i\mathcal F_{i,r}|^2=o(1).                \tag{12}
\]

Decompose the Gaussian input into the scalar `G_i=B_i Z` and its orthogonal
complement. One-dimensional conditional Gaussian Poincaré shows that
`Fcal_i,r` differs in L² by `o(1)` from a random variable measurable in the
orthogonal complement. Thus every finite vector of rooted chaoses is
asymptotically independent of `G_i`. No assertion that these chaoses themselves
have Gaussian limits is used.

## 5. Decoupling of the full rooted response

Expand `h` in even Hermite polynomials,

\[
 h(z)=b+\sum_{r\ge2,\ r\ \mathrm{even}}
          {d_r\over r!}H_r(z),\qquad
 d_r=\mathbb E h(Z)H_r(Z).
\]

For a fixed truncation use Section 4; then use (7) on the smooth remainder
and let the truncation degree tend to infinity. Consequently every
subsequential joint law of

\[
 (G_i,R_i,F_{i,2}),\qquad R_i=(Bv)_i-bG_i,
\]

has `G_i` independent of `(R_i,F_i)`, with

\[
 \mathbb ER_i^2=\nu:=\mathbb Eh(Z)^2-b^2,
 \qquad \mathbb EF_i^2=2,
 \qquad \mathbb ER_iF_i=d_2.                           \tag{13}
\]

Here a uniformly random coordinate may be used before passage to a
subsequence; all prior errors were uniform. The L² polynomial approximation
also justifies moment passage. For the cubic witness,

\[
 \mathbb EF_i^4\le3^6(\mathbb EF_i^2)^2=2916.           \tag{14}
\]

For the present smoothed indicator,

\[
 d_2=-{t a\over1+\tau^2}.                              \tag{15}
\]

## 6. Removing the own spin from the competing field

For fixed `i`, let `G_j^(i)=G_j-Bji S_i`. Taylor expansion gives

\[
 (Bu)_i=W_i+{S_i\over m}\sum_{j\ne i}f'(G_j^{(i)})+o_{L^2}(1),
 \quad W_i=\sum_jB_{ij}f(G_j^{(i)}).
\]

The empirical average of the even smooth function `f'` tends to `a` in L².
Two-dimensional Lindeberg replacement costs `O(n^-1/2)`, and the Gaussian
covariance of an even centered function is bounded by a constant times
`Qjk²`; its averaged sum tends to zero by (3). Removing coordinate `i`
changes this by `o(1)`. Thus

\[
 (Bu)_i=aS_i+W_i+o_{L^2}(1).                           \tag{16}
\]

Similarly, replacing all `G_j` by `G_j^(i)` in `(Bv)_i` changes it by `o(1)`
in L². Its leading difference is

\[
 {S_i\over m}\sum_{j\ne i}S_jh'(G_j^{(i)}).
\]

The summands have bounded second moments and off-diagonal covariances
`O(1/m)` by endpoint-spin extraction, so this normalized sum has L² norm
`O(m^-1/2)`. The Taylor remainder is also `o(1)`. The resulting replacement
`V_i` and `W_i` are jointly independent of `S_i`.

For every real `w,v`, convexity and evenness imply

\[
 {\max(|a+w|,|v|)+\max(|-a+w|,|v|)\over2}
                  \ge\max(a,|v|).
\]

Combining this with (1), (16), and Section 5 yields

\[
 \boxed{\quad\liminf {\ell\over n}
       \ge\mathbb E\max(a,|bZ+R|),\qquad Z\perp R.\quad} \tag{17}
\]

## 7. Uniform strictness, without a residual CLT

Put `j(r)=E max(a,|bZ+r|)`. For `a,b>0`, this is even and strictly convex:

\[
 j''(r)={1\over b}\left[
 \phi((a-r)/b)+\phi((a+r)/b)\right]>0.
\]

Thus `j(r)>j(0)` for every nonzero `r`. Nonzero variance alone would not give
a uniform gap without a tail argument. The cubic witness supplies one
directly. For all sufficiently small fixed `τ`, one has `ν≤1`,
`|d2|≥0.4`, and (13)--(14). For `E={|R|>0.1}`, Cauchy--Schwarz and Hölder give

\[
 0.4\le|\mathbb ERF|
 \le0.1\sqrt2+
       (\mathbb ER^2)^{1/2}(\mathbb EF^4)^{1/4}
          \mathbb P(E)^{1/4}.
\]

Using the deliberately weaker `0.4-0.1sqrt2>0.2` and `2916<4000`,

\[
 \mathbb P\{|R|>0.1\}\ge{0.2^4\over4000}=4\cdot10^{-7}. \tag{18}
\]

Also, for sufficiently small `τ`, `a∈[0.5,0.6]`, `b∈[0.5,0.7]`.
On `|r|≤0.1`, the displayed formula for `j''` is greater than `0.4`.
Therefore `j(0.1)-j(0)>0.002`, and conservatively (17)--(18) imply

\[
 \liminf {\ell\over n}\ge j(0)+4\cdot10^{-10}.         \tag{19}
\]

As `τ↓0`, with `t=t*` the one-probe optimizer, `a=tb` in the limit and

\[
 j(0)\longrightarrow
 a\,\mathbb P(|Z|\le t)+b\,\mathbb E[|Z|1_{|Z|>t}]
 =2ab=2c_*.
\]

The energy baseline (2) also tends to `c*`. The exact random partial-response
inequality

\[
 (1+p^2)q(B)\ge(1-p)^2e+p(1-p)\ell                  \tag{20}
\]

then gives a strict gain by choosing a sufficiently small fixed `p>0` first,
and a sufficiently small fixed `τ>0` afterwards. The gain can be extremely
small; no numerical improvement is claimed until the full proof is audited.

## 8. Audit obligations and nonclaims

The principal obligations are the uniform polynomial-growth version of the
factored Lindeberg estimate (5), the explicit distinct-index comparison in
Section 4, and simultaneous moment passage in Section 5. Each has a specified
finite-degree mechanism, rather than a generic AMP universality assumption.

No claim is made about a Gaussian limit of `R`, arbitrary nonlinear response
fields, spectral flatness, the value `1/2`, or convergence of `M_n/n^(3/2)`.
The elementary spectral bootstrap only gives `||B||op=O(n^(1/4))`, which is
used in exactly the vanishing normalized contractions displayed above.
