# Second rooted layer: independent falsification and exact covariance audit

Date: 2026-09-05. Initial sections record the counterexample-first stage;
Section 4 records the subsequently completed independent CLT proof audit.
Exact identities below are proved; simulations are finite diagnostics only.

Let `A` be symmetric hollow with off-diagonal entries `±1`, `m=n-1`,
`B=A/sqrt(m)`, and `Q=B²`. For independent uniform signs `S`, define

\[
G=BS,\qquad Y=\frac1{\sqrt2}B[S(G^2-1)],\qquad
W=\frac1{\sqrt2}B[S(Y^2-1)].
\]

Products and squares inside brackets are coordinatewise. The question is
whether `(G_i,Y_i,W_i)` tends to independent standard normals along every
low-cap sequence `q(A)=O(n^(3/2))`. That hypothesis implies
`||Q||op=O(sqrt(n))`, `Qii=1`, and `Tr Q²=O(n^(3/2))`; it does not imply a
bounded operator norm for `B`.

## 1. Exact first-layer variance and second-layer linear covariance

For each unordered triple `T={r,s,t}`, put

\[
c_j(T)=B_{jr}B_{rs}B_{rt}
      +B_{js}B_{sr}B_{st}
      +B_{jt}B_{tr}B_{ts}.
\]

The exact multilinear expansion is

\[
Y_j=\sqrt2\sum_{|T|=3}c_j(T)S_T.
\]

In particular `Y_j` is pure third Walsh chaos. Squaring and grouping equal
and unequal root choices gives the identity

\[
\boxed{\quad EY_j^2=1+\frac{2(Q^2)_{jj}-3}{m}.\quad}       \tag{1}
\]

For distinct `j,l`, only pairs of triples whose symmetric difference is
`{j,l}` contribute to `E Sj Sl Yj²`. Their common part is an unordered pair
`{a,b}` disjoint from `j,l`. Since
`c_j({j,a,b})=2B_ab/m`, direct multiplication gives

\[
\boxed{\quad
E[S_jS_lY_j^2]
 =\frac4m B_{jl}(B^3)_{ll}
  +\frac{8(n-4)}{m^2}Q_{jl}.
\quad}                                                   \tag{2}
\]

For clarity, the term with the root `l` uses

\[
\sum_{a<b:\ a,b\notin\{j,l\}}B_{ab}B_{la}B_{lb}
 =\tfrac12[(B^3)_{ll}-2B_{lj}Q_{jl}],
\]

while the other two root terms combine to
`(n-3)Qjl/m`. These two contributions give the coefficient `n-4` in (2).

Substitute (1)--(2) into `E Gi Wi`. The result is

\[
\boxed{\begin{aligned}
\sqrt2\,E G_iW_i={}&
 \frac{2\operatorname{Tr}Q^2-2(Q^2)_{ii}-3m}{m^2}\\
&+\frac4m\sum_l B_{il}Q_{il}(B^3)_{ll}\\
&+\frac{8(n-4)}{m^2}[(Q^2)_{ii}-1].
\end{aligned}}                                          \tag{3}
\]

Every term tends to zero uniformly under the low-cap hypothesis. Indeed,
`(Q²)ii≤||Q||op`, and Cauchy--Schwarz gives

\[
\left|\sum_l B_{il}Q_{il}(B^3)_{ll}\right|
\le \sqrt{\frac{(Q^2)_{ii}}m}\,
     \sqrt{\operatorname{Tr}Q^2},
\]

because `|(B³)ll|²≤(Q²)ll`. Consequently

\[
\sup_i|E G_iW_i|=O(n^{-1/2}).
\]

There is therefore no surviving linear feedback coefficient in the `G_i`
direction. This does not eliminate possible surviving nonlinear terms or
establish the variance or Gaussianity of `W_i`.

I checked all displayed identities by exhaustive spin averaging for random
signings of every order 3 through 10. Errors were at floating-point roundoff
(`≤5e-15`). Algebra independently derived the same formulas and additionally
computed the full first-chaos matrix in
`computations/fresh_limit_second_rooted_polynomial_audit.py`.

## 2. Raw polynomial tests

The script `computations/fresh_second_rooted_layer_falsify.py` records
means, second/fourth/sixth moments, pair covariances, mixed fourth moments,
tail probabilities, and a bounded characteristic-function test. Each is
averaged both over coordinates and separately at root zero. Monte Carlo
standard errors use independent whole-spin samples, not an independence
assumption between coordinates.

The initial run has 32 models at orders 128, 256, and 512, with 8192 samples
per model. The data are in
`computations/results/fresh_second_rooted_layer_falsify_2026_09_05.json`.
Random signs give bulk `EW²≈1.179,1.076,1.039` as the order doubles, and
Hadamard/bounded-seed tensor examples also trend toward variance one and
fourth moment three. All observed correlations with `G` and `Y` trend to
zero.

Low-rank Gram-sign matrices and critical planted twins show much slower
finite convergence. For Gram rank ratio 0.1, bulk `EW²≈4.29,2.03,1.38`.
For one planted twin group of size about `sqrt(n)`, bulk `EW²≈2.48,1.88,1.36`,
while the root inside that group has `EW²≈7.06,5.18,3.29`. Its fourth moments
are extremely large and noisy. These defects decrease with dimension, so
they are not asymptotic counterexamples. They do caution against inferring
uniform-root moment control from coordinate-average numerical evidence.

The planted-twin models remain within the intended low-cap regime: all
unplanted edges can be chosen from a low-cap random-sign family, and replacing
the edges incident to `O(sqrt(n))` vertices changes the quadratic cap by only
`O(n^(3/2))`. Inside the planted group, all rows to its complement are
identical random sign patterns. The model is a deliberate stress test of
the weakest allowed spectral behavior.

A follow-up raw-polynomial run at orders 1024 and 2048 is saved in
`computations/results/fresh_second_rooted_layer_polynomial_large_2026_09_05.json`.
At order 2048, the random-sign model gives bulk `EW²≈1.0075, EW⁴≈3.1374`,
and the Gram rank-ratio-0.1 model gives `1.0712,4.0878`. The planted-twin
model gives bulk `1.0650,15.9986`; its planted root gives `1.6909,397.91`.
The latter fourth moment remains highly variable at 8192 samples. The
polynomial tails are therefore a much harsher stress test than the bounded
smooth version, despite the downward variance and covariance trends.

## 3. Bounded smooth masks

To separate polynomial tails from a stable distributional defect, replace
`H2/sqrt2` at both layers by

\[
h_K(x)=\frac{K\tanh((x^2-1)/K)-\mu_K}{\sigma_K},
\quad \mu_K=E[K\tanh((Z^2-1)/K)],
\quad \sigma_K^2=\operatorname{Var}(K\tanh((Z^2-1)/K)).
\]

For fixed `K`, this is bounded smooth even, centered and normalized in
Gaussian measure, with bounded derivatives of every fixed order. The
normalizing integrals are numerical here, since the purpose is falsification,
not exact certification.

At `K=4`, the completed larger-size run has 14 models of orders 512, 1024,
and 2048, with 8192 samples per model. Data are in
`computations/results/fresh_second_rooted_layer_bounded_2026_09_05.json`.
At order 2048, random signs give bulk `EW²≈1.0048, EW⁴≈3.0397`; the Gram
rank-ratio-0.1 model gives `1.0295,3.2632`; the bounded-seed tensors give
values very near `1,3`. The planted-twin model gives bulk `1.0159,3.2580`
and planted-root `1.2052,7.0909`, down from planted-root `1.5780,26.6726`
at order 512. Correlations again decrease. No stable non-Gaussian term has
been identified, but uniform-root higher moments remain a demanding proof
obligation.

No derivative correction has been introduced: the only explicitly
identified linear feedback term vanishes by the exact calculation above.

## 4. Later proof audit: the raw second layer is Gaussian

The elementary injective-tree moment argument subsequently supplied a proof,
independently audited in full in
`artifacts/fresh_tree_energy_independent_audit_2026_09_05.md`. The exact raw
second-layer reduction is written in
`artifacts/fresh_limit_second_rooted_tree_2026_09_05.md`.

Its conclusions supersede the earlier open proof status in this note:
under `q(A)=O(n^(3/2))`, uniformly in the output root,

\[
(S_i,G_i,Y_i,W_i)\Longrightarrow(S,Z_0,Z_1,Z_2),
\]

with three independent standard normals, independent of the Rademacher
`S`. The finite planted-root tails in Sections 2--3 are compatible with the
slow convergence permitted by this theorem. They are not counterexamples.

The proof does not infer Gaussianity from the vanished linear projection.
It removes the own-spin cavity error directly, removes the lower Fourier
degrees using a restricted four-tree moment calculation and a root-
multiplication inequality, and identifies the remaining normalized
seven-vertex injective tree. Its rooted automorphism count is eight.

## 5. Alternative root-cross contraction inequality

Before the elementary tree proof was found, an operator argument also
closed the two root-cross cases in the tensor-contraction approach. It is
recorded because it may be useful beyond the odd-tree class. Let `P_j` be
symmetric rank-`d` tensors and consider the rooted tensor
`sum_j a_j e_j tensor P_j`. For a contraction of size `r`, put

\[
A_r=\sum_j a_j^2(P_j\otimes_{r-1}P_j),\qquad
V=\sum_j a_j^2\|P_j\|^2.
\]

Regard `A_r` as the positive operator obtained by flattening `P_j` with
`r-1` input slots and `d-r+1` output slots. For a one-root-cross pattern,
define the block operator `M(e_j tensor u)=a_j P_j u`. Then
`MM*=A_r`; the squared contraction norm is at most `||A_r||op V`.

For the two-root-cross pattern (`r≥2`), write
`R_[u;(j,k,t)]=a_j P_j[k,t,u]`, where `t` has length `r-2` and `u` length
`d-r+1`. The contraction is `R Swap(j,k) R*`. As the swap is an isometry,

\[
\|R\,\mathrm{Swap}\,R^*\|_{HS}^2
 \le\|RR^*\|_{op}\operatorname{Tr}(RR^*)
 =\|A_r\|_{op}V.
\]

Thus both root-cross patterns vanish whenever the corresponding root-root
contraction `A_r` vanishes in Hilbert--Schmidt norm and `V` is bounded. This
is an actual operator estimate, not an uncontrolled partial-trace bound.
The elementary moment proof made this alternative unnecessary for the
specific raw second layer.
