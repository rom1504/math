# Two-field masks: rigorous extension of the rooted response formula

Date: 2026-09-05. Status: direct corollary of the audited rooted CLT and
two-site cap argument. Root is independently optimizing the resulting
Gaussian functional; no numerical optimum is asserted here.

## 1. Multivariate endpoint formula

Keep `B=A/sqrt(n-1)`, `G=BS`, `Y=B[S h(G)]`, and the low-cap hypothesis
`q(A)=O(n^(3/2))`. Let `h` be fixed, even, bounded, and smooth with bounded
derivatives. Let `F(g,y),H(g,y)` be fixed bounded smooth functions with bounded
derivatives through the orders used below. Then

\[
 \boxed{\begin{aligned}
 {1\over n}\mathbb E F(G,Y)^\top B[S H(G,Y)]
 \longrightarrow\;&
  \mathbb E F_g(G_*,Y_*)\;\mathbb EH(G_*,Y_*)\\
 &+\mathbb EF_y(G_*,Y_*)\;
       \mathbb E[h(G_*)H(G_*,Y_*)],
 \end{aligned}}                                      \tag{1}
\]

where `(G*,Y*)` is centered Gaussian with covariance matrix

\[
 \begin{pmatrix}1&b\\b&v\end{pmatrix},\qquad
 b=\mathbb Eh(Z),\quad v=\mathbb Eh(Z)^2.
\]

Here are all changes needed in the one-field proof in
`fresh_limit_custom_response_2026_09_05.md`.

1. In the `Sj` central difference, `D_jGi=Bij` exactly and `D_jGj=0`
   exactly. Thus the direct term gains
   `m^-1 Σ_(i≠j) E F_g(Gi,Yi) H(Gj,Yj)`.
2. The `F_y` term has the previous direct factor `Bij h(Gj)` and the
   previous indirect Jacobian `Tij`. Its endpoint Taylor errors remain
   `O(1/m)` per pair because `Bij²=1/m` and `E|D_jYi|²=O(1/m)`.
3. The indirect off-diagonal term is again
   `Σ_k E Sk h'(Gk) Ck`, now with bounded cube vectors formed from
   `F_y(Gi,Yi)` and `H(Gi,Yi)`. Its discrete derivative obeys

   \[
   |D_kC_k|\le C{\beta(B)\over m}
       \left(\max_i|D_kG_i|+\max_i|D_kY_i|\right).
   \]

   The first maximum is at most `m^-1/2`; the expectation of the second is
   `o(1)` by the already proved maximal-influence lemma. Hence the entire
   indirect term is `o(n)`.
4. The own-coordinate indirect term contains only `H_y`, since `D_jGj=0`.
   The identity `Tjj=Z0−Sj h'(Gj)/m` and the bilinear cap bound apply
   unchanged, giving `o(n)`.
5. The rooted finite-chaos theorem gives joint convergence of
   `(Gi,Yi,Gj,Yj)` to independent coordinate pairs whenever `(B²)ij→0`.
   Their cross covariance is `(B²)ij` times the displayed covariance matrix,
   up to a uniform `o(1)` correction. Only a vanishing fraction of ordered
   pairs have non-small `(B²)ij`, because `Tr(B⁴)=o(n²)`. Both direct sums
   therefore factor to the right side of (1).

These steps include all `Sj→Gk→Yi` dependencies and add no generic AMP or
higher-depth universality assumption.

If `|F|+|H|≤1`, the two orientation means
`μ_i^σ=σF(Gi,Yi)+Si H(Gi,Yi)` lie in the cube. Independent conditional
rounding realizes them as Boolean vectors. The half difference of their
expected oriented energies is exactly the left numerator in (1), so it is
at most `q(B)`.

## 2. A scalar functional of an arbitrary jointly even Gaussian mask

Let `G,Z` be independent standard normals. Let
`H:ℝ²→[0,1]` be measurable and jointly even:
`H(-g,-z)=H(g,z)` almost everywhere. Define

\[
 k(g)=\mathbb E_ZH(g,Z),\qquad
 p=\mathbb EH(G,Z),\qquad
 c^2=\operatorname{Var}(k(G)).
\]

If `H` is smooth and `c>0`, choose

\[
 h(g)={k(g)-p\over c}.
\]

Then `h` is bounded, smooth, and even, with `Eh=0`, `Eh²=1`. Therefore the
limiting fields `(G*,Y*)` in (1) are independent standard normals. Moreover

\[
 \mathbb E[h(G)H(G,Z)]=c.
\]

Gaussian integration by parts in (1) now gives, for every admissible smooth
odd `F`,

\[
 p\,\mathbb EF_g+c\,\mathbb EF_z
          =\mathbb E[(pG+cZ)F(G,Z)].                  \tag{2}
\]

At fixed `H`, the pointwise optimum under `|F|≤1-H` is

\[
 F(g,z)=\operatorname{sign}(pg+cz)(1-H(g,z)).
\]

Approximate the sign by `2Φ((pg+cz)/ε)-1` before taking the dimension limit,
then let `ε↓0`. This yields the universal certificate

\[
 \boxed{\quad
 \liminf_n{M_n\over n^{3/2}}\ge
 J(H):=\mathbb E\bigl[|pG+cZ|(1-H(G,Z))\bigr].
 \quad}                                              \tag{3}
\]

## 3. Measurable masks and all quantifiers

The certificate is not restricted to smooth masks. For `0<r<1`, let

\[
 H_r(g,z)=\mathbb E H(rg+\sqrt{1-r^2}G',
                        rz+\sqrt{1-r^2}Z'),
\]

with independent Gaussian `G',Z'`. This is a jointly even smooth function
valued in `[0,1]`, with bounded derivatives of every fixed order. Its mean is
exactly `p`; its conditional mean is the one-dimensional Gaussian smoothing
of `k`. As `r↑1`, `H_r→H` and `k_r→k` in Gaussian L², so `c_r→c`.
Consequently `J(H_r)→J(H)` by Cauchy--Schwarz and the finite second moments
of `G,Z`.

For each fixed `r<1` with `c_r>0`, construct the fixed smooth `h_r` above,
fix the final sign smoothing `ε>0`, and take `n→∞` first. Then let `ε↓0`
and finally `r↑1`. No smoothing parameter is allowed to shrink during a
dimension limit. If `c=0`, use any fixed centered normalized bounded smooth
even `h` instead; since `k` is constant, `E hH=0`, and the same argument gives
`J(H)=E|pG|(1-H)`. The smoothed masks retain constant conditional means in
this case.

Thus (3) holds for every fixed measurable jointly even mask. Taking a supremum
over masks is legitimate because each is an independent finite choice:

\[
 \liminf_n{M_n\over n^{3/2}}\ge
   \sup_{\substack{0\le H\le1\\H(-g,-z)=H(g,z)}}J(H). \tag{4}
\]

This is a two-dimensional Gaussian variational lower certificate, not an
assertion that it equals the original asymptotic optimum or proves
convergence.

## 4. Exact reduction from masks to one conditional profile

Fix an even measurable profile `k:ℝ→[0,1]`. Its mean `p` and variance `c²`
fix the linear form `U=pg+cz`. Among all masks having conditional mean
`E_Z H(g,Z)=k(g)`, maximizing `J` is equivalent, separately at each `g`, to
minimizing `E_Z |pg+cZ| H(g,Z)`. The elementary lower-tail rearrangement
therefore gives an optimizer

\[
 H(g,z)=\mathbf1_{\{|pg+cz|\le T(g)\}},
 \qquad \mathbb P_Z\{|pg+cZ|\le T(g)\}=k(g).          \tag{5}
\]

For `c>0`, the conditional distribution is continuous, so thresholds exist
without a boundary randomization (allow zero and infinity at profile values
zero and one). Since `k` is even, `T` is even and the mask is jointly even.
This shows that the supremum in (4) can be taken over one even scalar profile
`k`, with its conditionally optimal interval mask (5). It does not prove that
the profile optimization is convex or that a stationary solution is global.

For orientation, a simple upper bound on this whole certificate family is
available. Since `c²≤p(1-p)`, `VarU=p²+c²≤p`. At fixed total mask mass `p`,
removing the smallest `|U|` values minimizes removed absolute energy. Hence

\[
 J(H)\le2\sqrt p\,
      \phi\!\left(\Phi^{-1}((1+p)/2)\right).          \tag{6}
\]

The numerical maximum of this scalar upper bound is about `0.449555`.
This is only a ceiling for the two-field certificate, not an upper bound for
the original Boolean optimum, and its displayed decimal is not certified.
