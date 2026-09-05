# A second rooted response and its positive oriented energy

Date: 2026-09-05. Status: analytic theorem independently audited by the
director and two reviewers. See `fresh_custom_response_independent_audit_2026_09_05.md`
and `fresh_custom_response_literature_audit_2026_09_05.md`.
The smooth-mask optimization and exact certificate are in
`fresh_limit_response_variational_2026_09_05.md`.

## 1. Exact asymptotic formula

Let `A` be a hollow symmetric sign matrix, `m=n-1`, `B=A/sqrt(m)`, and
`q(A)=max_x |xᵀAx|/2=O(n^(3/2))`. Let `S` be uniform independent signs,
`G=BS`, and, for a fixed even bounded smooth function `h`, set

\[
 v_i=S_i h(G_i),\qquad Y=Bv,\qquad
 b=\mathbb Eh(Z),\quad v=\mathbb Eh(Z)^2.
\]

Assume all derivatives used below are bounded. The intended examples are
Gaussian-smoothed indicators. For fixed bounded smooth `F` odd and `H` even,
the proved formula is

\[
 \boxed{\quad
 {1\over n}\mathbb E F(Y)^\top B[S H(Y)]
 \longrightarrow
 \mathbb EF'(Y_*)\;\mathbb E[h(G_*)H(Y_*)],\quad}       \tag{1}
\]

where `(G*,Y*)` is centered jointly Gaussian with

\[
 \operatorname{Var}G_*=1,\qquad
 \operatorname{Var}Y_*=v,\qquad
 \operatorname{Cov}(G_*,Y_*)=b.                        \tag{2}
\]

The indirect `Sj→Gk→Yi` paths are included explicitly below. They are not
discarded merely because each individual derivative is small.

For the response means

\[
 \mu_i^\sigma=\psi_\epsilon(\alpha S_i+\sigma Y_i)
               =\sigma F(Y_i)+S_iH(Y_i),
\]

the half difference of the two oriented expected quadratic energies is
exactly the numerator on the left of (1). Independent conditional rounding
realizes these means as actual Boolean vectors, so the expression is at most
`q(B)`. With nonnegative `h,H` and increasing `F`, the limit in (1) is
positive.

## 2. Tools already proved

The rooted Gaussian theorem from `fresh_limit_rooted_gaussian_2026_09_05.md`
gives the one-coordinate limit `(Gi,Yi)⇒(G*,Y*)`. Its finite-chaos proof also
gives joint convergence at two coordinates whenever `(B²)ij→0`.

Write `Q=B²` and `β(B)=max_(x,y signs)|xᵀBy|`. The low-cap hypothesis gives

\[
 \beta(B)=O(n),\quad \|Q\|_{\rm op}=O(\sqrt n),\quad
 \operatorname{Tr}Q^2=O(n^{3/2}),\quad Q_{ii}=1.        \tag{3}
\]

The even-function transport lemma gives, for each fixed smooth even `g`,

\[
 \sup_i\mathbb E\left|\sum_k B_{ik}S_kg(G_k)\right|^2=O_g(1).
                                                               \tag{4}
\]

We also repeatedly use endpoint-spin extraction: for any fixed smooth `g`,
not necessarily even, `E Sk g(Gk)=0` and, for `k≠l`,

\[
 |\mathbb E[S_kg(G_k)S_lg(G_l)]|\le C_g/m.             \tag{5}
\]

Formula (5) follows exactly by averaging the two endpoint signs, producing a
product of two increments of size `O(m^-1/2)`.

## 3. Uniformly small maximal single-spin influence

Define the indirect Jacobian entries

\[
 T_{ij}(S)=\sum_k B_{ik}B_{kj}S_kh'(G_k).               \tag{6}
\]

Their means are zero. Since their coefficient vectors have squared norm at
most `1/m` and absolute sum at most one, (5) gives

\[
 \sup_{i,j}\mathbb ET_{ij}^2=O_h(1/m).                 \tag{7}
\]

For maximum control, a weaker high-moment estimate suffices. Fix `i,j` and
write `a_k=Bik Bkj`, so `||a||2²≤1/m`, `||a||1≤1`. Let `δ_l` denote the
half difference under flipping `Sl`. Taylor expansion of `h'` gives the
pointwise gradient decomposition

\[
 \delta_l T_{ij}
 =-a_lS_lh'(G_l)
   -S_l\sum_k B_{lk}a_kS_kh''(G_k)+R_l,
 \qquad |R_l|\le C_h/m.
\]

Consequently, for every input `S`,

\[
 \sum_l|\delta_lT_{ij}|^2
       \le C_h{1+\|B\|_{\rm op}^2\over m}
       =O_h(n^{-1/2}).                                \tag{8}
\]

Here is an elementary moment consequence, avoiding any unproved Gaussian
concentration transfer. If a centered cube function `T` satisfies
`Σ_l|δ_lT|²≤w` pointwise, Efron--Stein gives `ET²≤w`. Apply the same inequality
to `T^k`. Since `|T(S^l)-T(S)|≤2sqrt(w)`,

\[
 \operatorname{Var}(T^k)
   \le C_k w\,\mathbb E(|T|+2\sqrt w)^{2k-2}.
\]

Also `|ET^k|²≤ET²·ET^(2k-2)`. Induction yields
`E|T|^(2k)≤C_k w^k`. With `k=4` in (8),

\[
 \sup_{i,j}\|T_{ij}\|_8=O_h(n^{-1/4}),\qquad
 \sup_j\mathbb E\max_i|T_{ij}|=O_h(n^{-1/8})=o(1).    \tag{9}
\]

Let `D_jY_i=(Y_i(Sj=+1)-Y_i(Sj=-1))/2`; it is independent of `Sj`.
Its exact expression is

\[
 D_jY_i=B_{ij}h(G_j)
    +\sum_{k\ne j}B_{ik}S_k
       {h(G_k^{(j)}+B_{kj})-h(G_k^{(j)}-B_{kj})\over2},
\]

where `Gk^(j)=Gk-Bkj Sj`. Taylor expansion gives

\[
 \sup_j\mathbb E\max_i|D_jY_i|=o(1).                 \tag{10}
\]

Indeed the second term differs from `Tij(S)` pointwise by `O(m^-1/2)`, which
already suffices for (10). More precisely, the restoration of `Sj` in `h'`
produces

\[
 {S_j\over m}\sum_{k\ne j}B_{ik}S_kh''(G_k^{(j)}).
\]

The sum has bounded L² norm by (4): replacing `Gk^(j)` with `Gk` costs at
most a constant pointwise by the Lipschitz bound. All other remainders are
`O_h(1/m)` pointwise. Thus the more accurate expansion is

\[
 D_jY_i=B_{ij}h(G_j)+T_{ij}(S)+R_{ij},
 \qquad \sup_{i,j}\|R_{ij}\|_2=O_h(1/m).              \tag{11}
\]

Equations (5) and the exact increment formula also show
`E|D_jYi|²=O_h(1/m)` uniformly. These increments are pointwise bounded by a
constant depending only on `h`, so `E|D_jYi|³=O_h(1/m)` as well.

## 4. The two-site expansion

The left numerator of (1) is

\[
 E_n=\sum_{i,j}B_{ij}\mathbb E[F(Y_i)S_jH(Y_j)].
\]

Average `Sj` first and Taylor-expand the resulting central difference in the
two arguments `(Yi,Yj)`. The previous second/third-moment estimates allow us
to evaluate derivatives at the actual fields and use (11), with error
`O_h,F,H(1/m)` per pair before multiplication by `Bij`. As
`Σij|Bij|=n sqrt(m)`, the total error is `O(sqrt(n))=o(n)`. We obtain

\[
 E_n=D_0+D_1+D_2+o(n),                                 \tag{12}
\]

where

\[
 \begin{aligned}
 D_0&={1\over m}\sum_{i\ne j}
    \mathbb E[F'(Y_i)H(Y_j)h(G_j)],\\
 D_1&=\sum_{i,j}B_{ij}\mathbb E[F'(Y_i)H(Y_j)T_{ij}],\\
 D_2&=\sum_{i,j}B_{ij}\mathbb E[F(Y_i)H'(Y_j)T_{jj}].
 \end{aligned}                                       \tag{13}
\]

### 4.1 All indirect off-diagonal paths vanish together

Reorder `D1` using (6):

\[
 D_1=\sum_k\mathbb E[S_kh'(G_k)C_k],\qquad
 C_k=\sum_{i,j}B_{ki}B_{ij}B_{jk}F'(Y_i)H(Y_j).
\]

Put `r_ki=sqrt(m)Bki∈{0,±1}`. Then `Ck=m^-1 p_kᵀBq_k`, where
`p_ki=r_ki F'(Yi)` and `q_ki=r_ki H(Yi)` are uniformly bounded vectors.
The bilinear cube cap and Lipschitz bounds give the exact discrete estimate

\[
 |D_kC_k|\le C_{F,H}{\beta(B)\over m}\max_i|D_kY_i|.
\]

Since `h'(Gk)` is independent of `Sk`, discrete integration by parts, (3), and
(10) yield

\[
 |D_1|\le C_h\sum_k\mathbb E|D_kC_k|=o(n).            \tag{14}
\]

This controls every nested indirect path inside `Yi,Yj`; none was discarded.

### 4.2 The own-coordinate indirect term vanishes

Flatness gives the exact identity

\[
 T_{jj}=Z_0-{1\over m}S_jh'(G_j),\qquad
 Z_0={1\over m}\sum_kS_kh'(G_k).
\]

By (5), `EZ0²=O(1/m)`. The bilinear cap bounds
`|Σ_j(BF(Y))j H'(Yj)|≤C β(B)=O(n)`. Therefore the `Z0` part of `D2` is
`O(sqrt(n))`; its removed diagonal part is at most `C β(B)/m=O(1)`. Thus

\[
 D_2=o(n).                                            \tag{15}
\]

## 5. Factoring the direct term

For a fixed rooted Hermite degree, the covariance between its outputs at
coordinates `i,j` is `r! Qij+o(1)` uniformly, by sandwiching the exact rooted
covariance matrix between rows of `B`. The finite-chaos multivariate Gaussian
theorem, followed by the existing L² tail transport, therefore gives:
whenever `Qij→0`, the triple `(Yi,Gj,Yj)` converges to `Yi*` independent of
`(G*,Y*)`, where `Yi*~N(0,v)`.

Only a vanishing fraction of ordered pairs violate `|Qij|≤η` for any fixed
`η>0`, since `TrQ²/n²→0`. All three functions in `D0` are bounded. Sequential
compactness of arbitrary coordinate-pair choices now gives

\[
 {D_0\over n}\longrightarrow
 \mathbb EF'(Y_*)\;\mathbb E[H(Y_*)h(G_*)].             \tag{16}
\]

Equations (12), (14)--(16) prove formula (1).

## 6. Hard-threshold limit and numerical consequence

Take `h` to be a Gaussian smoothing of `1{|g|≤t}`, and `F,H` the odd/even
channels of a smoothed threshold with own-spin bias `α>0`. The dimension
limit is taken while both smoothing parameters are fixed. Only then send
them to zero. Write `b=2Φ(t)-1`; the limiting law in (2) has `VarY=b`,
`Cov(G,Y)=b`. Formula (1) becomes

\[
 d(t,\alpha)=
 {2\phi(\alpha/\sqrt b)\over\sqrt b}
  \mathbb P\{|G|\le t,\ |Y|\le\alpha\}.                \tag{17}
\]

Hence

\[
 \liminf_n {M_n\over n^{3/2}}
            \ge\sup_{t,\alpha>0}d(t,\alpha).           \tag{18}
\]

This conclusion uses the second response alone and does not require the
original one-probe energy bound or a negative-energy penalty.

Uncertified numerical optimization gives
`t≈1.3018228135916594`, `α≈0.7607953158649846`, and
`d≈0.3638997630424741`. At `t=0.875`, `α=2φ(t)`, it gives
`d≈0.3244492208341546`; mixing that response with the first probe yields
approximately `0.35854439769`.

These optimized decimals are diagnostics, not interval certificates. The
separate smooth-mask optimization has an exact certificate above 0.38578.
All hard thresholds here are removed after the dimension limit; no exact
fixed unsmoothed-rule energy limit is asserted. No convergence or upper-bound
improvement is claimed.
