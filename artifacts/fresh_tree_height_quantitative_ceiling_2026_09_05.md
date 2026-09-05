# Height residuals and a quantitative ceiling for nonnegative tree masks

Date: 2026-09-05. Status: independent proof with an exact scalar certificate.
This is a limitation of a lower-bound method, NOT an upper bound on the
original Boolean minimax constant.

## 1. Exact height masses and approximate fixed points

Let `U` be the countable Gaussianization isometry from the even child-tree
Hermite basis to the first-chaos tree coordinates. Let
`g=sum_r b_r H_r/sqrt(r!)` be even with `sum b_r²=1`, and put

\[
K(q)=\sum_r b_r^2q^r.
\]

For a unit first-chaos Gaussian `V`, let `P_L` project onto trees of height
at most `L`, where the single-edge tree has height one and `P_0=0`. Write
`m_L=||P_LV||₂²`. Child-basis completeness gives the exact identity

\[
\boxed{\qquad \|P_L U g(V)\|_2^2=K(m_{L-1}).\qquad}  \tag{1}
\]

Equivalently, it is the squared norm of the conditional expectation of
`g(V)` onto Gaussian coordinates of height at most `L-1`. The identity has
no coefficient-sign assumption and applies to every unit first-chaos `V`.

Let `delta=||V-Ug(V)||₂`. The one-projection reverse triangle inequality
gives the residual barrier

\[
\delta\ge\sup_{0<q<1}
          (\sqrt q-\sqrt{K(q)})_+.
\]

There is a stronger version. Fix `q` with `K(q)<q`, and take the first height
`L` for which `m_L≥q`; it exists because `m_L↑1`. Then
`m_(L-1)<q`, so the corresponding output mass in (1) is at most `K(q)`.
Projection and complementary projection are orthogonal. Applying the
reverse triangle inequality to both therefore gives

\[
\boxed{\quad
\delta^2\ge
 (\sqrt q-\sqrt{K(q)})^2
 +(\sqrt{1-q}-\sqrt{1-K(q)})^2.
\quad}                                               \tag{2}
\]

The right side increases if the first mass is raised beyond `q` or the
second is lowered below `K(q)`; equivalently it is twice one minus the
affinity of two Bernoulli distributions. Taking a supremum over `q` is
legitimate. This proves a uniform residual gap for every supercritical
scalar response, not just nonexistence of an exact fixed point.

## 2. Central-interval kernel and monotonicity

For `alpha>0`, let `p=2Phi(alpha)-1` and
`g_alpha=1_{|Z|≤alpha}/sqrt(p)`. Its kernel is

\[
K_\alpha(q)=
 \frac{\Pr\{|X|\le\alpha,|Y|\le\alpha\}}p,
\]

where `(X,Y)` is standard Gaussian with correlation `q`. At fixed
`0≤q<1`, this kernel is nondecreasing in `alpha`. Here is a direct proof.
Write `J(alpha)` for the numerator and

\[
k_{\rm edge}=
 \Pr\{|q\alpha+\sqrt{1-q^2}Z|\le\alpha\}.
\]

Differentiating both square boundaries gives `J'=4phi(alpha) k_edge`, while
`p'=2phi(alpha)`. Thus

\[
K_\alpha'=\frac{2\phi(\alpha)}p(2k_{\rm edge}-K_\alpha).
\]

The conditional interval probability is maximal at conditioning value zero.
The interval `[-(1+q)alpha,(1-q)alpha]` contains `[-alpha,0]`, so its Gaussian
mass is at least half the maximum central-interval mass. Consequently
`2 k_edge≥max_x Pr(|q x+sqrt(1-q²)Z|≤alpha)≥K_alpha`, proving monotonicity.

The exact certificate described below uses `q=19/20` and `alpha=27/40`.
It proves that the residual in (2) exceeds `37/200=0.185`. Since this alpha
has interval mass greater than one half, monotonicity implies the same
residual bound for every central interval of mass at most one half.

## 3. Quantitative rearrangement stability

Let `0≤H≤1`, put `p=EH`, `v=EH²`, `W=UH`, and, when `v>0`, define the
unit Gaussian `V=W/sqrt(v)`. Let `t` satisfy `Pr(|V|≤t)=p`, and put
`I=1_{|V|≤t}`. Define

\[
e=E|H-I|,\qquad
d=E|V|H-E|V|I.
\]

Equal mask masses and rearrangement give

\[
d=E\bigl[\bigl||V|-t\bigr|\,|H-I|\bigr]\ge0.
\]

The folded-normal density is at most `2phi(0)`. The mass at distance at most
`u` from the boundary is therefore at most `4phi(0)u`. Layer-cake integration
yields

\[
\boxed{\qquad d\ge e^2/(8\phi(0)).\qquad}             \tag{3}
\]

The isometry relates the scalar interval residual to this mask distance:

\[
\left\|V-U\frac I{\sqrt p}\right\|_2^2
=2-\frac{2(p-e/2)}{\sqrt{vp}}.
\]

If a uniform residual floor `delta_0≤sqrt(2)` applies, this identity forces

\[
e\ge 2p-(2-\delta_0^2)\sqrt{vp}\ge p\delta_0^2.       \tag{4}
\]

Let `R(p)=2sqrt(p)phi(t)` be the scalar rearrangement envelope. The actual
functional `J(H)=E|UH|(1-H)` satisfies

\[
R(p)-J(H)=2\phi(t)(\sqrt p-\sqrt v)+\sqrt v\,d.
\]

Combining (3)--(4), and using
`p² delta_0^4/(8phi(0))≤2phi(t)`, gives

\[
\boxed{\quad
R(p)-J(H)\ge\frac{p^{5/2}\delta_0^4}{8\phi(0)}.
\quad}                                               \tag{5}
\]

This accounts for fractional-mask variance loss as well as rearrangement
loss; it does not assume the mask is already binary.

## 4. Exact numerical consequence

The exact script `computations/fresh_tree_mask_ceiling_certificate.py` uses
only the previously audited rational interval primitives. To enclose the
noise kernel it uses the positive series

\[
pK_\alpha(q)=p^2+4\phi(\alpha)^2
 \sum_{r\ge2\ {m even}}
    q^r H_{r-1}(\alpha)^2/r!.
\]

Truncation at degree 400 has upper tail
`q^402 p(1-p)`, by the full indicator Hermite variance. This proves (2)
has squared residual greater than `(37/200)²` at the stated rational test
parameters. The kernel upper bound is approximately `0.8394371382614183`.

The scalar envelope has a unique maximum, characterized by
`t(2Phi(t)-1)=phi(t)`. Exact intervals put the maximizing `t` between
`0.65730655` and `0.65730656`, and give envelope upper
`0.449554965850946512`.

Outside the mass window `[12/25,1/2]`, monotonicity on either side of that
unique maximum gives envelope upper bounds `0.449466756157967013` on the
left and `0.449431142839039345` on the right. Inside the window, use
`delta_0=37/200` and (5), obtaining a uniform deficit of at least

\[
0.000058585474836468466908744104217991227763193449465078533174.
\]

Thus the completed exact calculation proves

\[
\boxed{\quad
\sup_{0\le H\le1,\ H\ {\rm jointly\ even}}J(H)
\le0.449496380376110043072253622560721563771043534601666379760435
<0.4495.
\quad}                                               \tag{6}
\]

The full output is saved in
`computations/results/fresh_tree_mask_ceiling_certificate.json`.

## 5. Scope

Equation (6) limits the nonnegative one-mask hierarchical Gaussian lower
certificate. It is NOT an upper bound on `M_n/n^(3/2)` and does not establish
convergence of those numbers.

The earlier crude rearrangement ceiling also applies to signed masks using
`|H|`, but this sharpened proof does not automatically do so: `|H|` near a
central interval need not imply `H` near a scalar interval indicator when
its sign varies. No sharper signed-mask ceiling is asserted here.

The director independently checked the first-crossing projection argument,
conditional-noise monotonicity, fractional-mask stability, and every scalar
interval in the script, then reran it with identical output. The algebra
reviewer also independently audited the theorem and certificate.
