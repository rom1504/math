# Hierarchical tree energy and a stable Gaussian-mask recursion

Date: 2026-09-05. The polynomial graph argument was proposed by the root and
independently audited by the algebra and variational agents. Smooth transport
uses the separately proved own-free Sobolev lemma. The final stable-recursion
corollary is proved below without claiming generic AMP state evolution.

## 1. Tree fields and Gaussian child polynomials

Let `𝒯` be the countable set of rooted-isomorphism classes of finite trees
whose distinguished external root has degree one and whose other vertices
all have odd degree. For `T∈𝒯`, let `d(T)` be the number of nonroot vertices,
and `a(T)=|Aut_root(T)|`. For a hollow symmetric signing `A`, define the
normalized injective polynomial field

\[
 X_{T,i}={1\over\sqrt{a(T)}m^{d(T)/2}}
  \sum_{\substack{f:V(T)\hookrightarrow[n]\\f(o)=i}}
    \prod_{uv\in E(T)}a_{f(u)f(v)}
    \prod_{v\ne o}S_{f(v)},\qquad m=n-1.              \tag{1}
\]

Every field is exactly independent of its own spin `Si`. For fixed finite
families and `q(A)=O(n^(3/2))`, the tree-chaos theorem gives joint convergence,
uniformly in `i`, to independent standard normals `Z_T` indexed by the trees.

Delete the external-root edge of `T`. Its first child has an even number of
remaining branches. Regarding that first child as the external root of each
branch gives child types `τ∈𝒯` with multiplicities `m_τ(T)`. Define

\[
 h_T(Z)=\prod_\tau
       {\operatorname{He}_{m_\tau(T)}(Z_\tau)
          \over\sqrt{m_\tau(T)!}}.                    \tag{2}
\]

For the single-edge tree `T1`, the empty product is `h_T1=1`. The root
automorphism identity is

\[
 a(T)=\prod_\tau m_\tau(T)!\,a(\tau)^{m_\tau(T)}.     \tag{3}
\]

All finite field families below are enlarged to contain the child types of
their members. Since children are strictly smaller, this closure is finite.

## 2. Hierarchical direct-energy theorem

For every fixed finite tree family and fixed polynomials `F,H` in its fields,

\[
 \boxed{\quad
 {1\over n}\mathbb E F(X_i)^\top B[S H(X_i)]
 \longrightarrow
 \sum_{T}\mathbb E\partial_TF(Z)\,
                \mathbb E[h_T(Z)H(Z)],
 \quad}                                              \tag{4}
\]

where `B=A/sqrt(m)`, the vector notation means evaluation at each root `i`,
and only finitely many derivatives in the sum are nonzero. The theorem
extends to fixed bounded smooth `F,H` with bounded derivatives.

### Polynomial proof: all powers and pairings

It suffices to take monomials `F,H`. In the expansion of
`Σij Bij E F(X_i) Sj H(X_j)`, let `D` be the total number of nonroot vertex
positions in all tree copies. Add the bridge edge `ij` and the explicit
spin `Sj`. After division by `n`, the normalization is asymptotic to
`n^(-(D+3)/2)`.

A nonzero spin expectation pairs or groups into larger even blocks the
`D+1` spin positions, including the explicit `Sj`. There are at most
`(D+1)/2` distinct spin labels. The external root `i` can contribute one more
label. Therefore the critical number of distinct labels is
`V=(D+3)/2`, and it is attained only when all spin positions are paired and
`i` occurs nowhere among those positions. Every lower-label pattern vanishes
by absolute counting. If `D+1` is odd the expectation is zero.

Reduce edge multiplicities in a leading quotient graph modulo two. If the
resulting graph is nonempty, pick any surviving edge. All its vertex labels
are summed here, including `i,j`; no Eulerian condition is needed. Fixing
the other labels leaves a bilinear cube sum bounded by `β(A)`. Thus its
normalized contribution is `O(β(A)/n²)=o(1)`.

If every edge multiplicity is even, the connected quotient graph has `V`
vertices and `D+1=2(V-1)` edge occurrences. Connectedness forces at least
`V-1` distinct edges; even multiplicity forces at most that many. Hence the
underlying graph is a tree and every edge occurs exactly twice.

The bridge's partner must be the top edge of one tree copy from `F`.
Indeed, its endpoint `i` occurs nowhere internally, so no other edge incident
to `i` can come from the interior of a tree or a tree rooted at `j`. The top
vertex of this selected `F` tree is `j`, paired with the explicit `Sj`.
Its child branches therefore must pair with tree copies from `H`: two child
branches cannot pair with each other, since they lie inside one injectively
labeled `F` tree. Whole-copy pairing propagation applies on all other
branches, just as in the tree-chaos theorem.

Selecting the `F` copy gives its formal derivative. The other `F` copies
produce Gaussian moments at root `i`. At `j`, the selected tree's child
branches must pair with `H` copies but not internally among themselves;
this is precisely the Wick/Hermite polynomial in (2). The automorphism
factors cancel according to (3), leaving one factor `1/sqrt(m_τ!)` per
child type. Thus the total surviving contribution is the right side of (4).
There are finitely many patterns at each fixed monomial degree, so every
remainder is uniform over low-cap signings and roots as required. □

### Smooth extension with ordered limits

The independent proof in
`fresh_ownfree_sobolev_transport_audit_2026_09_05.md` applies exactly to finite
collections of (1). Coefficient counting and hypercontractivity give
`||D_kX_i||p≤C_p/sqrt(n)` at every fixed finite `p`. For every fixed smooth
remainder `g` of polynomial growth, that note proves

\[
 \limsup_n\max_i\mathbb E|B[Sg(X)]_i|^2
   \le\|g\|_{L^2(\gamma)}^2
          +C_X\|\nabla g\|_{L^4(\gamma)}^2.           \tag{5}
\]

The underlying exact identity is
`E SjSk g(X_j)g(X_k)=E[D_k g(X_j) D_j g(X_k)]`, valid because the fields are
own-free. Cauchy--Schwarz then makes the normalized energy continuous in
its first argument in Gaussian `L²` and its second argument in Gaussian
`L²+W^(1,4)`. Polynomials are dense in the latter norm, with a self-contained
Fourier/Taylor proof in the cited note. Taking the dimension limit first for
each fixed polynomial approximation, and only afterwards taking the
approximation limit, proves (4) for the claimed smooth functions. No degree
or smoothing parameter grows during a dimension limit.

## 3. Gaussianization as an isometry

Let `Z=(Z_T)_(T∈𝒯)` be the countable independent Gaussian family. The
functions `(h_T)` are an orthonormal basis of the jointly even subspace of
`L²(Z)`. They are distinct normalized Hermite monomials by (2), and every
finite even Hermite monomial occurs: its multiset of child types defines
exactly one larger rooted tree.

Define the linear isometry `𝒰` on this even subspace by

\[
 \mathcal U H=\sum_{T\in\mathcal T}
                  \mathbb E[h_T(Z)H(Z)]\,Z_T.         \tag{6}
\]

The series converges in L², is a centered Gaussian linear form, and satisfies

\[
 \mathbb E[(\mathcal U H)(\mathcal U K)]=\mathbb E HK.
                                                               \tag{7}
\]

For every jointly even measurable `H(Z)∈[0,1]`, the direct-energy theorem
implies the universal certificate

\[
 \boxed{\quad
 \liminf_n{M_n\over n^{3/2}}
      \ge J_\infty(H):=\mathbb E|\mathcal U H|(1-H).
 \quad}                                              \tag{8}
\]

Here is the finite realization argument. First suppose `H` is smooth and
depends on finitely many tree coordinates. Choose a finite enlarged set `J`
containing those coordinates and their required children, and let
`V_J` be the orthogonal projection of `𝒰H` onto `span{Z_T:T∈J}`. Use the
bounded odd response
`F=ψ_ε(V_J)(1-H)`. Its orientation means `±F+S_iH` lie in the cube.
Gaussian integration by parts turns (4) into
`E V_J ψ_ε(V_J)(1-H)`. Let `J` increase and then `ε↓0` to get (8).

For general `H`, conditional expectations on finite Gaussian coordinate
sets preserve range and joint parity and converge in L². Gaussian smoothing
then gives bounded smooth finite-coordinate masks. The functional in (8)
is L²-continuous on `[0,1]`:

\[
 |J_\infty(H)-J_\infty(K)|\le2\|H-K\|_2,
\]

by the isometry (7) and Cauchy--Schwarz. Thus every infinite-coordinate mask
is approximated by legitimate finite constructions. This is a variational
closure, not an infinite algorithm run on the original matrix.

## 4. A stable nonlinear Gaussian recursion

Let `f:ℝ→[0,1]` be fixed, smooth, and even, with bounded derivative. Set

\[
 v=\mathbb Ef(Z)^2>0,\qquad D=\mathbb Ef'(Z)^2,
 \qquad D<v.                                        \tag{9}
\]

Define on the countable Gaussian space

\[
 H_0=\sqrt v,\quad V_t=\mathcal U H_t,\quad
 H_{t+1}=f(V_t/\sqrt v).                              \tag{10}
\]

Every `Ht` is jointly even and valued in `[0,1]`, and every `Vt` is centered
Gaussian of variance `v`. Write `q_t=E V_tV_(t-1)/v` for `t≥1`. Initially
`q1=Ef(Z)/sqrt(v)∈[0,1]`. The isometry gives the exact recursion

\[
 q_{t+1}=K(q_t),\qquad
 K(q)={\mathbb Ef(Z)f(Z')\over v},                    \tag{11}
\]

where the standard Gaussian pair has correlation `q`. The Hermite expansion
of `f` has nonnegative squared coefficients. Hence, on `[0,1]`,

\[
 K(1)=1,\qquad 0\le K'(q)\le D/v<1.
\]

Therefore `1-q_t` tends to zero geometrically. Moreover

\[
 \|H_t-f(V_t/\sqrt v)\|_2^2
     =2v[1-K(q_t)]\longrightarrow0.
\]

Substituting `Ht` into (8) and using Cauchy--Schwarz proves

\[
 \boxed{\quad
 \liminf_n{M_n\over n^{3/2}}
 \ge\sqrt{\mathbb Ef(Z)^2}\;
       \mathbb E[|Z|(1-f(Z))],
 \quad \mathbb Ef'(Z)^2<\mathbb Ef(Z)^2.
 \quad}                                              \tag{12}
\]

All operations in (10) are Gaussian function-space operations. For every
desired accuracy, choose a finite `t` first and approximate its bounded mask
by the finite constructions in Section 3. Then take the matrix dimension
limit. Only afterwards let the finite approximation accuracy improve and
`t` increase. Thus (12) does not assume a closed nonlinear AMP recursion on
the original matrix, and it does not assert convergence of its minima.

## 5. Unbounded intermediate responses and an exact Gaussian fixed point

The range restriction in Section 4 is unnecessary for the intermediate
Gaussian response. Only the final mask used in (8) must lie in `[0,1]`.
Let `g` be any even function in Gaussian `W^(1,2)` such that

\[
 \mathbb Eg(Z)^2=1,\qquad \mathbb Eg(Z)>0,\qquad
 D_g:=\mathbb Eg'(Z)^2<1.                              \tag{13}
\]

Set `V0=Z_T1=𝒰1` and `V_(t+1)=𝒰[g(Vt)]`. These are well-defined linear
Gaussian forms of variance one, even when `g` is unbounded. For
`q_t=E V_tV_(t-1)`, `q1=Eg∈(0,1]` and

\[
 q_{t+1}=K_g(q_t),\qquad
 K_g(q)=\mathbb Eg(Z)g(Z'),\qquad
 0\le K_g'(q)\le D_g<1\quad(0\le q\le1).             \tag{14}
\]

Thus `1-q_t≤D_g^(t-1)(1-q1)` and the adjacent distances
`||V_t-V_(t-1)||2=sqrt(2(1-q_t))` are summable. Their limit `V` belongs
to the closed Gaussian linear space, is standard normal, and satisfies

\[
                         V=\mathcal U[g(V)].          \tag{15}
\]

To pass to (15), use Hermite expansion: for any fixed Gaussian `L²`
function, its values at jointly standard Gaussian arguments converge in
`L²` when their correlation tends to one. No pointwise Lipschitz assumption
on `g` is required, and no intermediate response is rounded on the matrix.

For `α>0`, define `H=1{|V|≤α}`, `W=𝒰H`, and

\[
 p=2\Phi(\alpha)-1,\quad
 w=\mathbb E[g(Z)1_{\{|Z|\le\alpha\}}],\quad
 s^2=p-w^2.                                         \tag{16}
\]

By the isometry and (15), `(V,W)` is a centered jointly Gaussian pair with
`EV²=1`, `EW²=p`, and `EVW=w`. Therefore the universal certificate is

\[
 \begin{split}
 J_\infty(H)
 &=\mathbb E|wV+sZ|1_{\{|V|>\alpha\}}\\
 &=2w\phi(\alpha)
       \left[2\Phi\left({w\alpha\over s}\right)-1\right]
   +4\sqrt p\,\phi(0)
       \left[1-\Phi\left({\alpha\sqrt p\over s}\right)\right],
                                                               \tag{17}
 \end{split}
\]

where the last expression is interpreted continuously if `s=0`. Formula
(17) follows by integrating the conditional folded-normal expectation over
`V>α` and integrating its linear term by parts. It is even in `w`.

This proves a finite-checkable lower-bound family: any fixed even
polynomial `g` satisfying (13), together with any rational `α>0`, gives a
valid value of (17). The construction of the infinite Gaussian mask is
legitimate by the finite-mask approximation in Section 3. In particular,
polynomial growth of `g` creates no new matrix-transport requirement.

### Finite Hermite certificates

Write `g=Σ_(r even) u_r He_r/sqrt(r!)`, with finitely many nonzero
coefficients. Then the three hypotheses in (13) are exactly

\[
 \sum_r u_r^2=1,\qquad u_0>0,\qquad
 \sum_r r u_r^2<1.                                  \tag{18}
\]

For the central indicator, its Hermite coefficients are

\[
 c_0=p,\qquad
 c_{2j}=-{2\phi(\alpha)\operatorname{He}_{2j-1}(\alpha)
                     \over\sqrt{(2j)!}},\quad j\ge1,
 \qquad w=\sum_r u_r c_r.                            \tag{19}
\]

Thus all feasibility constraints can be checked rationally by taking
`g=(Σ v_r He_r)/sqrt(Σ r!v_r²)` with rational `v_r`. The remaining value
uses only rational arithmetic, square roots, `φ`, and `Φ`; exact outward
interval series supply a rigorous decimal lower certificate.

### Exact one-parameter optimization at a fixed threshold

For fixed `α>0`, maximizing the covariance `w` under the relaxed stability
constraint `Eg'^2≤1` is a convex Hilbert-space problem. Replacing
`Σu_r²=1` by `Σu_r²≤1` does not change its optimum: if this constraint had
slack, increasing `u0` would improve the objective without spending any
derivative energy. Its exact dual value is

\[
 w_*(\alpha)^2=
 \inf_{t\ge0}(1+t)\sum_{r\text{ even}}{c_r^2\over1+tr}.
                                                               \tag{20}
\]

Indeed, weighted Cauchy--Schwarz proves the upper bound for every `t`.
For `t>0`, put
`u_r(t)=C_t c_r/(1+tr)` with normalization `Σu_r(t)²=1`.
Its derivative energy is the weighted average of `r` with weights
`c_r²/(1+tr)²`. Differentiation expresses its derivative as minus twice the
covariance of the strictly increasing functions `r` and `r/(1+tr)`.
Consequently it is strictly decreasing, because `c0,c2` are nonzero.
It tends to zero as `t→∞`. For the central indicator it tends to infinity
as `t↓0`, since this nonconstant discontinuous function is not in Gaussian
`W^(1,2)`. Thus there is a unique `t>0` with derivative energy exactly one.
At that value weighted Cauchy--Schwarz is an equality, proving (20).

The boundary optimizer is approximated by responses satisfying strict
stability, by slightly increasing `t`; finite Hermite truncation and
renormalization can then preserve strictness. Finally, the value in (17)
is increasing in `w≥0`: its covariance derivative is

\[
 {d\over dw}\mathbb E|W|1_{\{|V|>\alpha\}}
   =2\phi(\alpha)
       \left[2\Phi\left({w\alpha\over\sqrt{p-w^2}}\right)-1\right]
   \ge0.                                                       \tag{21}
\]

One may obtain (21) either by differentiating (17), or by Gaussian
covariance differentiation followed by smoothing the two absolute-value
and indicator corners. Therefore (20) exactly optimizes the stable
fixed-point covariance certificate at each threshold; it is not a claim
of optimality among all Gaussian masks or all Boolean constructions.

## 6. Exact numerical certificate above 0.426

Take `α=37/50` and `a=97/10`. Let `c_r` be the central-indicator Hermite
coefficients in (19), and take the finite even polynomial

\[
 g(z)=C\sum_{\substack{0\le r\le200\\r\text{ even}}}
       {c_r\over a+r}{\operatorname{He}_r(z)\over\sqrt{r!}},
 \qquad
 C^{-2}=\sum_{\substack{0\le r\le200\\r\text{ even}}}
                         {c_r^2\over(a+r)^2}.        \tag{22}
\]

The exact-rational interval computation
`computations/fresh_limit_hierarchical_fixed_point_certificate.py`, with
stdout saved in the matching `computations/results/` JSON file, verifies

\[
 \mathbb Eg'^2\in
 [0.997762540769427182369956738328581413730480359774925265750549,
  0.997762540769427182369956738328581413730480359774925265751941]
 \subset(0,1),
\]

and `Eg>0.8493`. Its exact value of (17) is enclosed by

\[
 [0.426090354752424224090748824426603635740460994016689775818265,
  0.426090354752424224090748824426603635740460994016689775819420].
\]

Consequently the proved hierarchical certificate gives

\[
 \boxed{\displaystyle
 \liminf_{n\to\infty}{M_n\over n^{3/2}}>0.426.}
                                                               \tag{23}
\]

The script uses no floating-point arithmetic. Hermite recurrences and all
three weighted coefficient sums are exact rational sums; `φ(α)` and
`Φ(α)` use the existing outward Fraction interval primitive. The two larger
Gaussian CDF arguments, less than 2.5, use their exact rational endpoint
integrated Taylor sums through degree 256. Their alternating tails decrease
from that point onward for every argument at most 8, and the first omitted
term is an explicit error bound. Every remaining operation rounds outward
on a 60-decimal rational grid. This is a lower-bound result only; it does
not resolve convergence of the original normalized minima.

## 7. A ceiling for this energy-conversion mechanism

The entire certificate (8), not merely its scalar fixed-point subclass,
has a strict ceiling below `1/2`. If `p=EH` and `v=EH²`, then `v≤p` and
`𝒰H` is Gaussian of variance `v`. Among masks of fixed mass `p`, the
quantity `E|𝒰H|(1-H)` is maximized by filling the lowest absolute Gaussian
values. Therefore, with `α_p=Φ^(-1)((1+p)/2)`,

\[
 J_\infty(H)\le 2\sqrt v\,\phi(\alpha_p)
              \le 2\sqrt p\,\phi(\alpha_p).         \tag{24}
\]

The maximum of the last expression is approximately `0.44955496156`
(this decimal is numerical, not an interval certificate). Its unique
interior stationary point is characterized by
`α[2Φ(α)-1]=φ(α)`; the left side minus the right side has strictly positive
derivative `2Φ(α)-1+3αφ(α)` for `α>0`. This is only a ceiling on the single
orientation-averaged cross-energy conversion, not an upper bound on the
original optimization problem or on other possible Boolean constructions.
