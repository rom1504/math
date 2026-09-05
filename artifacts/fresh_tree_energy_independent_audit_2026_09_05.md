# Independent audit of injective-tree response energies

Date: 2026-09-05. Status: the combinatorial polynomial identity and its
smooth-function extension pass independent audit. The formulation below
does not assert a general recursively defined AMP state evolution.

## 1. Tree family and normalization

Let `mathcal T` be a fixed finite family of rooted trees, each with external
root of degree one and all nonroot degrees odd. Include all child-branch
types recursively; equivalently, adjoin any missing child coordinates and
let the response functions ignore them. Write `d(T)` for the number of
nonroot vertices and `a(T)=|Aut_o(T)|`.

For `m=n-1`, define the normalized injective field

\[
U_{T,i}=\frac{m^{-d(T)/2}}{\sqrt{a(T)}}
 \sum_{\substack{f:V(T)\hookrightarrow[n]\\f(o)=i}}
 \prod_{uv\in E(T)} A_{f(u)f(v)}
 \prod_{v\ne o}S_{f(v)}.
\]

Under `beta(A)=o(n²)`, the elementary tree moment theorem gives joint
convergence at each root to independent standard Gaussian coordinates
`Z_T`, uniformly in the root, with convergence of all fixed moments. The
fields are exactly independent of their own input spin.

I independently checked the tree CLT proof in
`fresh_limit_second_rooted_tree_2026_09_05.md`: leading moment partitions
pair all nonroot positions; odd degrees make the reduced parity graph
Eulerian; a nonempty parity graph has an edge away from the fixed root and
is killed by `beta(A)/n²`; an empty parity graph has connected quotient
`E=V`, so it is a tree with each edge used twice. Whole-copy pairing then
propagates from the first root edge because each nonroot block contains
exactly the two already paired copies. This establishes the claimed rooted
automorphism variance and joint Gaussian moments.

Let the child of the external root in `T` have child-branch multiplicities
`m_tau(T)`. Their total is even. Define

\[
h_T(Z)=\prod_{\tau\in\mathcal T}
        \frac{H_{m_\tau(T)}(Z_\tau)}{\sqrt{m_\tau(T)!}},
\]

using probabilists' Hermite polynomials and `H_0=1`. Examples are `h=1`
for the single-edge tree, `h=H_2(Z_edge)/sqrt(2)` for the three-vertex
binary branch, and `h=H_2(Z_binary)/sqrt(2)` for the seven-vertex tree.

## 2. Polynomial energy identity

For any fixed polynomials `F,H` in the field coordinates, the identity is

\[
\boxed{\quad
\frac1n E F(U)^T B[S H(U)]
\longrightarrow
\sum_{T\in\mathcal T}
 E[\partial_TF(Z)]\,E[h_T(Z)H(Z)].
\quad}                                                   \tag{1}
\]

To check it, expand each polynomial into monomials of tree fields. A term
has `D` nonroot spin positions from those fields, one additional spin `Sj`,
and the bridge edge `Bij`. Its total normalization is
`n^(-1) m^(-(D+1)/2)`.

Rademacher averaging requires the `D+1` spin positions to occur in even
blocks. At leading order every block has size two, and the unspun root `i`
is a separate label. This gives `(D+3)/2` free labels, including both roots.
Any coincidence between `i` and a spin block loses at least one free label
and vanishes by absolute counting. The label `j` is paired with exactly
one nonroot position in an `F`-tree; own-freeness excludes the `H`-trees.

Reduce edge multiplicities modulo two. Any nonempty reduced graph can be
bounded by selecting one edge, fixing all remaining labels, and using the
bilinear cut cap for its two free endpoints. This costs `beta(A)/n²=o(1)`.
Here both root labels are summed, so there is no fixed-root exception.
Distinctness exclusions are implemented by zeroing the two bounded test
vectors, with their mutual coincidence removed automatically by `A_rr=0`.

If the parity graph is empty, the quotient is connected, has `(D+3)/2`
vertices, and has `D+1` original edge occurrences. Every distinct edge occurs
at least twice. Connectedness forces the quotient to be a tree and every
edge to occur exactly twice.

At `i`, the bridge must pair with the top edge of exactly one `F`-tree `T`;
all other `F` copies pair by rooted isomorphisms. The first child of `T` is
the position paired with the explicit `Sj`. Every child branch of that
vertex must then pair with one `H`-tree rooted at `j`; injectivity prevents
pairing two branches within `T`. All remaining `H` copies pair among
themselves. Pairing propagates through the quotient just as in the tree CLT.

The choice of the distinguished `F` copy gives `E partial_T F`. At `j`,
the forced branch matchings and remaining Gaussian pairings give the
Hermite coefficient of `H` stated in (1). The normalization is important:

\[
a(T)=\prod_\tau m_\tau(T)!\,a(\tau)^{m_\tau(T)}.
\]

The child-isomorphism factors cancel the powers of `a(tau)` from the
normalized child coordinates, leaving precisely
`prod_tau 1/sqrt(m_tau!)`. This verifies the coefficient in `h_T`.

Only the odd part of `F` and the even part of `H` can survive this empty-
parity pattern. Formula (1) is nevertheless valid for general polynomials,
since its right side automatically discards the other parity components.

## 3. Smooth response functions and Boolean rounding

The exact own-free transport and polynomial-density proof in
`fresh_ownfree_sobolev_transport_audit_2026_09_05.md` extends (1) to fixed
bounded smooth response functions with bounded derivatives. The first
argument is continuous in Gaussian `L²`; the second is continuous in
Gaussian `L²` plus the gradient `L^4` norm. The dimension limit is always
taken with both polynomial degree and all smoothing parameters fixed.

If `|F|+|H|≤1`, the conditional means
`mu_i^sigma=sigma F(U_i)+Si H(U_i)` lie in `[-1,1]`. Independent conditional
rounding turns these means into actual Boolean vectors. The half difference
of their expected quadratic energies is exactly the numerator on the left
of (1), since `B` is hollow. Thus (1) is a valid lower certificate for the
quadratic cap.

The standalone raw second-layer statement also passes independent audit:
the exact cavity decomposition removes own-spin terms at `o(L²)` cost,
the restricted four-copy moment calculation shows that only the sixth
Fourier degree of the squared cavity field survives, and the root-
multiplication inequality transports its vanishing remainder. The leading
seventh-degree polynomial is the normalized injective seven-vertex tree,
whose rooted automorphism count is `2·2·2=8`. This proves that particular
module; the general energy identity above avoids needing a recursive
approximation theorem for every nonlinear response.

## 4. A useful variational formulation and its scope

For a jointly even mask `0≤H≤1`, put
`c_T(H)=E h_T(Z)H(Z)`. Gaussian integration by parts rewrites (1) as

\[
E\left[F(Z)\sum_T c_T(H)Z_T\right].
\]

At fixed smooth `H`, the optimal bounded odd channel is the smoothed limit
of `sign(sum_T c_T Z_T)(1-H)`. Therefore the family gives the certificate

\[
J_{\mathcal T}(H)
 =E\left[\left|\sum_T c_T(H)Z_T\right|(1-H(Z))\right].
\]

The child-multiplicity multiindices are distinct for distinct rooted trees,
so the `h_T` form an orthonormal collection in Gaussian space. In particular
`sum_T c_T(H)^2≤E H²`. This is a property of this response certificate,
not a characterization of the original Boolean optimum or a convergence
theorem for `M_n/n^(3/2)`.

### An upper ceiling for this response family

Put `p=EH` and `v=sum_T c_T²≤EH²≤p`. The linear form
`X=sum_T c_T Z_T` is Gaussian with variance `v`, regardless of its
dependence on the mask. Among all `[0,1]` masks of mean `p`, rearrangement
maximizes `E|X|(1-H)` by assigning `H=1` to the smallest absolute values
of `X`. Thus, with `t_p=Phi^{-1}((1+p)/2)`,

\[
J_{\mathcal T}(H)\le 2\sqrt v\,\phi(t_p)
                 \le 2\sqrt p\,\phi(t_p).
\]

The scalar upper envelope has maximum approximately `0.4495549615640635`.
Its maximizing threshold solves
`t(2Phi(t)-1)=phi(t)`, giving `t≈0.65730655`, `p≈0.48901618`. These decimals
are diagnostic evaluations of an exact analytic upper formula, not an
interval-certified numerical theorem. The same formula bounds signed
masks if `p` is replaced by `E|H|`, because then `v≤EH²≤E|H|` and the odd
channel budget is `1-|H|`.

This ceiling applies only to the stated one-shot tree-response certificate.
It is not an upper bound on the original Boolean minimax constant.

### Countably infinite indexing observation

Over all finite admissible rooted trees, the functions `h_T` form the
complete jointly-even Hermite basis on the countable Gaussian tree
coordinates. Indeed every finite multiindex with even total degree is a
finite multiset of child trees and hence defines a unique parent tree;
the inverse operation recovers that multiset. The output tree is strictly
larger than each child tree.

Consequently the formal map
`H -> sum_T <H,h_T> Z_T` extends as an isometry from the jointly-even
Gaussian `L²` space to first Gaussian chaos. Its coefficients have a strict
grading constraint through the child-tree structure. This is a useful
organization of possible infinite-depth variational limits, not a proof
that the scalar envelope above is attained or that finite tree certificates
characterize the original optimization problem.

## 5. Independent audit of the stable Gaussian fixed point

Let `U` denote the isometry from the complete even Hermite basis to first
Gaussian chaos described above. Let `g` be even, belong to Gaussian
`W^{1,2}`, and satisfy

\[
E g(Z)^2=1,\qquad D:=E g'(Z)^2<1.
\]

Intermediate `g` need not be bounded or smooth: it will never be evaluated
on the original matrix. This slightly strengthens the initially proposed
bounded-smooth formulation.

Consider the unit sphere `S` of first Gaussian chaos, with its `L²` metric.
It is a complete metric space, being a closed subset of a Hilbert space.
Each `V∈S` is a standard Gaussian linear form in the same countable
underlying coordinate family. Hence

\[
\mathcal R(V)=U[g(V)]
\]

belongs to `S`. The input is jointly even and has norm one. If `V,W∈S`
have correlation `q`, write
`K(q)=E g(V)g(W)`. The Hermite expansion gives an absolutely and uniformly
convergent derivative series on `[-1,1]`, with `|K'(q)|≤D`. Since `K(1)=1`,

\[
\begin{aligned}
\|\mathcal R(V)-\mathcal R(W)\|_2^2
 &=2[1-K(q)]\\
 &\le 2D(1-q)=D\|V-W\|_2^2.
\end{aligned}
\]

Banach's contraction theorem therefore gives a unique fixed point `V∈S`
satisfying

\[
\boxed{\qquad V=U[g(V)]\quad\hbox{in }L².\qquad}
\]

This is an equality on one common probability space, not merely an identity
of marginal distributions. All iterates are jointly Gaussian because all
remain in first chaos. No separate positive-mean assumption on `g` is
needed for existence, although the optimizing examples have positive mean.

### Boolean certificate from the fixed point

Fix `alpha>0`, set `H=1_{|V|≤alpha}` and `W=UH`. Then `(V,W)` is jointly
Gaussian, and

\[
\operatorname{Var}V=1,\quad
\operatorname{Var}W=p:=2\Phi(\alpha)-1,\quad
\operatorname{Cov}(V,W)
 =E[g(Z)1_{|Z|\le\alpha}]=:w.
\]

The covariance follows from the actual fixed point and the isometry:
`E VW=<Ug(V),UH>=E g(V)H`. The countable-mask certificate is valid by
finite-coordinate conditional approximation, fixed Gaussian smoothing, and
`L²` continuity; it does not run an infinite iteration on a matrix. It gives

\[
\liminf_n M_n/n^{3/2}\ge E|W|1_{|V|>\alpha}.
\]

For `w≥0`, write `sigma=sqrt(p-w²)`. Direct integration yields

\[
\boxed{\quad
2w\phi(\alpha)
 [2\Phi(w\alpha/\sigma)-1]
 +2\sqrt{2p/\pi}\,
   \overline\Phi(\alpha\sqrt p/\sigma).
\quad}                                                \tag{2}
\]

The formula extends continuously to `sigma=0`; it is even in `w`. At fixed
`p`, its derivative with respect to `w≥0` is exactly
`2phi(alpha)[2Phi(w alpha/sigma)-1]≥0`. Thus maximizing the threshold
correlation `w` under the norm and derivative constraints is the correct
scalar optimization.

### Finite-polynomial certificate target

An independently computed simple parameter choice is `alpha=37/50` and
`a=97/10`. Let `beta_r` be the normalized Hermite coefficient of
`1_{|Z|≤alpha}`, retain even degrees `0≤r≤200`, and take

\[
g_r=\frac{\beta_r/(a+r)}
          {\sqrt{\sum_{s\le200}\beta_s^2/(a+s)^2}}.
\]

This is a finite polynomial, so the Gaussian `W^{1,2}` requirement is
automatic. Independent floating-point calculations give
`D≈0.99776254077<1` and (2) approximately `0.4260903547524`.
The mathematical fixed-point argument applies directly to this polynomial
and needs no bounded intermediate-mask approximation.

### Completed independent exact-arithmetic audit

I read every line of
`computations/fresh_limit_hierarchical_fixed_point_certificate.py`, reread
its imported exact interval primitives, and independently reran the script.
All checks pass. The normalized derivative energy lies in

\[
\begin{split}
0.997762540769427182369956738328581413730480359774925265750549
\le D\\
\le
0.997762540769427182369956738328581413730480359774925265751941<1.
\end{split}
\]

The resulting lower certificate lies in

\[
\begin{split}
0.426090354752424224090748824426603635740460994016689775818265
\le L\\
\le
0.426090354752424224090748824426603635740460994016689775819420.
\end{split}
\]

In particular, `L>213/500=0.426` exactly. The calculation uses no floating
point and no infinite Hermite tail: its intermediate function is the
specified finite polynomial of degree 200.

The new CDF routine evaluates an exact rational integrated Taylor sum
through degree 256 at the two interval endpoints. Its terms need not
decrease at the start, but the first omitted term and every later term do
decrease for the asserted domain `0≤x≤8`: the successive magnitude ratio
is at most `(x²/2)/(degree+2)<1`. Therefore the first omitted alternating
term rigorously bounds the remaining tail. Monotonicity of `Phi` then gives
the interval hull. The argument bounds in this run are approximately 2.378
and 2.490, inside the checked domain.

The finite sums correctly use
`beta_0²=p²` and `beta_r²=4phi(alpha)² H_(r-1)(alpha)²/r!` for even `r≥2`.
The separate sums for norm, derivative energy, and threshold covariance
have denominators `(a+r)²`, `(a+r)²` with numerator `r`, and `(a+r)`,
respectively. Normalizing by the exact norm produces the actual candidate
used in the proved Gaussian fixed-point construction, not a relaxed
surrogate.
