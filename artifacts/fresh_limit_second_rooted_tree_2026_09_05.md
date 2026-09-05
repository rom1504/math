# The second cubic rooted layer: an elementary tree-chaos theorem

Date: 2026-09-05. The tree moment argument was proposed by the director and
independently checked by all three agents. The variational and literature
audits are saved separately. Exact finite covariance identities below were
also checked by exhaustive cube evaluation at forty signings; those tests
are diagnostics of the analytic identities.

## 1. Rooted odd-degree tree theorem

Let `T` be a fixed finite tree with a distinguished external root `o` of
degree one. Assume every other vertex of `T` has odd degree. Write
`d=|V(T)|-1`; necessarily `d` is odd. For a hollow symmetric sign matrix `A`
of order `n`, fixed coordinate `i`, and independent Rademachers `S`, define

\[
 X_{T,i}=m^{-d/2}
  \sum_{\substack{f:V(T)\hookrightarrow[n]\\f(o)=i}}
      \prod_{uv\in E(T)}a_{f(u)f(v)}
      \prod_{v\ne o}S_{f(v)},\qquad m=n-1.            \tag{1}
\]

The map is injective on all vertices, including the external root. Put
`β(A)=max_(x,y signs)|xᵀAy|`. If `β(A)=o(n²)`, then, uniformly in `i`,

\[
 X_{T,i}\Longrightarrow
 N(0,|\operatorname{Aut}_o(T)|).                       \tag{2}
\]

For any fixed finite family of such rooted trees, their joint moments tend
to Gaussian Wick moments; covariance is the number of root-preserving
isomorphisms between the two trees. In particular, nonisomorphic rooted trees
have independent Gaussian limits. The fixed sign `Si` is independent of all
polynomials in (1), exactly, since the external-root label is excluded from
the Rademacher factors.

### Proof by moment partitions

Expand a joint moment of `l` tree polynomials, with total `D` nonroot vertex
positions. A nonzero Rademacher expectation requires every free label to
occur an even number of times. Thus there are at most `D/2` free labels.
Patterns with fewer free labels contribute `o(1)` by absolute counting after
the normalization `m^(-D/2)`. If `D` is odd the expectation is exactly zero.

In a leading pattern, every free label identifies exactly two vertex
positions. Within-copy injectivity forbids identifying two vertices in the
same tree copy. Form the quotient multigraph with one fixed root and these
`V=D/2` free vertices. Reduce all edge multiplicities modulo two, obtaining
a simple graph `H`.

Every free quotient vertex comes from two odd-degree vertices, so its degree
in `H` is even. The root has degree `l`, which is even in every leading
pattern. Therefore `H` is Eulerian.

If `H` is nonempty, it contains an edge `uv` with both endpoints different
from the fixed root: a nonempty simple Eulerian graph contains a cycle, and
a cycle has at least three vertices. Fix all remaining free labels. The sum
over the two labels at `u,v` has the form

\[
 \sum_{r,s}a_{rs}f(r)g(s),\qquad |f|,|g|\le1,
\]

where zero values enforce exclusion of the other labels and the fixed root.
The exclusion `r≠s` is automatic since `arr=0`. Its absolute value is at most
`β(A)` by multilinearity on the two cubes. Summing the remaining labels and
normalizing bounds this entire pattern by `O_T(β(A)/n²)=o(1)`, uniformly in
the root coordinate. All even edge powers have canceled to one.

It remains to classify patterns with `H` empty. Every underlying quotient
edge then occurs at least twice. There are `D` original edge occurrences, so
the number `E` of distinct quotient edges is at most `D/2=V`. On the other
hand the underlying quotient graph is connected, since all tree copies share
the root, and it has `V+1` vertices. Therefore `E≥V`. Equality follows: the
underlying quotient is itself a tree, and every edge occurs exactly twice.

Each original copy has only one edge incident to the external root. Its
paired occurrence belongs to another copy. At the first nonroot vertex,
the pair-partition block contains precisely those two copies; no third copy
can meet that vertex. Since every incident edge occurs twice, the two copies
must match on every incident edge. Propagate this statement away from the
root through the quotient tree. The entire two copies are paired by a
root-preserving isomorphism. Thus all surviving patterns are exactly whole-
copy Gaussian pairings. For each pairing and each choice of rooted
isomorphisms, the number of injective free labelings is `(n-1)_V`, whose
ratio to `m^V` tends to one. This proves the stated Wick moments.

There are finitely many partition patterns at every fixed moment order.
Their bounds also provide uniform moment bounds and tightness. Gaussian
moment determinacy, or the same moment calculation for every finite linear
combination, proves joint convergence. This proves (2). □

The assumption that all nonroot degrees are odd is material: it guarantees
that a nonempty parity graph has a nonroot edge. No arbitrary-tree or full
AMP assertion is being made.

## 2. Exact first-layer covariance and second-layer linear correction

Now assume the stronger low-cap hypothesis `q(A)=O(n^(3/2))`, put
`B=A/sqrt(m)`, `Q=B²`, and define

\[
 G=BS,\qquad
 Y={1\over\sqrt2}B[S(G^2-1)],\qquad
 W={1\over\sqrt2}B[S(Y^2-1)].                          \tag{3}
\]

All products and squares inside brackets are coordinatewise. The bilinear
cap and spectral bootstrap give

\[
 \beta(A)\le4q(A)=O(n^{3/2}),\quad
 \|Q\|_{\rm op}=O(\sqrt n),\quad Q_{ii}=1.            \tag{4}
\]

The following covariance identities are exact at every order:

\[
 \operatorname{Cov}(Y)
        ={m-3\over m}Q+{2\over m}Q^2,                 \tag{5}
\]

and, for `j≠l`,

\[
 \mathbb E[S_jS_lY_j^2]
    ={4\over m}B_{jl}(B^3)_{ll}
       +{8(m-3)\over m^2}Q_{jl}.                     \tag{6}
\]

For example, the degree-three coefficient of `Yj` on a triple containing
`j` is `(2sqrt(2)/m) Bab`; multiplying the two cubic expansions and summing
their common two-set gives (6). Consequently the complete first-chaos
coefficient matrix `L_il=E Wi Sl` is

\[
 \boxed{\begin{aligned}
 L_{il}={1\over\sqrt2}\bigg\{&
 {B_{il}[2(Q^2)_{ll}-3]\over m}
 +{4Q_{il}(B^3)_{ll}\over m}\\
 &+{8(m-3)\over m^2}[(B^3)_{il}-B_{il}]\bigg\}.
 \end{aligned}}                                      \tag{7}
\]

Since `(Q²)ll≤||Q||op`, `|(B³)ll|≤sqrt((Q²)ll)`, and
`||B³_i||2≤||Q||op`, (4) gives

\[
 \sup_i\sum_l L_{il}^2=O(1/n).                        \tag{8}
\]

Thus there is no surviving linear Onsager correction in this explicit second
layer. This calculation alone would not establish a Gaussian limit; the
remaining argument removes all lower-degree pieces.

The script `computations/fresh_limit_second_rooted_polynomial_audit.py`
exhaustively verifies (5)--(7) on forty deterministic/seeded signings at
orders three through ten. It is a floating-point smoke test of exact
identities, not the proof of them.

## 3. A root-multiplication bound that avoids operator amplification

If `P_j` is a Boolean polynomial of degree at most `r` independent of `Sj`,
then for arbitrary real coefficients `a_j`,

\[
 \boxed{\quad
 \left\|\sum_j a_jS_jP_j\right\|_2^2
       \le(r+1)\sum_j a_j^2\|P_j\|_2^2.
 \quad}                                              \tag{9}
\]

Indeed, the output Fourier coefficient on a set `T` is
`Σ_(j∈T) a_j P̂_j(T\{j})`. Cauchy--Schwarz costs at most `|T|≤r+1`;
summing coefficients gives (9). In particular, a uniformly vanishing L²
remainder can be transported after multiplying by its own independent root
spin without paying `||B||op`.

## 4. Removing the own-spin dependence in the first layer

There is an exact decomposition

\[
 Y_j=P_j+S_jU_j,
 \qquad U_j={\sqrt2\over m}
       [S^\top BS-2S_jG_j],                           \tag{10}
\]

where both `Pj` and `Uj` are independent of `Sj`. The polynomial `Pj` is
exactly `X_(T3,j)/sqrt2`, where `T3` is the tree with an external root, one
child, and two leaves. Its rooted automorphism count is two. Formula (10)
follows by expanding `(Gk^(j)+Bkj Sj)²` in (3).

For every fixed finite `p`, `||Uj||p=O(n^-1/2)` uniformly: the quadratic
form `SᵀBS` has L² norm `sqrt(2n)`, and fixed-degree hypercontractivity controls
its higher moments. The tree theorem bounds all fixed moments of `Pj`
uniformly and makes its second/fourth moments tend to one/three.

From (5),

\[
 \operatorname{Cov}(BY)
       ={m-3\over m}Q^2+{2\over m}Q^3,
\]

so `||(BY)i||4=O(n^(1/4))` by (4) and degree-three hypercontractivity.
Also `||(BP)i||4=O(n^(1/4))`, since
`||B[S U]_i||4≤||B_i||1 max_j||Uj||4=O(1)`.

Let `Wcav=B[S(P²-1)]/sqrt2`. Expanding (10), the difference `W-Wcav`
contains `B[S U²]/sqrt2`, of L² norm `O(n^-1/2)`, and `sqrt2 B[P U]`.
Writing `Uj=U0−2sqrt2 SjGj/m`, with `U0=sqrt2 SᵀBS/m`, the latter has a
common-scalar term bounded by

\[
 \|U_0\|_4\,\|(BP)_i\|_4=O(n^{-1/4}),
\]

and a remainder bounded by
`C m^-1 ||B_i||1 max_j||Pj Gj||2=O(n^-1/2)`. Therefore

\[
 \sup_i\|W_i-W_i^{\rm cav}\|_2=O(n^{-1/4}).           \tag{11}
\]

## 5. Only the degree-seven tree survives

Let `Vj` be the degree-six Fourier projection of `Pj²`. Then

\[
 \sup_j\|(P_j^2-1)-V_j\|_2\longrightarrow0.           \tag{12}
\]

Here is a moment-counting proof, requiring no new contraction theorem.
The tree theorem gives `E(Pj²-1)²→2`. The second moment of `Vj` is an expansion
of four `T3` copies with the additional restrictions that copies 1 and 2
have disjoint nonroot label sets, as do copies 3 and 4. Repeat the proof of
Section 1 with these restrictions. All lower-label and nonempty-parity
patterns vanish. Of the three Gaussian whole-copy pairings, pairing 1 with
2 and 3 with 4 is prohibited. The other two survive, each with weight one
after the `1/sqrt2` normalization. Hence `E Vj²→2`, uniformly in `j`.
Orthogonality of Fourier degrees proves (12).

The remainders in (12) have degree at most four and are independent of `Sj`.
Apply (9) with `a_j=Bij`; since `Σ_jBij²=1`,

\[
 \sup_i\left\|W_i^{\rm cav}
       -{1\over\sqrt2}\sum_jB_{ij}S_jV_j\right\|_2\to0. \tag{13}
\]

The surviving degree-seven polynomial is precisely
`X_(T7,i)/sqrt8`, except for terms in which one of its seven distinct
nonroot labels equals the fixed output label `i`. Here `T7` is the full
binary depth-two tree beneath one external-root edge: the first child has
two children, each with two leaves. The excluded-output-label terms have
only `O(n^6)` distinct monomials, each coefficient `O(m^-7/2)` (only finitely
many tree-vertex assignments contribute to a fixed monomial). Their L² norm
is `O(n^-1/2)`.

The automorphism count is
`|Aut_o(T7)|=2·2·2=8`: exchange the two main branches and independently
exchange the two leaves in each. Combining (11)--(13),

\[
 \sup_i\left\|W_i-{X_{T7,i}\over\sqrt8}\right\|_2\to0. \tag{14}
\]

Likewise `G_i=X_(T1,i)`, and (10) gives
`Y_i-X_(T3,i)/sqrt2→0` in L² uniformly. The three rooted trees are
nonisomorphic. Section 1 therefore proves the explicit state evolution

\[
 \boxed{\quad
 (S_i,G_i,Y_i,W_i)\Longrightarrow
       (S,Z_0,Z_1,Z_2),
 \quad}                                              \tag{15}
\]

uniformly in the output coordinate, where `S` is Rademacher and the three
`Z` variables are independent standard normals, independent of `S`.

This proves the requested raw `H2` second-layer module. It is not a theorem
for arbitrary multilayer nonlinearities, nor a generic AMP universality
claim. Extending the cavity error control or the allowed tree family requires
new proofs, even though the injective odd-degree tree CLT itself applies to
every fixed tree satisfying its stated conditions.
