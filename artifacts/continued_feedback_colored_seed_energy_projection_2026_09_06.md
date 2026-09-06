# Finite-vector independent-seed coherent energy projection

Date: 2026-09-06. Status: direct extension of the retained-coherent scalar
proof, submitted for independent audit. This is not an assertion for
all canonical old-tree coordinates or all finite-query algorithms.

## 1. Exact statement

Fix `d` before matrix order. Let `S^1,...,S^d` be independent Rademacher
seed vectors, let `B=A/sqrt(n-1)` be an arbitrary symmetric hollow
signing with fixed `||B||op<=L`, and set

```math
G^a=BS^a,\qquad Q=B^2,\qquad
G_i=(G_i^1,\ldots,G_i^d).
```

Let `f:R^d→R` be bounded, Gaussian-a.e.-continuous, and globally odd:
`f(−g)=−f(g)`. Let `H` be bounded, Gaussian-a.e.-continuous, and globally
even. Write its product normalized Gaussian Hermite expansion as

```math
f(g)=b\cdot g+r(g),\qquad
r(g)=\sum_{\substack{\alpha\in\mathbb N^d\\|\alpha|\ge3\ {m odd}}}
 f_\alpha h_\alpha(g),\qquad
b_a=\mathbb E[G_a f(G)].
```

Retain the literal first-chaos return and the actual nonlinear channel:

```math
V_i=\sum_{a=1}^d b_a(QS^a)_i,\qquad
Z=Br(G),\qquad
T=B\left[\sum_{\alpha}f_\alpha^2Q^{\circ|\alpha|}\right]B,
\qquad\sigma_i^2=T_{ii}.
```

Let `psi` be odd `C²` with `psi,psi',psi''` bounded, and define

```math
C_i=H(G_i)\psi(V_i+Z_i),
\quad c_i^0=H(G_i)\mathbb E_N\psi(V_i+\sigma_iN),
\quad a_i=\mathbb E_{S^1,\ldots,S^d,N}
                      H(G_i)\psi'(V_i+\sigma_iN).
```

Then the proposed conclusion is exactly

```math
\frac{\mathbb E[C^{\mathsf T}BC]}{2n}
=\frac{\mathbb E[(c^0)^{\mathsf T}Bc^0]}{2n}
 +\frac{\mathbb E[(c^0)^{\mathsf T}B D_aZ]}{n}
 +\frac{\operatorname{Tr}(B D_aT D_a)}{2n}+o(1).                (1)
```

The local tested identity is

```math
\frac1n\mathbb E\sum_i H(G_i)(V_i+Z_i)\psi(V_i+Z_i)
=\frac1n\sum_i\mathbb E_{S^1,\ldots,S^d,N}
 H(G_i)(V_i+\sigma_iN)\psi(V_i+\sigma_iN)+o(1).                 (2)
```

Every displayed expectation containing the `S^a` keeps their actual
Boolean law. In particular, `V_i` is not an independent Gaussian.

If `H>=0`, `|f|+H<=1`, `||psi||infinity<=1`, and `y psi(y)>=0`, the
two feasible endpoints `±f(G)+C` consequently satisfy

```math
\Lambda(B)\ge j_n+
\left|\frac{\|b\|_2^2}{2}\frac{\operatorname{Tr}(B^3)}n
       +e_n^0+k_n+t_n\right|-o(1),                            (3)
```

where `j_n` is the right side of (2), and the remaining three terms
are the three terms on the right of (1). As in the scalar theorem,
no uniform positive sign or magnitude is asserted for the common term.

## 2. Colored kernel verification

View the input coordinates as pairs `(a,u)`, with color `a` and vertex
`u`. The linear old row of color `a` is row `b_j` of `B`, embedded in
its color block of `R^{dn}`. Distinct colors are orthogonal exactly.

For a Hermite multiindex `alpha`, let `p=|alpha|`. Its Gaussian kernel
is the normalized symmetrized product of `alpha_a` old rows of each
color. Its covariance between roots `i,j` is `Q_ij^p`; distinct
multiindices have zero covariance, even when their total degrees agree.
The factorial normalization is the standard product-Hermite one.
After transport by `B`, the channel covariance is `B Q^{circ p}B`.

Every proper fixed-root flattening is a finite sum of maps of the form

```math
M_{\mathrm{left}}^{\mathsf T}
\operatorname{diag}(b_i)M_{\mathrm{right}}.
```

Each group Gram is a Schur power of `Q` (or zero if its color matching
is impossible), so its operator norm is bounded. The middle diagonal
has norm `1/sqrt(n-1)`. Thus the flattening bound remains
`O_{d,p,L}(n^{-1/2})`. Repeated marked slots can occur only in the same
color; the scalar collision proof applies to each color choice, with
only finitely many choices at fixed degree.

When a whole transported leaf of total degree `q>=3` matches a selected
set of center slots, its colors must match. If they do, the scalar
covariance factor is still `Q_ij^q`; otherwise the term is zero.
Consequently its root-summed matrix is still `Q^{circ q}B`, which is
entrywise `O(n^{-1/2})` by the same absolute-row-sum estimate.

## 3. Coherent Walsh factors and the surviving star

The coherent fields consist of finitely many block matrices on the
joint seed space: one embedded copy of `B` for each old color, and
embedded copies of `Q` for its coherent returns. These matrices have
bounded operator norms and bounded row and column Euclidean norms.
Their finite entrywise products obey the same row/column
Cauchy--Schwarz bounds as in the scalar proof.

For a full nonlinear contraction into coherent row factors, restrict
each factor matrix to the relevant color block. Its contraction with
an old row is then a matrix `B N^T` with bounded operator norm. The
flat leading root factor `B_ik` and two column Euclidean bounds again
give an entrywise `O(n^{-1/2})` bare-star bound. All impossible color
matchings vanish exactly.

Thus the scalar proof's entire classification is unchanged: partial
non-star components and multiple stars have two gains; one star with
at least one transported leaf is small in operator norm; a bare star
with only coherent leaves is the exact retained `c0`--`D_aZ` cross;
and any extra positive-degree coherent bridge makes it energy-negligible.
Multiple whole covariance factors also have vanishing flat-`B` trace.
This proves the polynomial version of (1)--(2).

Approximate bounded `H` in the universal product Gaussian marginal of
the `d` old coordinates; independence of the seed colors makes that
marginal standard. The scalar bounded-response argument, variance bins,
and exact-coefficient Gaussian coupling then apply in fixed dimension.
No dimension grows with matrix order. The old input self-energy is
`||b||² Tr(B³)/(2n)+o(1)` because each nonconstant Hermite multiindex
has covariance `Q^{circ|alpha|}`, and all degrees at least three vanish
in the flat-`B` trace. This proves (3).

## 4. Precise scope boundary at canonical old-tree fields

The independent-color extension is genuinely multicoordinate, but it
does not encompass an arbitrary finite canonical old-tree family from
the rich response architecture. Those old coordinates are higher
degree polynomials of a COMMON seed. Their one-root Gaussian law and
pair covariance `delta_TU Q` do not make their marked kernels equal to
independent colored linear rows.

Two proof steps would need new verification there:

1. The literal coherent return is `sum_T b_T B X_T`, or approximately
   `Q[S K_f(X)]`, not a bounded-op LINEAR transformation of independent
   seeds. The exact coherent Walsh root-map lemma from products of
   matrix rows is not automatically available for nonlinear functions
   of this return and the old trees.
2. In a mixed bare-star contraction, marked slots can follow nontrivial
   old-tree branches with different original input degrees. Whole local
   Gaussian covariance `Q^{circ k}` alone does not determine those
   cross-branch contractions. The scalar/colored factorization yielding
   `Q^{circ q}B` must be proved from the actual tree tensors or replaced
   by another small-matrix bound.

The previously proved full nonlinear old-forest theorem supplies proper
flattenings and several linear-tested contractions. That is useful
input, but it does not by itself discharge these two obligations for
an arbitrary nonlinear function of the coherent return. This note
does not claim growing-depth universality or a general GFOM law.
