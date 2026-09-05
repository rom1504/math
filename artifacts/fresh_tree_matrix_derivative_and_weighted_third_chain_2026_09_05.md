# Matrix derivative bounds and the difficult weighted third-chain term

Date: 2026-09-05. This module is proved for Gaussian input under a fixed
normalized operator bound. It supplies a component of the larger weighted
unmarked-energy identity; it does not assert that every other component
or the Rademacher transfer has been completed here.

Let `B=A/sqrt(n-1)` be a symmetric hollow sign matrix with `||B||op≤L`.
Put `Q=B²`, write `b_j` for its rows, and let `N` be an independent standard
Gaussian input vector. The fields `X_T,i` are the normalized injective
fully marked Gaussian tree chaoses. Trees and polynomial degrees are fixed.

## 1. Global tree coefficient flattenings

Regard the external root as another tensor index. For a tree having `d`
nonroot vertices, let

\[
 K_T(i,a_1,\ldots,a_d)
\]

be its coefficient tensor before symmetrization. Normalization constants
depending only on the tree are suppressed in this section.

**Global flattening bound.** Every nontrivial bipartition of these `d+1`
indices has matrix operator norm `O_{T,L}(1)`.

First allow all labels to coincide. The crossing edges of a bipartition
form a bipartite graph. For each incident vertex, copying its label into
all its crossing-edge slots is an isometry. The crossing-edge operator is
therefore a compression of the tensor product of the corresponding `B`
matrices, with operator norm at most `L^(number of crossing edges)`.
Each vertex isolated in that crossing graph contributes a factor
`sqrt(n)`, while each internal edge multiplies a row or a column by a
function of modulus at most `1/sqrt(n-1)`.

If `I` counts crossing-isolated vertices and `E_int` counts internal
edges, the bound is

\[
 L^{E_{\rm cross}} n^{I/2}(n-1)^{-E_{\rm int}/2}.
\]

The whole vertex graph is a connected tree. Every component of the
internal-edge forest contains a crossing-incident vertex, because both
sides of the cut are nonempty. Thus `I≤E_int`, proving the bound.

Deleting all collisions changes the tensor only on `O_T(n^d)` entries,
each of squared modulus at most `(n-1)^(-d)`. Its total Hilbert--Schmidt
change is `O_T(1)`. Every flattening operator change is at most this.
Symmetrizing the fixed number of nonroot slots likewise preserves the
bound. This proves the claim for the actual injective tensor.

Notice the difference from fixing the external root: fixing it yields
an extra `n^-1/2` in every proper flattening. The global statement needs
only `O(1)` and has a correspondingly larger collision Hilbert norm.

## 2. A precise primary matrix-chaos theorem and its application

The imported result is the iterated noncommutative Khintchine inequality
and its square-free decoupling step, Theorems 2.1 and 2.4 of
[Bandeira--Lucca--Nizic-Nikolac--van Handel, *Matrix Chaos Inequalities and
Chaos of Combinatorial Type*](https://web.math.princeton.edu/~rvan/chaosconf241224.pdf).
For a fixed-order square-free Gaussian matrix chaos, they bound the expected
operator norm by a fixed power of `log n` times the largest oriented
coefficient flattening. Here an oriented flattening keeps the two matrix
indices on opposite sides and assigns each stochastic index to either
side. We independently inspected both the theorem hypotheses and its
iteration proof in Appendix A.2; no independence beyond the input
coordinates is assumed.

For `1≤r≤d`, form the matrix

\[
 J_T^{(r)}(i,j)=D_{b_j}^{,r}X_{T,i}.                 \tag{1}
\]

It is a homogeneous square-free Gaussian matrix chaos of order `d-r`.
The derivative contracts `r` input slots with `b_j^(tensor r)`. The row
tensors of this contracting map have Gram matrix `Q^(circ r)`, whose
operator norm is at most `||Q||op≤L²`: Schur multiplication by a correlation
matrix contracts the self-adjoint operator norm.

Every oriented flattening of the new coefficient tensor arises by putting
all `r` contracted slots on the side that contains the matrix index `j`
in the original global tree flattening, and composing with this bounded
contracting map. Section 1 therefore bounds all of them uniformly.
The remaining random indices stay pairwise distinct, so the cited
square-free theorem applies literally; no non-square-free Wick extension
is necessary. If `r=d`, the assertion is just the deterministic flattening
bound.

Consequently, for every fixed finite `p`,

\[
 \boxed{\quad
 \|\,\|J_T^{(r)}\|_{op}\,\|_p
       \le C_{T,L,p}(\log(n+2))^{C_T}.
 \quad}                                                    \tag{2}
\]

For clarity, the passage from expectation to fixed moments is
dimension-free. The Gaussian noise semigroup and Jensen give the usual
hypercontractive bound for the norm of a Banach-valued homogeneous chaos.
It bounds the fourth moment by a fixed multiple of the second. Applying
Paley--Zygmund to its squared norm then bounds the second moment by a
fixed multiple of the first; higher fixed moments follow. Thus the
expectation theorem gives (2), not merely a probability bound too weak
for subsequent Holder estimates.

## 3. Directional derivative matrices of child polynomials

Let `P` be any fixed polynomial in a finite family of exact old fields at
root `l`, and define, for `r≥1`,

\[
 K_r(l,j)=D_{b_j}^{,r}P(X_l).
\]

The chain rule expresses this matrix as a finite sum of random diagonal
matrices times Hadamard products of the matrices (1), each differentiated
factor receiving at least one derivative. The elementary inequality

\[
 \|M\circ N\|_{op}\le\|M\|_{op}\|N\|_{op}
\]

follows by diagonal compression of `M tensor N`. All random diagonal
coefficients are maxima of fixed-degree Gaussian polynomials of bounded
moments. Hypercontractivity at order proportional to `log n`, followed by
a union bound, gives polylogarithmic bounds for their fixed moments.
Hence

\[
 \boxed{\quad
 \|\,\|K_r\|_{op}\,\|_p\le C(\log(n+2))^C.
 \quad}                                                    \tag{3}
\]

If `P(X_l)` is a child feature `h_T(X_l)`, it excludes the own Gaussian
coordinate `N_l` exactly. Directional derivatives preserve this exclusion:

\[
 \partial_{N_l}K_r(l,j)=0\quad\hbox{for every }j.       \tag{4}
\]

This is a rowwise assertion. Other rows of `K_r` can depend on `N_l`.

## 4. A controlled one-level recursive representation

Let `h_T` be the normalized Hermite product of the child-tree coordinates
for a non-edge tree `T`. Define

\[
 Y_{T,i}=\sum_l B_{il}N_l h_T(X_l),\qquad Y_{edge}=G.   \tag{5}
\]

Here the child fields are the actual injective fields, not recursively
substituted approximations. Then for every fixed finite `p`, uniformly in
the output root,

\[
 \boxed{\qquad \|Y_{T,i}-X_{T,i}\|_p=O_{T,L,p}(n^{-1/2}).\qquad} \tag{6}
\]

To see this, expand the child Hermite product by the Gaussian multiplication
formula. Full matching child contractions are canceled by the displayed
Hermite polynomials. Any surviving proper child contraction has Hilbert
norm `O(n^-1/2)` by the fixed-root proper-flattening theorem. Full
contractions of distinct child types, and the variance normalization
errors of identical types, are also `O(n^-1/2)`: their two-tree moment
graphs either match or have a nonempty parity graph, bounded using
`beta(A)≤n||A||op=O(n^{3/2})`. A fixed multiplication expansion therefore
leaves the top Wick tensor product plus `O_Lp(n^-1/2)` error. Both terms
exclude `N_l`.

For polynomials `e_l` independent of `N_l`, of degree at most `D`, the
Gaussian own-root divergence inequality is

\[
 E\left|\sum_l c_l N_l e_l\right|^2
 \le(D+1)\sum_l c_l^2 E|e_l|^2.                       \tag{7}
\]

Expand in the normalized Gaussian Hermite basis. A resulting coefficient
can arise only from a coordinate with exponent exactly one, and it has
at most `D+1` such root choices. Cauchy--Schwarz proves (7).
Apply it with `c_l=B_il` to transport the child-product error. Finally the
top transported tensor differs from the desired injective parent only
on cross-child collisions and collisions with the fixed external root;
their fixed-root Hilbert norm is `O(n^-1/2)`. Hypercontractivity proves
the full statement (6).

All fixed unit-direction derivatives of the error in (6) satisfy the
same order of `L^p` bound, by the Gaussian derivative formula for a
fixed-degree polynomial. The same holds for ordinary coordinate
derivatives. In particular replacing the response inputs `X` by the
one-level fields `Y` preserves their uniform coordinate-influence bounds.

For fixed polynomial responses, this replacement is legitimate in a
weighted normalized energy such as
`n^-1 E M(X) Z dot B F(X)`: bounded operator norm, Holder, (6), and the
uniform fixed moments of `Z` control the difference by `o(1)`.
One may therefore perform the following calculation with responses
evaluated at `Y`, while retaining exact child exclusion in (5).

## 5. The single third-derivative chain vanishes with a local weight

This is the term not controlled by a naive pointwise directional bound.
Let `M,f` be fixed polynomial local responses evaluated at the one-level
field family `Y`, set `D_M=diag(M(Y_j))`, `D_f=diag(f(Y_i))`, and put
`C=B D_M B`. We prove

\[
 \boxed{\quad
 \frac1n E\sum_{i,j} C_{ij}f(Y_i)D_{b_j}^3Y_{T,i}=o(1).
 \quad}                                                    \tag{8}
\]

The assertion is zero for the edge tree. For a non-edge tree use (5) and
differentiate the explicit top Gaussian coordinate:

\[
 [D_{b_j}^3Y_{T,i}]_{i,j}
   =B\operatorname{diag}(N)K_3+3B(B\circ K_2),        \tag{9}
\]

where `K_r(l,j)=D_{b_j}^r h_T(X_l)`. The term with a derivative on the
top coordinate is small in Frobenius norm:

\[
 \|B(B\circ K_2)\|_F
 \le\frac L{\sqrt{n-1}}\|K_2\|_F
 \le L\sqrt{n/(n-1)}\|K_2\|_{op}.
\]

It therefore has polylogarithmic fixed moments by (3). The other matrix
in its energy pairing, `D_f C`, has Frobenius moments bounded by
`sqrt(n)` times a polylogarithm. Dividing by `n` proves this contribution
is `O(polylog(n)/sqrt(n))`.

For the root-not-differentiated term, rearrange the sum exactly as

\[
 \frac1n\sum_l E\left[
 N_l\{B D_f B D_M B K_3^T\}_{ll}\right].              \tag{10}
\]

Integrate `N_l` by parts. The derivative of the last matrix contributes
zero: in this diagonal entry it differentiates only row `l` of `K_3`,
and (4) applies exactly. The two remaining contributions differentiate
`D_f` or `D_M`. Uniform coordinate influences and hypercontractivity give

\[
 \left\|\max_i|\partial_{N_l} f(Y_i)|\right\|_p
 +\left\|\max_j|\partial_{N_l} M(Y_j)|\right\|_p
 \le\frac{C(\log(n+2))^C}{\sqrt n},                 \tag{11}
\]

uniformly in `l`. Undifferentiated diagonal maxima, and `||K_3||op`, have
polylogarithmic moments. Bound each differentiated diagonal entry of
(10) by the product of the corresponding operator norms and use Holder.
After the normalized sum over `l`, the bound is again
`O(polylog(n)/sqrt(n))`. This proves (8).

The calculation uses an exact own-coordinate exclusion, not a claim that
the whole third-derivative matrix has small operator norm. That latter
claim is false even for a star: it can contain `B diag(N_l G_l) Q^(circ3)`.
Polylogarithmic operator control plus the own-coordinate integration by
parts is the appropriate combination.

## 6. Scope of the closure

Equation (8) closes the all-derivatives-on-one-response third-chain class
in a weighted cubic unmarked-energy calculation. The all-derivatives-on-one-
response second-chain classes admit simpler Frobenius estimates: one
non-edge directional second derivative together with an edge derivative
has matrix Frobenius moments `O(1)`, because `sum Q_ij²=O(n)`; pairing it
with a weighted `B diag(M) B` costs only `O(n^-1/2)` after normalization.

The full weighted theorem additionally requires auditing derivative
splits between the two response roots, the surviving all-edge term, and
the complete sign-input replacement. Those steps belong to the separate
weighted-energy proof. They are not silently inferred from this module.
