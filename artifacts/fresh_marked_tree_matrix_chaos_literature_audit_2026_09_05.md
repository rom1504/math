# Exact marked-tree derivative matrices: independent primary-source audit

Date: 2026-09-05. This note audits the operator-norm ingredient for the
bounded-operator-norm local-weight calculation. It does not assume AMP.

## 1. Precisely imported result

Primary source: A. S. Bandeira, K. Lucca, P. Nizic-Nikolac, and R. van
Handel, *Matrix Chaos Inequalities and Chaos of Combinatorial Type*,
December 24, 2024 version, [author PDF](https://web.math.princeton.edu/~rvan/chaosconf241224.pdf).

Theorem 2.1 compares a homogeneous square-free matrix chaos to its
independent-copy decoupling, up to a constant depending only on its degree.
Theorem 2.4 bounds the latter's expected operator norm by a degree-dependent
constant times a logarithmic dimension factor and the largest coefficient
sigma-flattening norm. A sigma-flattening assigns each random-variable
index to one of the two sides while keeping the original matrix row and
column indices on opposite sides. Standard Gaussian inputs satisfy the
stated distributional hypothesis. Theorems 2.1 and 2.4, the definition of
sigma-flattening, and the square-free hypothesis were read directly in the
primary PDF. No stronger matrix-chaos theorem is needed here.

## 2. Global tree coefficient tensor

Let `B=A/sqrt(n-1)` be an actual hollow symmetric signing with
`||B||op <= L`. A fixed fully marked rooted tree has `d` free vertices;
include its external root as an additional tensor index. Its
collision-allowed coefficient tensor has `d+1` indices and entries equal
to products of the `d` edge weights `B_uv`.

Every nontrivial flattening of this global tensor has operator norm
`O_(T,L)(1)`. Here is the independent direct check. Split the vertices
into two nonempty sets. Crossing edges give a compression of a tensor
product of copies of `B`, with norm at most `L^(E_cross)`. A vertex not
incident to a crossing edge contributes a `sqrt(n)` constant-vector
factor. Internal edges give diagonal multipliers of modulus at most
`(n-1)^(-1/2)`. If `I` denotes the number of crossing-isolated vertices,
the resulting bound is

`L^(E_cross) n^(I/2) (n-1)^(-E_internal/2)`.

Since the original graph is a connected tree, each component of the
internal-edge forest meets a crossing edge. Consequently
`I <= E_internal`, giving the claimed bound. This argument includes
flattenings that place the external root on either side.

Deleting all label collisions changes the global tensor by Hilbert norm
`O_T(1)`: there are `O_T(n^d)` forbidden tuples, each squared coefficient
is `O_T(n^-d)`. Although this global error is not vanishing, a bounded
error is all that the present operator-norm lemma requires. Every
flattening of the injective tensor therefore remains bounded. Fixed
symmetrization and tree-automorphism constants do not change the order.

## 3. Exact directional-derivative matrix

Write `X_T,i(N)` for the actual injective Gaussian tree polynomial, and
let `b_j` be row `j` of `B`. For `1 <= r <= d`, set

`J_T^(r)[i,j] = D_(b_j)^r X_T,i(N)`.

This matrix is a homogeneous square-free Gaussian chaos of degree
`q=d-r`; its remaining Gaussian coordinates are distinct because the
original tree kernel is injective. Its coefficient tensor is obtained
from the global tree tensor by contracting `r` slots with `b_j^tensor r`.

Let `V_r` be the matrix whose rows are `b_j^tensor r`. Then

`V_r V_r^T = Q^(circ r)`, where `Q=B^2`,

and `||V_r||op <= L`. Indeed `Q` is a correlation matrix, and Schur
multiplication by a correlation matrix contracts the operator norm of
self-adjoint matrices; iterate this to bound every Schur power by
`||Q||op <= L^2`.

For any sigma-flattening of the derivative coefficient tensor, place the
`r` contracted slots on the side that contains `j` before contraction.
This is one of the bounded global tree flattenings from Section 2.
The contraction is multiplication by `V_r` on that side, tensored with
identities on other indices. Its operator norm is bounded. Thus all the
sigma-flattenings required by the imported theorem are `O_(T,L)(1)`.

For `q >= 1`, the imported result now gives

`E ||J_T^(r)||op <= C_(T,L) log(2n)^(q/2)`.

For `q=0`, the same conclusion follows directly from the bounded
deterministic flattening. This is a statement about the exact injective
fields, not recursively substituted approximations.

## 4. Fixed operator moments

The expectation bound suffices to obtain every fixed operator moment.
For completeness, this upgrade is dimension-free. If `J` lies in one
homogeneous Gaussian chaos of degree `q`, the Ornstein--Uhlenbeck
semigroup satisfies `T_rho J=rho^q J`. Convexity of a norm gives

`||T_rho J||op <= T_rho(||J||op)` pointwise.

Scalar Gaussian hypercontractivity therefore implies

`|| ||J||op ||_p <= (p-1)^(q/2) || ||J||op ||_2`, for `p >= 2`.

In particular its fourth moment is bounded by a degree-dependent
constant times the square of its second moment. Paley--Zygmund applied
to `||J||op^2` then bounds the second-moment norm by a degree-dependent
constant times `E||J||op`. Combining the two inequalities proves a
polylogarithmic bound in every fixed `L^p`. In particular, for every
fixed `epsilon > 0`,

`|| ||J_T^(r)||op ||_p = O_(T,L,p,epsilon)(n^epsilon)`.

## 5. Derivatives of polynomial child responses

For a fixed polynomial `h` of a finite collection of exact child fields,
form `K_r[i,j]=D_(b_j)^r h(X_i)`, with `r >= 1`. The chain rule writes
this as a fixed sum of diagonal coefficient matrices times Hadamard
products of the matrices from Section 3. Every summand has at least one
positive-order directional derivative. The elementary compression
identity for a Hadamard product gives

`||M circ N||op <= ||M||op ||N||op`.

All fixed moments of each rooted child field are bounded; consequently
the maximum over `n` roots of a fixed polynomial coefficient is
`O_Lp(n^epsilon)` for every fixed positive epsilon, by taking a sufficiently
high but fixed moment before the union bound. Holder's inequality and
Section 4 thus give the same `n^epsilon` operator bound for `K_r`.

This closes the derivative-matrix norm ingredient. It does not alone
identify an energy limit: the separate local-weight proof must still use
own-root independence when integrating the uncontracted top spin, and
must retain the mixed cubic transport channel rather than Gaussianize it
away.

## 6. Audit of the one-level representation and weighted third chain

The completed proof in
`fresh_tree_matrix_derivative_and_weighted_third_chain_2026_09_05.md`,
Sections 4--5, also passes this independent audit.

The one-level field uses the actual child fields and thus preserves their
exact exclusion of their own Gaussian coordinate. Gaussian multiplication
expands the normalized child Hermite product into its highest Wick tensor
plus proper contractions. Full equal-child contractions are precisely
those canceled by the Hermite correction; full distinct-child contractions
vanish at the rooted moment rate. Proper contractions are bounded by the
fixed-root flattening estimate. Cross-child collision deletion has the
stated `O(n^-1/2)` Hilbert norm.

The transport of this input error is justified by the Gaussian creation
inequality in that artifact. A coefficient of the resulting Hermite
polynomial can receive contributions only from root coordinates with
exponent exactly one, since each input error excludes its own root.
There are at most `D+1` such choices. This proves the inequality with no
independence assumption between different roots. Fixed-degree Gaussian
derivative bounds then preserve the same rate in any fixed unit direction.

In the hard third-chain term, the stated trace order is correct. The
diagonal entry involving the last factor `K_3^T` sees exactly row `l` of
`K_3`, so differentiation with respect to `N_l` kills that contribution
exactly. When it instead hits a response diagonal, the maximum coordinate
derivative has size `polylog(n)/sqrt(n)` by uniform coordinate influence
and fixed-degree hypercontractivity at logarithmic moment order. All other
operators have polylogarithmic moments. This proves a vanishing normalized
trace, while avoiding the false claim that the original third-derivative
matrix itself has vanishing operator norm.

The independent conclusion remains scoped: this verifies the complete
single-third-derivative chain, not by itself all derivative splits or the
final Rademacher weighted-energy theorem.
