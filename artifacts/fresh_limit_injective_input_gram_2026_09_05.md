# Operator orthogonality of injective marked tree inputs

Date: 2026-09-05. This is an exact-kernel statement. It does not identify a
raw product of previously computed fields with an injective kernel in
operator norm.

Audit status: the variational agent independently checked Sections 8--9,
including arbitrary raw-product collision patterns, the doubled-Frobenius
bound, and the weighted Gaussian-Sobolev approximation. The director
independently reconstructed the operator-Gram mechanism. This note keeps
the stronger smooth-response result distinct from a full unmarked CLT or
the still-unproved general odd-response `Q`-trace lower bound.

## 1. Statement and normalization

Let `A` be a symmetric hollow sign matrix, `m=n-1`, and suppose

`beta(A):=max_{x,y in {+-1}^n}|x^T A y| = O(n^{3/2})`.

This follows from the original low-cap hypothesis. Let `R` be a finite
rooted tree whose root has even degree and every nonroot vertex has odd
degree. The single-vertex tree is allowed. Its vertex count `d` is odd.
Let `a(R)` be the number of root-preserving automorphisms. Define

`V_{R,i} = [sqrt(a(R)) m^{(d-1)/2}]^{-1}`

`  * sum_{f injective, f(root)=i}
      [product_{uv in E(R)} A_{f(u),f(v)}]
      [product_{v in V(R)} S_{f(v)}]`.                            (1)

Unlike an output tree field, the root in (1) itself carries a spin.
Its even degree is intentional: adjoining a new external root and its
edge produces exactly an odd-degree fully marked tree of the previous
hierarchy. The input in (1) is own-spin times an own-spin-free polynomial.

For any fixed finite list of pairwise nonisomorphic rooted trees `R`,

`|| E[V_R V_U^T] - 1_{R isomorphic U} I_n ||op = O(n^{-1/4})`.    (2)

The constant depends on the fixed trees and the bound on
`beta(A)/n^{3/2}`, not on `n`. In particular the stacked family has
asymptotically identity covariance in operator norm. The variables have
mean zero, since all their homogeneous degrees are odd.

## 2. Pairing is exact, not a Gaussian approximation

Each summand in (1) is a multilinear monomial with `d` distinct spin
indices. Two such monomials have nonzero Rademacher inner product if and
only if their sets of indices coincide. Consequently inputs with different
values of `d` are exactly orthogonal. For equal `d`, covariance is a finite
sum over bijections between the two formal vertex sets. There are no
lower-order repeated-spin partitions: within-copy injectivity makes each
label occur at most twice, so every surviving label occurs exactly twice.

Overlay the two trees using such a bijection. Delete every edge appearing
twice; the remaining simple graph is the parity-edge graph `H`. Edge signs
in the covariance summand are exactly `product_{uv in E(H)} A_{uv}`.
The deleted paired edges still supply their normalization powers of `m`;
they are not removed from the denominator.

## 3. Diagonal entries

At equal output roots `i`, the root must be matched to the root. Every
nonroot quotient vertex combines two odd degrees, and the quotient root
combines two even degrees. Thus every vertex of `H` has even degree.

If `H` is empty, both trees have exactly the same edge set on the quotient
vertex set. The bijection is therefore a rooted isomorphism. This is
possible exactly when `R` and `U` are rooted-isomorphic. There are `a(R)`
such bijections. Each contributes `(m)_{d-1}/m^{d-1}=1+O_R(1/n)` after
normalization, where `(m)_k` denotes a falling factorial. The factors
`a(R)` cancel the automorphism normalization in (1).

If `H` is nonempty, its Eulerian parity prevents it from consisting only
of edges incident to the fixed root: a nonempty simple star has odd-degree
leaves. Thus `H` has an edge with both endpoints free. Average over these
two labels after fixing all the other labels. Since `H` is simple, all
other incident edge factors split into one bounded function of each
endpoint. The absolute average is at most `beta(A)/n^2`. Dropping
injectivity changes the normalized average by `O_R(1/n)`. Hence, uniformly
in `i`, every nonempty-parity diagonal contribution is `O(n^{-1/2})`.

This proves the diagonal part of (2), with the stronger rate `n^{-1/2}`.

## 4. Off-diagonal entries have exactly two odd parity vertices

Now let the two output labels be distinct, `i != j`. The root of the first
copy must match a nonroot vertex in the second, and vice versa. Therefore
in `H` the vertices labeled `i` and `j` have odd degree. Every other vertex
has even degree. In particular `H` is nonempty.

There are `d-2` free labels after fixing `i,j`, while the product
normalization is `m^{-(d-1)}`. Up to a fixed normalization constant and an
entrywise error `O_R(n^{-2})`, each pattern matrix therefore has the form

`K_H(i,j) = n^{-1} E_{free labels} product_{uv in E(H)} A_{uv}`.  (3)

The labels in the expectation in (3) are unrestricted independent uniform
labels. Indeed injectivity exclusions have probability `O_R(1/n)` and
the prefactor is `O(1/n)`. Replacing `m` by `n` has the same error. An
entrywise `O(n^{-2})` matrix has operator norm `O(n^{-1})`, so it remains
to bound the matrices in (3). Set their diagonals to zero if needed;
this alters operator norm by `O(1/n)`.

There are two cases.

### 4.1 Only the edge `ij` survives

If the only edge of `H` is `ij`, all other formal vertices are isolated,
and (3) is simply `A/n`. The elementary spectral bootstrap gives

`||A||op^2 <= beta(A)`.

To check this inequality, take an eigenvector `v` for an eigenvalue of
maximum magnitude and use a randomized sign vector of mean
`v/||v||infinity`. Then
`beta(A)>=|lambda| ||v||1/||v||infinity`. The eigenvector equation and
`|A_ij|<=1` give `|lambda| ||v||infinity<=||v||1`.
Therefore `||A/n||op=O(n^{-1/4})`.

### 4.2 Some other edge survives

Square the Frobenius norm of (3). It is exactly the normalized signed
homomorphism density of two copies of `H` with only their distinguished
vertices `i,j` glued together. All other formal vertices in the copies
remain separate and are independently summed. If the edge `ij` is present,
its two factors cancel. Every edge other than `ij` is present exactly once
in this doubled graph: it meets at least one free vertex that is not glued.
Thus the parity graph of this Frobenius-square expression is nonempty.

All labels, including `i,j`, are now averaged. Choose any surviving edge
and fix the other labels. Exactly the same bilinear bound as in Section 3
gives

`||K_H||F^2 <= beta(A)/n^2 = O(n^{-1/2})`.                       (4)

If one first enforces `i!=j`, the excluded diagonal contributes only
`O(1/n)`; either convention gives the same conclusion. Thus
`||K_H||op <= ||K_H||F=O(n^{-1/4})`.

There are finitely many bijection patterns, depending only on the fixed
trees. Summing their bounds and combining with Section 3 proves (2).

## 5. What this does and does not yet give

The result is stronger than rowwise moment orthogonality: it controls any
deterministic linear combination of output roots, uniformly in that
coefficient vector. It uses no bounded-operator assumption beyond the
spectral bootstrap supplied by the low Boolean cap, and no Gram or
all-power coherence.

The exact fully injective output tree field is not quite `B V_R`: terms
where an internal descendant has the output label must be excluded.
Those collision fields have small rowwise `L^2` norm, but (2) by itself
does not establish a small covariance operator for that error. Likewise
the raw product `S_i h_R(X_i)` must not be substituted for (1) solely on
the basis of rowwise Gaussian convergence. These are separate operator
comparison problems. The purpose of (2) is to bank the exact injective
input Gram statement before attempting those extensions.

## 6. The output-root collision has small covariance operator

There is a direct comparison for the exact output tree itself. Let `T`
be obtained from `R` by adjoining a new external root and its single edge,
so the number of marked vertices of `T` is `d=|V(R)|`. Let `X_T` be its
fully injective output field. Then

`B V_R = X_T + E_T`,

where `E_T,i` consists of precisely those terms in which a nonroot vertex
of `R` has label `i`. The raw root of `R` cannot have label `i`, since its
outer factor would be `B_ii=0`. Since `R` itself is injectively labeled,
at most one of its vertices can have the output label. There are only
`d-1` possible choices of that vertex.

Each choice remains a homogeneous multilinear polynomial with `d`
distinct spin indices, including `i`. Its coefficient is a product of
`d` matrix factors, and thus has absolute value at most a fixed constant
times `m^{-d/2}`. Identifying the output root with its chosen descendant
can create a double edge or a longer cycle; neither changes this absolute
bound or the distinctness of the spin indices.

At equal roots, nonzero covariance again requires equality of the two
spin sets. There are at most `C_T n^{d-1}` assignments after the matching
pattern is fixed. Hence

`E E_{T,i}^2 <= C_T n^{d-1} m^{-d} = O_T(n^{-1})`.

At distinct roots `i,j`, equality of the spin sets forces both labels to
occur in each set, leaving at most `d-2` free labels. Therefore

`|E E_{T,i} E_{T,j}| <= C_T n^{d-2} m^{-d} = O_T(n^{-2})`.

These statements include the finite sum over the possible collided
vertices. Maximum absolute row sum consequently gives

`||E[E_T E_T^T]||op = O_T(n^{-1})`.                              (5)

No sign cancellation or operator hypothesis was used in (5).

If in addition `||B||op<=L` is fixed, (2), (5), and covariance
Cauchy--Schwarz yield

`||E[X_T X_U^T] - 1_{T isomorphic U} Q||op = o(1)`.              (6)

Indeed `Cov(BV_T,BV_U)=B Cov(V_T,V_U) B`, and each collision cross term
has operator norm at most the geometric mean of the two corresponding
covariance operator norms. The latter are bounded for `BV_T` and
`O(1/n)` for `E_T`. Thus the collision cross terms are `O_L(n^{-1/2})`,
while the main Gram error is `O_L(n^{-1/4})`.

Equation (6) is an actual covariance-operator theorem for the old output
fields, not an inference from their rootwise central limit theorem. It
still does not identify covariance operators of nonlinear raw products
of those output fields; those require a separate Wick-product comparison.

## 7. Independent audit of the two-fixed-root moment formula

There is also a uniform finite-dimensional statement that does not follow
from operator covariance alone. Fix two distinct output labels `i,k` and
a fixed list of old output trees at either root. Their joint moments are
those of a Gaussian family with within-root covariance `I` and cross-root
covariance `Q_{ik} I`, up to

`O_T(beta(A)/n^2 + 1/n)`.                                       (7)

The moment order and trees are fixed. No convergence of `Q_{ik}` is
required; the statement compares with the indicated possibly degenerate
Gaussian at each order. Here is the graph classification independently
checking the proposed formula.

Let `D` be the total number of marked vertex positions across the copies.
The normalization is `m^{-D/2}`. Nonzero Rademacher terms have every spin
label repeated an even number of times. Terms with fewer than `D/2` free
spin labels cost `O(1/n)` by counting. This includes a spin block assigned
either fixed output label: it consumes a pairing block without a free
label. Thus a leading pattern pairs all marked positions into `D/2`
distinct free labels, all different from `i,k`.

Each free parity-graph vertex combines two odd degrees, so is even. If
there is any free-free parity edge, fixing the other labels and using the
bilinear cap kills the contribution at cost `O(beta(A)/n^2)`, uniformly
in both fixed roots. If there is no such edge, every free parity vertex
is adjacent to neither root or to both: its possible degree is at most
two and must be even. Therefore the parity graph consists of `p` disjoint
length-two paths `i--v--k`.

For clarity, the full quotient graph can have doubled edges in addition
to those paths. It has `V=D/2+2` vertices and
`E=D/2+p` distinct edges. If `p>0`, it is connected: all original copies
touch one of the two roots and the paths join the roots. Its cycle rank
is consequently `E-V+1=p-1`. The `p` displayed paths already have exactly
that cycle rank, so no additional doubled-edge cycle or connecting route
can occur. If `p=0`, all edges are doubled and `E=V-2`; the quotient has
exactly two tree components, each containing one root.

More directly, every displayed path edge incident to a fixed root is a
top edge of one original copy. Its midpoint contains exactly two marked
positions, so it identifies the top vertices of one copy from each root.
Every child edge at that midpoint must be doubled, and the second copy
is already determined by those two positions. Within-copy injectivity
therefore propagates this same pair through all descendants, forcing a
rooted isomorphism of the entire two trees. Root-top edges that are
doubled instead identify two copies at the same root and have the same
whole-copy propagation. The cycle count above excludes any extra route
that could reconnect the resulting paired components.

Thus the surviving patterns are precisely Wick pairs of isomorphic tree
copies. Each cross-root pair contributes
`m^{-1}sum_v A_{iv}A_{kv}=Q_{ik}`; a same-root pair contributes one.
All remaining paired branch edges square to one, and their normalized
injective label counts contribute `1+O(1/n)`. The usual automorphism
normalizations cancel the isomorphism counts. Distinctness constraints
between different paired components alter the bounded normalized sum by
`O(1/n)`. This proves (7), including the potentially hostile many-cross-pair
case `p>=2` and mixed nonisomorphic trees.

## 8. Direct covariance-operator theorem for raw even responses

The same parity argument proves a stronger statement without any
operator-level Wick-product substitution. Let `a` be a fixed jointly even
polynomial in a fixed old tree family, and put

`U_i=S_i a(X_i)`.

Then, under the same low-cap hypothesis,

`|| E[U U^T] - E[a(Z)^2] I ||op = O_a(n^{-1/4})`,                (8)

where `Z` is the independent standard Gaussian old-tree family. More
generally, for two even polynomials `a,b`, their cross covariance is
`E[a(Z)b(Z)] I+o_op(1)`.

Expand a covariance entry into polynomial monomials, and then into their
old tree copies. Consider distinct roots `i,j`. Add the two explicit spin
positions `S_i,S_j`. If `D` is the total number of internal marked tree
positions, the normalization is `m^{-D/2}`. Every nonzero Rademacher term
has even label multiplicities. A leading term pairs all `D+2` spin
positions, including the two explicit endpoint positions. The endpoint
`S_i` must be paired with an internal marked position of a tree at root
`j`: every tree at root `i` is exactly own-root-free. Likewise `S_j` is
paired with an internal position of a tree at root `i`. The other pairs
have distinct free labels, leaving exactly `D/2-1` free labels.

Since each monomial of `a` and `b` contains an even number of old-tree
factors, the total graph degree at each of its fixed roots is even.
Identifying that root with the marked internal position from the other
root adds an odd degree. Therefore the parity-edge quotient has odd
degree exactly at `i,j`, and even degree everywhere else. Pairings within
one response, including collisions between two of its different tree
copies, are allowed in this argument: their degree sum is still even.
They do not need to be removed or approximated.

As in Section 4, each leading off-diagonal pattern has scale `1/n` and is
`1/n` times a fixed signed graph density. The quotient has a nonempty
parity graph with odd vertices `i,j`. If only the edge `ij` survives,
its matrix has operator norm `O(n^{-1/4})`. Otherwise the doubled
Frobenius argument in (4) gives that same rate. If fewer spin-label blocks
occur, or an additional block is assigned one of the two fixed labels,
there is one fewer free sum. Its entry is `O(n^{-2})`, hence its operator
norm is `O(n^{-1})`. There are only finitely many partition patterns for
the fixed polynomials.

On the diagonal, the explicit own spin squares to one, and uniform rooted
moment convergence gives
`E a(X_i)b(X_i)=E a(Z)b(Z)+O(beta(A)/n^2+1/n)`.
Together these facts prove (8) and its cross-covariance version. In
particular this proof already includes every raw polynomial collision;
it never substitutes an injective input for a raw product.

## 9. Smooth extension by a weighted own-free Sobolev estimate

Let `g` be a fixed smooth function of the old finite family, with bounded
derivatives or polynomially bounded derivatives and all the fixed moments
used below finite. The exact own-free endpoint identity gives, for `i!=j`,

`E[S_i g(X_i) S_j g(X_j)] = E[D_j g(X_i) D_i g(X_j)]`.

Consequently, for every deterministic unit vector `w`,

`E (sum_i w_i S_i g(X_i))^2`

`<= max_i E g(X_i)^2 + n max_{i!=j} E[D_j g(X_i)]^2`.            (9)

Indeed `sum_i |w_i| <= sqrt(n)`, so the absolute off-diagonal weight sum
is at most `n`. Direct coefficient counting and hypercontractivity for
the fixed old-tree kernels give
`max_ij ||D_j X_i||_4 <= C_X/sqrt(n)`. The chain rule, with its discrete
Taylor remainder, therefore yields

`limsup n max_ij E[D_j g(X_i)]^2
 <= C_X (E ||grad g(Z)||^4)^{1/2}`.                              (10)

For precision, each old field is affine in the varied spin. Replace the
gradient along its central-difference segment by the gradient at the
actual old field. The segment displacement has every fixed moment
`O(n^{-1/2})`; its contribution to the derivative is `O_{L^2}(n^{-1})`
for fixed `g`. The leading product is bounded by fourth-moment
Cauchy--Schwarz. Uniform rooted moment convergence then gives (10).
All approximation functions are fixed before `n` tends to infinity.

Taking the supremum over `w` in (9) proves

`limsup || Cov(S g(X)) ||op
 <= ||g||_{L^2(Gaussian)}^2 + C_X ||grad g||_{L^4(Gaussian)}^2`. (11)

Now let `a` be a fixed smooth even bounded response with bounded
derivatives. Approximate it by even polynomials in Gaussian `W^{1,4}`
(which also controls `L^2`). Apply (11) to the error between `a` and a
fixed polynomial. Covariance Cauchy--Schwarz controls the cross errors in
operator norm, while the polynomial case is (8). Letting the
approximation error tend to zero proves

`|| Cov(S a(X)) - E[a(Z)^2] I ||op -> 0`.                        (12)

The same holds for a fixed family of smooth even responses, with its
Gaussian Gram matrix tensored with `I_n`. If `||B||op<=L` is fixed, their
raw transported outputs consequently have covariance operators
`E[a(Z)b(Z)] Q+o_op(1)`. This is still a covariance statement, not a claim
of a complete locally Gaussian state evolution for arbitrary unmarked
outputs.

## 10. Why the response must be even

The evenness in (8)--(12) cannot be removed. For `a(G)=G`, direct
calculation gives, at distinct roots,
`E[S_i G_i S_j G_j]=B_{ij}^2=1/m`; hence the covariance has a nonvanishing
rank-one all-ones channel.

A more informative example comes from rerooting the four-vertex star.
Let `F_leaf=X_{T3}` be the old degree-three tree (external root is a leaf),
and let `F_center=h3(G)` (root is the center of the same unrooted star).
After omitting its harmless degree-one correction, exact monomial
counting gives

`sum_i S_i F_leaf,i = sqrt(3) sum_i S_i F_center,i`.

The omitted term is `O_{L^2}(n^{-1/2})` after dividing either sum by
`sqrt(n)`, under a fixed operator bound. The limiting diagonal covariance
of these two responses is the identity. The off-diagonal all-ones-channel
coefficient matrix is

`[[2, sqrt(3)], [sqrt(3), 0]]`,

whose eigenvalues are `3` and `-1`. In particular an odd mixture can
have a *negative* rank-one covariance correction: the combination
`F_leaf-sqrt(3)F_center` has an asymptotically vanishing all-ones mode.
The graph reason is that for odd responses the two endpoint parity
degrees are even, so the leading parity graph may be empty. The proof
of Section 8 explicitly excludes that case by its even-response
hypothesis. This obstruction does not contradict the separate weighted
`Q`-trace projection theorem for the unmarked Hermite directions.
