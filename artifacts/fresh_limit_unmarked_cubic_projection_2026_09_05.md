# A bounded-operator-norm unmarked cubic module

Date: 2026-09-05. This is a new finite-field result, not a general AMP assertion.
All limits below are along actual symmetric hollow sign matrices. No entrywise
coherence assumption on powers of the matrix is made.

Primary Gaussian-chaos reference: Noreddine--Nourdin, *On the Gaussian
approximation of vector-valued multiple integrals*,
<https://arxiv.org/pdf/1009.1310>, Theorems 1.1 and 1.3 and Proposition 1.4.
The fixed-dimensional vector here consists of one multiple integral per
tree, together with the cubic integral; their covariance limits and
vanishing self-contractions verify the fourth-cumulant hypotheses.
Theorem 1.3 permits a singular limiting covariance, including a root whose
cubic variance tends to zero. Degree-one Gaussian coordinates can be
adjoined using the same derivative-covariance bound; their nontrivial
self-contractions are absent. The input replacement is separately proved
below, not attributed to that Gaussian theorem.

Independent audits: `fresh_unmarked_projection_literature_audit_2026_09_05.md`
and the variational agent's direct Section 1--5 audit both checked the
normalizations, arbitrary-cut flattening, exceptional `BCQ` contraction,
and weighted replacement. Sections 7--8 were subsequently checked by the
director; a further independent audit was requested separately.

## 1. Definitions and conclusions

Let `m=n-1`, `B=A/sqrt(m)`, `Q=B^2`, and suppose `||B||op <= L`, with `L`
fixed. Let `S` have independent Rademacher coordinates. Write `G=BS` and
`h3(t)=He_3(t)/sqrt(6)`. Let `X_i=(X_{T,i})` be any fixed finite family of
the fully marked, odd-degree, own-root-free tree fields defined in
`fresh_limit_hierarchical_tree_energy_2026_09_05.md`, including the single
edge `X_edge,i=G_i`. Its limiting vector `X` consists of independent standard
Gaussians, with edge coordinate denoted by `G`.

Define the genuinely unmarked field `Z=B h3(G)`, coordinatewise. Then:

1. Uniformly in the output root, `(S_i,X_i,Z_i)` is asymptotically distributed
   as `(S,X,sqrt(v_i) N)`, where `S`, `X`, and `N` are independent and
   `v_i=(B Q^{circ 3} B)_{ii}`. This is Gaussian approximation with a possibly
   root-dependent variance, not a claim that the variances converge.
2. The variances satisfy

   `0 <= v_i <= L^2`,

   `n^{-1} sum_i v_i = n^{-1} sum_{j,k} Q_{jk}^4 >= 1`.                 (1)

3. For every fixed polynomial `F` in the old tree coordinates,

   `n^{-1} E F(X_i)^T Q h3(G)`

   `= E[F(X)h3(G)] * n^{-1} sum_{i,j} Q_{ij}^4 + o(1)`.             (2)

   Here the vector in the left side has coordinates `F(X_i)`. The error is
   uniform over all the indicated matrices with the fixed bound `L`.
   The same formula holds for bounded continuous `F`, and for bounded
   Gaussian-almost-everywhere continuous `F` by the usual approximation.

Equation (2) is the precise positive-direction bilinear projection: no deeper
old-tree coordinate survives except through the Gaussian conditional
projection of `F` onto `G`. In particular its sign is the sign of
`E[F(X)h3(G)]`, up to an error tending to zero.

These statements alone are not yet a Boolean energy improvement. A final
response depending on `Z_i` creates weighted correlations with `B F(X)`;
those correlations need a separate transport/endpoint argument.

## 2. Proper flattening bound for a marked tree

Use the Gaussian tensor convention `I_d(K)` with variance `d! ||K||^2`.
Every fixed tree kernel has bounded Hilbert norm. Every proper flattening of
a non-edge marked tree kernel has operator norm `O_{T,L}(n^{-1/2})`.

Here is a direct proof, including arbitrary cuts. Initially allow collisions
among all free vertex labels. The tensor entry is the product of the `B`
entries on the tree edges, with the output root fixed. Divide the free
vertices into two nonempty parts. Crossing edges form a bipartite graph.
For vertices incident to crossing edges, copying a vertex label into its
incident edge slots is an isometry. Thus the crossing-edge matrix is a
compression of the tensor product of the crossing copies of `B`, and has
operator norm at most `L^{E_cross}`. Each vertex isolated in the crossing
graph instead contributes one summation/constant-vector factor `sqrt(n)`.
The internal free-free edges multiply rows or columns by signs of modulus
at most `m^{-1/2}`. The edge incident to the fixed external root supplies
one more such diagonal factor. The resulting bound is

`L^{E_cross} n^{I/2} m^{-(E_internal+1)/2}`.

Since the free-vertex graph is a connected tree and the cut is nonempty,
every component of its internal-edge forest contains a crossing-incident
vertex. In each such component the number of crossing-isolated vertices
is at most its number of edges. Consequently `I <= E_internal`, proving
the required bound. Symmetrizing a fixed number of slots does not change
the order of this bound.

Finally delete collisions between free labels or with the fixed root.
There are `O_T(n^{d-1})` forbidden tuples, each coefficient has squared
magnitude `O_T(m^{-d})`, and the kernel change therefore has Hilbert norm
`O_T(n^{-1/2})`. Every flattening operator norm is bounded by that norm.
This proves the statement for the actual injective tree kernels.

It follows, by Gaussian hypercontractivity, that for any unit vector `b`
and every non-edge tree of degree `d`,

`||D_b^r X_{T,i}||_p = O_{T,L,p}(n^{-1/2})`, `1 <= r <= min(3,d)`. (3)

For `r=d`, use any proper flattening to bound the full scalar contraction;
for `r<d`, flatten at the `r` contracted slots. In contrast,
`D_{b_j}G_i=Q_{ij}` and its higher derivatives vanish.

## 3. The unmarked cubic central limit theorem

Let `b_j` be row `j` of `B` and put

`K_i = sum_j B_{ij} b_j^{tensor 3}`.

For Gaussian input, the cubic Wick field is `I_3(K_i)/sqrt(6)`, with
variance `||K_i||^2=(B Q^{circ 3} B)_{ii}`. Put
`D_i=diag(B_{i1},...,B_{in})` and `R=Q^{circ 2}`. The two nontrivial
contractions have squared norms

`||K_i tensor_1 K_i||^2 = Tr[(D_i Q D_i) R (D_i Q D_i) R]`,

`||K_i tensor_2 K_i||^2 = Tr[(D_i R D_i) Q (D_i R D_i) Q]`.          (4)

All Schur powers of the correlation matrix `Q` have operator norm at most
`||Q||op <= L^2`: Schur multiplication by a correlation matrix is a unital
positive map and is operator-norm contractive on self-adjoint matrices.
Also `||D_i||op <= m^{-1/2}` and `||Q||F^2 <= L^2 n`. Thus each expression
in (4) is `O_L(n^{-1})`.

For each input coordinate `k`,

`|| e_k contracted K_i ||^2`

`= sum_{j,l} B_{ij}B_{jk} B_{il}B_{lk} Q_{jl}^2`

`<= ||R||op sum_j (B_{ij}B_{jk})^2 = O_L(n^{-1})`.                (5)

This controls every input influence, including the own input `k=i`.
The third-chaos central limit theorem and its multivariate version now
apply, including subsequences on which `v_i` converges. The contraction
bounds also give uniform approximation without selecting such a
subsequence. Formula (1) follows from cyclicity of trace and `Q_{ii}=1`.

For completeness, the replacement of repeated-coordinate Gaussian Wick
kernels by the Rademacher multilinear polynomial costs `o(1)` uniformly.
The row identity is

`He_3(G_j) = sum_{a,b,c distinct} B_{ja}B_{jb}B_{jc} S_a S_b S_c
             - (2/m)G_j`.

Transport of the last term is `-(2/m)Q S`, of uniformly vanishing `L^2`
 norm. In the Gaussian cubic kernel, repeated slots have squared norm
 `O_L(1/n)`: explicitly

 `K_i(a,a,b) = m^{-1}(Q_{ib}-B_{ia}B_{ab})`,

 so summing its square over `a,b` is at most
 `2m^{-2}(n (Q^2)_{ii}+1)=O_L(1/n)`. Union over the three possible
 coincident slot pairs gives the assertion. Third-chaos invariance for the
remaining multilinear kernels follows from (5).

All non-edge old trees have odd degrees at least three. Their Gaussian
chaoses are orthogonal to the cubic unless their degree is three. There is
only one degree-three tree. Before collision deletion its kernel is

`(1/sqrt(2)) Sym sum_k B_{ik} e_k tensor b_k tensor b_k`.

Set `C=B circ Q^{circ 2}`. Its cross covariance with `h3(G_j)` is
`sqrt(3) (BC)_{ij}`. After transport by `B`, covariance with `Z_i` is
`sqrt(3) (BCB)_{ii}+o(1)`. Since

`max_j sum_k |C_{jk}| <= m^{-1/2} sum_k Q_{jk}^2 <= L^2/sqrt(m)`,

we have `||C||op=O_L(n^{-1/2})`; hence that covariance vanishes uniformly.
Collision deletion adds `o(1)` by Hilbert-norm Cauchy--Schwarz. All the old
tree self-contractions vanish by Section 2, and their mutual Gaussian
limit is the already proved independent tree family. The multivariate
chaos theorem proves joint independence of `Z_i` from that entire family.
Removing its own-spin dependence using (5) proves independence from `S_i`.

## 4. Gaussian-input proof of the weighted projection

All fields in this section are the injective Gaussian Wick tree fields.
Gaussian integration by parts gives

`E F(X_i) h3(G_j) = (1/sqrt(6)) E D_{b_j}^3 F(X_i)`.              (6)

Write `f_a,f_ab,f_abc` for the partial derivatives of the fixed polynomial
`F`, evaluated at `X_i`. The chain rule separates the right side into

`sum_{a,b,c} f_abc (D X_a)(D X_b)(D X_c)`;

`3 sum_{a,b} f_ab (D^2 X_a)(D X_b)`;

`sum_a f_a D^3 X_a`,                                            (7)

where `D=D_{b_j}`. Hypercontractivity bounds all fixed moments of the
derivative factors `f`. We weight by `Q_{ij}`, sum `i,j`, and divide by
`n`. The useful uniform sums are

`sum |Q_{ij}| <= n^{3/2}L`, `sum Q_{ij}^2 <= L^2 n`,

`sum |Q_{ij}|^r <= L^2 n` for every fixed `r>=2`.

In the first line of (7), if exactly one differentiated coordinate is a
non-edge tree, the other two supply `Q_{ij}^2`, and (3) makes its normalized
contribution `O(n^{-1/2})`. If exactly two are non-edge trees, their product
is `O_{L^p}(n^{-1})` and the remaining edge supplies another `Q_{ij}`;
the contribution is `O(n^{-1})`. If all three are non-edge trees, their
product is `O(n^{-3/2})` and the absolute-sum bound makes the contribution
`O(n^{-1})`. The all-edge term is retained.

In the second line the coordinate with two derivatives must be non-edge.
If the other one is edge, (3) and the extra `Q_{ij}` give `O(n^{-1/2})`.
If both are non-edge, their product is `O(n^{-1})`, giving `O(n^{-1/2})`
after the absolute-sum bound.

The third line must be summed over `j` before it is bounded. Define

`S_i^{(3)} = sum_j Q_{ij} b_j^{tensor 3}`.

Its squared Hilbert norm is
`Q_i^T Q^{circ 3} Q_i <= L^2 ||Q_i||^2 <= L^6`. For a tree of degree
greater than three, flattening its kernel at three slots therefore gives

`|| sum_j Q_{ij} D_{b_j}^3 X_{T,i} ||_p = O(n^{-1/2})`.

For the degree-three tree, the collision-free leading expression from
Section 3 instead gives a fixed multiple of `(BCQ)_{ii}`, whose absolute
value is `O_L(n^{-1/2})`; collision deletion has the same error by the
bound on `S_i^{(3)}`. Hence the whole third line contributes `o(1)`.

Only the all-edge term remains:

`(1/sqrt(6)n) sum_{i,j} Q_{ij}^4 E F_{GGG}(X_i)`.

Uniform rooted moment convergence replaces the expectation by
`E F_{GGG}(X)`. Gaussian integration by parts in the limiting vector
identifies `E F_{GGG}(X)/sqrt(6)=E[F(X)h3(G)]`, proving (2) in this model.

## 5. Weighted Rademacher replacement and extension

Use the same injective multilinear kernels for hybrid independent inputs,
each input being either a standard Gaussian or a Rademacher variable.
For every fixed tree and every input index `k`, direct coefficient
counting and fixed-degree hypercontractivity give

`||D_k X_{T,i}||_p <= C/sqrt(n)`.

The field is exactly affine in each input. Consequently the fourth
derivative in that input of the product `F(X_i) h3(G_j)` has bounded
`L^1` norm `C/n^2`, uniformly in the hybrid replacement stage. This follows
by the ordinary chain/product rule: every one of the four derivatives
contributes either one `D_k X` or one `B_{jk}`, and all fixed moments of
the undifferentiated polynomial factors are bounded. Evaluation along the
Taylor segment adds only fixed polynomial moments of the replaced input.
Gaussian and Rademacher variables match their first three moments.
Fourth-order Taylor replacement therefore bounds the difference of each
pair expectation by `C/n`. Summing with weights `Q_{ij}/n` costs at most
`L sqrt(n)`, so the complete weighted difference is `O_L(n^{-1/2})`.

Finally, for a bounded continuous `F`, approximate it in the limiting
Gaussian `L^2` space by a fixed polynomial `P`. Rooted moment convergence
and uniform integrability give the corresponding averaged matrix-input
`L^2` approximation. On either input model,

`n^{-1} |E (F-P)^T Q h3(G)|`

`<= L^2 [n^{-1}E||F-P||^2]^{1/2}
           [n^{-1}E||h3(G)||^2]^{1/2}`.

The second factor stays bounded and tends to one. Apply the polynomial
result first, then let the approximation error tend to zero. This proves
the stated extension, with the same argument for bounded functions whose
discontinuity set has zero limiting Gaussian measure.

## 6. Scope and unresolved energy step

The theorem is valid for an arbitrary fixed old tree family and a fixed
bounded operator norm. Spectral diagonal regularization already reduces
universal Boolean lower bounds to this regime, but that reduction does
not itself complete the next step.

Equation (2) yields the exact positive direction

`n^{-1} E [B F(X)]^T [B h3(G)]
    = b_3 n^{-1}sum Q_{ij}^4 + o(1)`,

where `b_3=E[F(X)h3(G)]`. It gives a nonnegative, nonzero correlation when
`b_3>0`. Turning this into an admissible bounded response requires control
of the local weights introduced by a perturbation of the current Boolean
decision rule. We have not asserted those weighted identities here.

## 7. Extension of the projection to every fixed Hermite degree

The weighted projection has a stronger conclusion than the cubic local
central limit theorem. For every fixed integer `r>=2`, put
`h_r=He_r/sqrt(r!)`. Then, for every fixed old-tree polynomial `F`,

`n^{-1} E F(X)^T Q h_r(G)`

`= E[F(X)h_r(G)] * n^{-1}sum_{i,j} Q_{ij}^{r+1} + o(1)`.          (8)

The proof in Sections 4--5 extends as follows. In the `r`th Gaussian chain
rule, a term with at least two differentiated non-edge fields has a factor
`O(n^{-1})`; after the absolute weighted sum its contribution vanishes.
A term with exactly one differentiated non-edge field and at least one
edge derivative retains a factor `Q_{ij}` in addition to the outer weight,
so the square-sum bound makes its contribution vanish. The all-edge term
is precisely the right side of (8). The sole remaining case is

`sum_T E F_T(X_i) sum_j Q_{ij} D_{b_j}^r X_{T,i}`.                 (9)

If `d(T)<r` it vanishes. If `d(T)>r`, flatten at `r` slots against
`sum_j Q_{ij} b_j^{tensor r}`, whose squared norm is at most `L^6`, exactly
as before. Hence (9) is `O(n^{-1/2})`. If `d(T)=r`, only odd `r` can occur.
The following explicit recursive contraction closes this final case.

Ignore the fixed, tree-dependent normalization constants temporarily and
allow collisions among the free labels. Write

`phi_T^j(i)=<K_{T,i},b_j^{tensor d(T)}>`.

For the edge, `phi_edge^j=Q e_j`. If the top vertex of `T` has child trees
`T_1,...,T_k`, where `k>=2` is even, then, up to a fixed normalization
factor,

`phi_T^j = B diag(b_j) [ product_{a=1}^k phi_{T_a}^j ]`.         (10)

The product is coordinatewise. Inductively, every fixed non-edge tree has
`||phi_T^j||_2=O(n^{-1/2})`, uniformly in `j`. For a star, this follows from
`||B diag(b_j)||op <= L/sqrt(m)` and
`||(Q e_j)^k||_2 <= ||Q e_j||_2 <= L`.
For a non-star, at least one child is non-edge. Its `l2` norm is
`O(n^{-1/2})`, all other child `linfinity` norms are bounded, and (10)
therefore improves to

`||phi_T^j||_2=O(n^{-1})` for each fixed non-star tree.            (11)

Consequently the matrix with columns `phi_T^j` has Frobenius norm
`O(n^{-1/2})` for a non-star. Each of its rows has the same upper bound,
and its scalar product with a row of `Q` is `O(n^{-1/2})`.

For a star with `r-1` leaves, the corresponding matrix is a fixed multiple
of `B [B circ Q^{circ(r-1)}]`. Since `r-1>=2` is even and `|Q_{ij}|<=1`,
the middle matrix has absolute row sums at most
`m^{-1/2} sum_j Q_{ij}^2 <= L^2/sqrt(m)`. Thus its weighted diagonal
product with `Q` is also `O(n^{-1/2})`. Finally collision deletion changes
the contraction with `sum_j Q_{ij} b_j^{tensor r}` by `O(n^{-1/2})` in
view of its bounded Hilbert norm. This proves (9) is negligible in every
case and completes the Gaussian proof of (8).

The fourth-order replacement from Section 5 applies without change to any
fixed polynomial `h_r`: its fourth and lower derivatives have bounded
fixed moments, and every input derivative still supplies a flat column
entry. Bounded a.e.-continuous `F` is covered by the same `L^2` extension.

## 8. A positive direction for a general old-hierarchy response

Let `F` be a bounded, smooth, odd function of the fixed old Gaussian tree
family, and set

`f(g)=E[F(X)|G=g]`, `a=E[G f(G)]`,

`r_0^2=E f(G)^2-a^2`.

Suppose `r_0>0`. Let `h(g)=(f(g)-ag)/r_0`, and assume, for the convenient
matrix-input extension below, that the positive-order derivatives of `F`
are bounded. Then `h` is odd, has at most linear growth, has bounded
positive-order derivatives, and is orthogonal to the constant and linear
Gaussian functions. If `f_s=E[f(G)h_s(G)]`, its Hermite coefficients are
`h_s=f_s/r_0` for odd `s>=3` and zero at degree one.

Applying (8) term by term and then taking an `L^2` limit gives

`n^{-1} E F(X)^T Q h(G)`

`= sum_{s>=3, s odd} (f_s^2/r_0)
           [n^{-1}sum_{i,j}Q_{ij}^{s+1}] + o(1) >= r_0-o(1)`.   (12)

The series is uniformly bounded and converges uniformly with respect to
the matrices in its tail: each bracket lies between `1` and `L^2`, since
`s+1` is even and `sum Q_{ij}^2<=L^2n`. The Gaussian and sign-input
bilinear truncation errors are controlled by `||Q||op<=L^2` and the
limiting one-row Gaussian `L^2` approximation to `h`; all rows have the
same normalized Rademacher-sum marginal distribution. Smoothness and
linear growth provide the needed moment bounds.

Thus the exact Hermite-positive residual construction previously proved
only for `F(G)` also works for an arbitrary finite old-hierarchy response,
provided its conditional edge projection has a nonzero nonlinear part.
Combined with the influence, truncation, orientation-average, and exact
own-spin endpoint argument in
`fresh_bounded_op_positive_first_variation_2026_09_05.md`, this gives a
strict finite-step gain over each fixed smooth slack baseline at each
fixed operator bound. The size of that gain can depend on the bound, the
response, and the available slack. This statement does not claim that it
already dominates a vanishing smoothing error or a spectral-deletion
loss uniformly as those auxiliary parameters change.
