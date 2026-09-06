# Every fixed odd Hermite direction: weighted energy and sign input

Date: 2026-09-05. This generalizes the audited cubic weighted theorem to
each fixed odd degree `r>=3`. Degrees, old-tree families, responses, and
operator bounds are fixed before the matrix order tends to infinity.

Let `B=A/sqrt(m)`, `m=n-1`, be a symmetric hollow signing with
`||B||op<=L`. Put `Q=B^2`, `h_r=He_r/sqrt(r!)`, and

`G=BS`, `Z_r=B h_r(G)`, `v_(r,i)=(B Q^(circ r) B)_ii`.

## 1. A hybrid weighted-polynomial variance estimate

Let the independent input coordinates be arbitrary mixtures of standard
Gaussians and Rademachers. For any fixed real polynomial `f`, any
deterministic real vector `w`, and `T_w=sum_j w_j f(G_j)`,

`Var(T_w) <= C_(f,L) ||w||_2^2`.                              (1)

The product Poincare/Efron--Stein inequality bounds the variance by the
sum of Gaussian-coordinate derivative energies and Rademacher central
difference energies. The ordinary gradient is exactly

`grad T_w = B [w circ f'(G)]`.

Uniform moments of each normalized hybrid sum give

`E||grad T_w||_2^2 <= L^2 sum_j w_j^2 E f'(G_j)^2`
`                                      <= C_(f,L)||w||_2^2`.

For a Rademacher coordinate `k`, compare its central difference with the
ordinary derivative evaluated at the actual input. Taylor expansion of
the fixed polynomial, with `|B_jk|<=1/sqrt(m)`, bounds the difference in
every fixed `L^p` norm by

`C_f sum_j |w_j| B_jk^2 <= C_f ||w||_1/m`.

The sum over `k` of its squared `L^2` bounds is at most

`C_f n ||w||_1^2/m^2 <= C_f n^2 ||w||_2^2/m^2`.

This is uniformly bounded by a fixed multiple of `||w||_2^2` and proves
(1), after the elementary squared-sum inequality. No exact Gaussian
Hermite covariance identity is assumed for a hybrid model.

Finite-degree product hypercontractivity upgrades the centered estimate
to every fixed finite moment:

`||T_w-E T_w||_p <= C_(f,L,p)||w||_2`.                         (2)

For odd `f`, each `E f(G_j)=0` exactly. For an even normalized Hermite
`h_s`, `s>=2`, the scalar fourth-order replacement bound gives

`max_j |E h_s(G_j)| <= C_s/m`.                               (3)

Indeed each normalized row has variance one, the coordinate laws match
the first three Gaussian moments, and the sum of fourth powers of its
coefficients is `1/m`. Polynomial Taylor remainders have bounded moments
uniformly in the replacement stages.

These observations imply two estimates needed below:

`||sum_j B_ij h_s(G_j)||_p <= C_(s,L,p)` for odd `s>=1`,

`||sum_j B_ij B_jk h_s(G_j)||_p <= C_(s,L,p)/sqrt(n)`
`                                                   for even s>=2`. (4)

Use `||b_i||_2=1` in the first, and
`||b_i circ b_k||_2<=1/sqrt(m)`, `||b_i circ b_k||_1<=1` in the second.

## 2. Exact derivatives of every odd unmarked transport

For an ordinary input derivative `partial_k` of order `a<=r`,

`partial_k^a Z_(r,i)`
` = sqrt(r!/(r-a)!) sum_j B_ij B_jk^a h_(r-a)(G_j)`.          (5)

When `2<=a<r` is even, flatness writes (5), apart from its fixed constant,
as

`m^(-a/2) [B h_(r-a)(G)-B_ik h_(r-a)(G_k)]_i`.

The remaining Hermite degree is odd, so (4) bounds its pointwise fixed
moments by `O(n^(-a/2))`. When `a<r` is odd, (5) is instead

`m^(-(a-1)/2) sum_j B_ij B_jk h_(r-a)(G_j)`.

The remaining even degree is at least two, and the second part of (4)
gives exactly the same pointwise bound. The fully differentiated case
is the deterministic vector

`partial_k^r Z_(r,i)=sqrt(r!) m^(-(r-1)/2) Q_ik`,

whose Euclidean norm is `O_L(n^((1-r)/2))`. Derivatives above `r` vanish.

Consequently, for every `0<=a<=4`,

`|| ||partial_k^a Z_r||_2 ||_p <= C_(r,L,p) n^((1-a)/2)`.     (6)

At `a=0`, use the first estimate in (4). For `r=3,a=3`, use the exact
vector estimate instead of the false pointwise `n^-3/2` assertion.
For `r>=5`, every derivative through order four is covered pointwise.

Freezing the replacement coordinate at `t` adds only a fixed polynomial
factor in `1+|t|`. For precision, first make that coordinate Gaussian,
expand the fixed-degree polynomial in its one-dimensional Hermite basis,
and use coefficient projection bounds and (6). Evaluation of the finite
Hermite expansion at `t` gives the assertion uniformly in all other
hybrid coordinates. This avoids inferring polynomial coefficients merely
from their values at the two Rademacher endpoints.

## 3. Whole-functional Rademacher replacement

Let `M,F` be fixed old-tree polynomial responses and let `d_i` be bounded
deterministic numbers. Consider

`T(s)=n^-1 [d circ M(X(s)) circ Z_r(s)]^T B F(X(s))`.

The injective old fields are affine in each input coordinate and have
uniform coordinate derivative moments `O(n^-1/2)`. Chain rule supplies
the vector scale `n^((1-a)/2)` for derivatives of each response through
order four. Products with `Z_r` preserve this scale by (5)--(6). For the
fully differentiated deterministic vector, multiply its entries by the
pointwise response derivative moments, exactly as in the cubic proof.

Thus the fourth derivative of the complete normalized functional has
`L^1` norm `O(n^-2)`, with the corresponding polynomial-in-`t` bound on
every Taylor segment. The input laws match three moments. Replacing all
`n` coordinates proves

`E T(S)=E T(N)+O(n^-1)`.                                     (7)

The constants can depend on the fixed degree `r`, but not on the input
hybrid stage, root labels, deterministic weights within their bound, or
matrix order.

## 4. Complete Gaussian weighted projection at odd degree r

The exact result is

`n^-1 E[d circ M(X) circ Z_r]^T B F(X)`

` = E M(gamma) * E[F(gamma)h_r(gamma_edge)]`
`                               * n^-1 sum_i d_i v_(r,i) + o(1)`. (8)

Here `gamma` is the independent old-tree Gaussian family. We give every
change needed from `fresh_full_weighted_unmarked_projection_2026_09_05.md`.

Gaussian integration by parts now differentiates
`M(Y_i)F(Y_k)` exactly `r` times in direction `b_j`, with prefactor
`1/sqrt(r!)`. The one-level fields `Y` and all their fixed-order
directional derivative matrices have the same polylogarithmic operator
bounds proved for arbitrary order in the matrix-derivative module.

Every mixed allocation `a>=1,b>=1,a+b=r` is exactly a Frobenius pairing

`n^-1 <B circ D^a(d M), B D^b F>_F`

and is `o(1)` by flatness and those operator bounds.

In the one-response `r`th chain rule, retain its all-edge term and its
single nonedge `r`th derivative. Every other term has Frobenius moments
`O(polylog(n))`. If at least two differentiated factors are nonedge,
their pointwise product is `O(n^-1)` and the squared entry sum is
bounded. If exactly one factor is nonedge but another factor is present,
that other factor is an edge derivative and supplies at least one
`Q_ij`; squared entry summation then uses `sum_ij Q_ij^2=O(n)`.
Higher edge derivatives vanish. This exhausts the chain partitions.

The sole `r`th derivative of a nonedge one-level tree has exact form

`D^r Y_T=B diag(N)K_r+rB(B circ K_(r-1))`.

The root-hit term is small in Frobenius norm. In the root-not-hit term,
ordinary-coordinate integration by parts differentiates the local
diagonal response factors; the derivative of row `l` of `K_r` is zero
by exact own-coordinate exclusion. The all-on-F trace is controlled by
polylogarithmic operator norms and the `n^-1/2` coordinate influence.
The all-on-M trace is controlled by the coordinate-Jacobian Frobenius
bound for `M_T(Y_i)(BF)_i`, exactly as in Section 6 of the cubic assembly.
Neither argument uses the value `r=3`.

For the all-on-M all-edge term, `r>=3` gives

`|sum_j B_ij Q_ij^r| <= m^-1/2 sum_j Q_ij^2=O(n^-1/2)`.

It vanishes. The all-on-F all-edge term has weights

`P_ik=d_i B_ik (B Q^(circ r))_ik`.

The two-root Gaussian polynomial moment theorem factors
`E[M(Y_i)F_(g...g)(Y_k)]` into its two means, with error
`O(|Q_ik|)+o(1)`. The normalized absolute weight sum is bounded, and
the additional `|Q_ik|` sum is `O(n^-1/2)`, since
`||Q^(circ r)||op<=L^2`. Finally the row weight sum is `d_i v_(r,i)`
and Gaussian integration by parts identifies the coefficient as
`E[F h_r(G0)]`. This proves (8) for Gaussian input; (7) proves it for
Rademacher input.

## 5. Local Gaussian independence for the higher odd transports

For completeness, the fixed-degree local joint limit needed for bounded
cutoffs also extends. With Gaussian input, the normalized kernel of
`Z_(r,i)` is

`K_(r,i)=sum_j B_ij b_j^(tensor r)`.

The Gaussian field is `I_r(K_(r,i))/sqrt(r!)`, so its variance is
`||K_(r,i)||_HS^2=v_(r,i)`. Any proper flattening splits the `r` slots
into nonempty sets of sizes `a,b` and has the form

`V_a^T diag(b_i) V_b`,

where the row Gram matrices of `V_a,V_b` are `Q^(circ a)` and
`Q^(circ b)`. Its operator norm is at most `L^2/sqrt(m)`. Hence every
nontrivial self-contraction vanishes uniformly, and the fixed-chaos
Gaussian central limit theorem applies, including subsequences on which
the variance tends to zero.

Different-degree old fields are orthogonal. For an old tree of degree
`r`, use the exact contraction recursion from Section 7 of
`fresh_limit_unmarked_cubic_projection_2026_09_05.md`. First ignore
collisions in the old-tree kernel, and form its columns
`phi_T^j(i)=<K_T,i,b_j^(tensor r)>` obey

`phi_T^j=B diag(b_j) product_children phi_child^j`.

For a nonstar their column norms are `O(n^-1)`, hence the matrix has
Frobenius norm `O(n^-1/2)`; its row scalar product with `b_i` vanishes
uniformly. For a star the matrix is a fixed multiple of
`B[B circ Q^(circ(r-1))]`. The middle matrix has absolute row sums
`O(n^-1/2)` because `r-1>=2`, so its diagonal product with `B` also
vanishes. Restoring old-tree injectivity costs `O(n^-1/2)` by pairing
the fixed-root collision kernel directly with the already summed kernel
`K_(r,i)`, whose Hilbert norm is bounded. One must not add separate
per-column collision errors over `j`. These are exactly the covariances
of `Z_(r,i)` with the old degree-`r` fields. The old marginal contractions
already vanish, so the
multivariate chaos theorem gives independence from the entire finite old
family.

Here is a uniform collision check for transfer to sign input. A prescribed
pair of repeated kernel slots has coefficients

`K_(r,i)(a,a,c_3,...,c_r)`
` = m^-1 [K_(r-2,i)(c_3,...,c_r)`
`                          -B_ia product_(s>=3) B_(a,c_s)]`.

The first term has squared norm at most `n L^2/m^2`; the second at most
`1/m^2`, because each row of `B` is a unit vector. Thus deleting all
repeated slots changes the fixed-root kernel in Hilbert norm by
`O(n^-1/2)`. On sign input, the scalar Hermite-versus-multilinear
recurrence has row error `O_Lp(1/m)` for each fixed degree. Transport by
`B` changes any individual root in `L2` by at most
`L [sum_j O(m^-2)]^(1/2)=O(n^-1/2)`. The remaining multilinear kernels
have influences `O(1/n)` by their proper flattenings. Fixed-degree
invariance therefore proves the uniform joint local limit

`(X_i,Z_(r,i)) -> (gamma,sqrt(v_(r,i)) N0)`,

where `N0` is independent of `gamma`. Uniform fixed moments give the
required uniform integrability. This also justifies extension of (8) to
fixed bounded a.e.-continuous responses by the same polynomial
approximation argument as in the cubic theorem.

## 6. Generic strict escape from scalar Gaussian certificates

For every odd `r`, the spectral inequality proved separately is

`n^-1 sum_i sqrt(v_(r,i))>=1`.

Its Schur map is `Psi(X)=Q^(circ((r-1)/2)) circ X`, with intermediate
matrix `P=Q^(circ((r+1)/2))`; square-root Jensen and trace duality give
`tr[Q sqrt(Q^(circ r))]>=tr(P^(3/2))>=n`.

Thus the variance normalization, unit Gaussian cutoff, localized slack,
`c/log L` gain, and principal-deletion transfer of
`fresh_normalized_unmarked_gain_and_limit_order_2026_09_05.md` apply
unchanged whenever the target final response has a nonzero fixed higher
odd edge coefficient. All constants can depend on that chosen finite
degree; they remain independent of `L` in the required covariance
neighborhood.

A useful sufficient condition needs no search over individual Hermite
coefficients. Let `F` be a bounded odd Gaussian response and put
`f(g)=E[F|G0=g]`. If `E[G0 F]!=0`, then some odd coefficient
`E[F h_r(G0)]` with `r>=3` is nonzero. Otherwise Hermite completeness
would imply `f(g)=a g` almost everywhere, with `a=E[G0 F]!=0`; this
contradicts boundedness of the conditional expectation `f`.
More quantitatively, if `|F|<=1`, its nonlinear conditional residual
satisfies

`sum_(r>=3, odd) |E[F h_r(G0)]|^2`
` >= E (|a G0|-1)_+^2 >0`.

Consequently any scalar-mask certificate satisfying the nondegenerate
local slack hypotheses and `E[G0 F]!=0` admits a strict unrestricted
unmarked improvement. This is an obstruction to optimality of that
certificate family. It does not prove that the enlarged variational
family exhausts the Boolean optimum, that it has a uniform improvement
over its own supremum, or that the sequence `M_n/n^(3/2)` converges.

Independent audit record: literature checked Sections 1--3, including
the hybrid endpoint bound and the frozen-coordinate Hermite projection;
algebra checked Section 5, including summed-kernel collision control and
Gaussian normalization. The Gaussian energy proof uses the already
independently audited derivative matrix estimates at arbitrary fixed
order. A separate compactness argument in
`fresh_uniform_scalar_hierarchy_escape_2026_09_05.md` does upgrade this
pointwise escape to a uniform gap above the entire scalar central-mask
subfamily; it makes no such claim about the enlarged unmarked family.
