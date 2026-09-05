# Complete weighted cubic projection for a bounded-operator signing

Date: 2026-09-05. This assembles the Gaussian derivative classes and the
independently audited whole-functional input comparison. The theorem is
for a fixed finite old marked-tree construction and a fixed operator
bound. No parameter below grows with the matrix order.

## 1. Statement

Let `B=A/sqrt(n-1)` be a symmetric hollow signing with `||B||op<=L`,
where `L` is fixed. Write `Q=B^2`, `G=BS`, and

`Z=B h3(G)`, with `h3(g)=(g^3-3g)/sqrt(6)`.

Let `X_i` be a fixed finite family of normalized injective old marked-tree
fields, including its edge coordinate `G_i`. Let `M,F` be fixed polynomial
responses. Put

`vbar_n = n^-1 sum_ij Q_ij^4`.

Then

`n^-1 E [M(X) circ Z]^T B F(X)`

` = E M(gamma) * E[F(gamma) h3(gamma_edge)] * vbar_n + o(1)`. (1)

Here `gamma` is the independent standard Gaussian old-tree family. In
particular `1<=vbar_n<=L^2`, since `Q_ii=1`, `|Q_ij|<=1`, and
`sum_ij Q_ij^2<=L^2 n`.

The proof does not require `M` even or `F` odd. Those symmetries are useful
when applying (1) to cube-feasible paired responses, not for the identity.

There is also a deterministic weighted version. If `d_i` are arbitrary
deterministic numbers with `sup_i|d_i|<=D` for a fixed `D`, then

`n^-1 E [d circ M(X) circ Z]^T B F(X)`

` = E M(gamma) * E[F(gamma) h3(gamma_edge)]`
`                         * n^-1 sum_i d_i v_i + o(1)`,       (1a)

where `v_i=(B Q^(circ3) B)_ii`. The error is uniform over this bounded
choice of deterministic weights. In the proof every `M(Y_i)` is simply
replaced by `d_i M(Y_i)`; no derivative hits `d_i`. Every norm estimate
changes by at most a factor depending on `D`, and the main term has
`P_ik=d_i B_ik(BR)_ik`. Summing `k` gives `d_i v_i` exactly.

## 2. Audited inputs and notation

The derivative estimates and the one-level representation below are
proved in `fresh_tree_matrix_derivative_and_weighted_third_chain_2026_09_05.md`.
They rely on the primary matrix-chaos theorem quoted there, with every
coefficient flattening checked for the exact injective kernels.

Work first with independent standard Gaussian input `N`. For each nonedge
tree, replace its exact field by

`Y_T,i = sum_l B_il N_l h_T(X_l)`; set `Y_edge=G`.

Uniformly in the root, `||Y_T,i-X_T,i||_p=O(n^-1/2)` for every fixed finite
`p`. Bounded operator norm and fixed polynomial moments allow this
replacement in the whole functional in (1), at cost `o(1)`.

Write `D_j=D_{b_j}` for differentiation in the unit row direction `b_j`.
For every fixed polynomial response `a`, every fixed derivative order
`r>=1`, and every fixed finite `p`, the matrices

`A_r(i,j)=D_j^r a(Y_i)`

have operator-norm `L^p` bounds by a fixed power of `log(n+2)`. This also
holds for the one-level fields: their exact derivative matrices are

`D^r Y_T = B diag(N) K_r + r B(B circ K_(r-1))`,               (2)

where `K_r(l,j)=D_j^r h_T(X_l)`. When `r=1`, the root-hit term in (2) is
`B diag(h_T(X)) B`; for higher `r` use the Hadamard-product operator
bound. Chain rule then gives the asserted response bounds.

For a nonedge tree, all fixed row-direction derivatives have pointwise
`L^p` size `O(n^-1/2)`. For the edge, `D_j G_i=Q_ij` and its higher
derivatives vanish. Every ordinary coordinate derivative of a local
polynomial response has pointwise size `O(n^-1/2)`.

We abbreviate a bound by a fixed power of `log(n+2)` as `polylog(n)`.
All bounds in this proof also hold in every fixed higher moment needed
by Holder. No random maxima are discarded without such a moment bound.

## 3. Gaussian integration by parts and every mixed split

Gaussian Hermite integration by parts gives the exact identity

`n^-1 E [M(Y) circ Z]^T B F(Y)`

` = (sqrt(6)n)^-1 sum_ijk B_ij B_ik E D_j^3[M(Y_i)F(Y_k)]`.   (3)

If `a>=1`, `b>=1`, and `a+b=3`, set

`A_a(i,j)=D_j^a M(Y_i)`, `C_b(k,j)=D_j^b F(Y_k)`.

The corresponding sum before its fixed binomial coefficient is exactly

`n^-1 E <B circ A_a, B C_b>_F`.

Flatness and the derivative operator bounds imply

`||B circ A_a||_F <= ||A_a||_F/sqrt(n-1) = O_Lp(polylog(n))`,

`||B C_b||_F = O_Lp(sqrt(n) polylog(n))`.

Both mixed derivative allocations therefore contribute
`O(polylog(n)/sqrt(n))=o(1)`. This treats their entire chain rules at once.

## 4. The third chain rule: a useful remainder bound

For any polynomial response `a`, its third row-direction derivative has
the matrix decomposition

`[D_j^3 a(Y_i)] = diag(a_ggg(Y)) Q^(circ3) + R_a`

`                         + sum_(T nonedge) diag(a_T(Y)) [D_j^3 Y_T,i]`. (4)

The remainder `R_a` has Frobenius `L^p` norm `O(polylog(n))`.
Here is the count behind that assertion. In the product of three first
derivatives, a term with one nonedge factor and two edge factors has
squared entry bounds `C n^-1 Q_ij^4`; summing gives `O(1)`. A term with
two nonedge factors and one edge factor has squared entry bounds
`C n^-2 Q_ij^2`, and a term with three nonedge factors has squared entry
bounds `C n^-3`. Both sums are bounded. In a second-derivative times a
first-derivative term, the twice-differentiated tree is nonedge. An edge
first derivative gives squared bounds `C n^-1 Q_ij^2`; a nonedge first
derivative gives `C n^-2`. Again both sums are bounded. Polynomial chain
coefficients are diagonal maxima with polylogarithmic fixed moments.
This proves the assertion, including every nonsingle-third-chain term.

## 5. All three derivatives on F

Summing the local-weight root in this part of (3) gives
`C=B diag(M(Y)) B`, with Frobenius norm `O_Lp(sqrt(n) polylog(n))`.
Thus the `R_F` part of (4) is negligible by Frobenius Cauchy--Schwarz.

For each nonedge single-third-chain term, the exact result proved in
Section 5 of the derivative module is

`n^-1 E sum_kj C_kj F_T(Y_k) D_j^3 Y_T,k = o(1)`.

Its proof uses (2), own-coordinate exclusion of the entire row `l` of
`K_3`, and ordinary-coordinate Gaussian integration by parts. It does not
assume the third-derivative matrix itself is small in operator norm.

Only the first term of (4) remains. Set `R=Q^(circ3)` and
`P_ik=B_ik (B R)_ik`. Its contribution is

`(sqrt(6)n)^-1 sum_ik P_ik E[M(Y_i) F_ggg(Y_k)]`.              (5)

The two-root Gaussian moment theorem for the old fields, followed by the
uniform `Y-X` approximation, gives

`E[M(Y_i) F_ggg(Y_k)] = E M(gamma) E F_ggg(gamma)`

`                                    + O(|Q_ik|) + o(1)`     (6)

uniformly for distinct roots. The diagonal is harmless because `B_ii=0`.
The constants in (6) depend only on the fixed responses; a joint Gaussian
polynomial moment with cross covariance `q I` is a polynomial in `q`,
whose deviation from its value at zero is bounded by `C|q|` on `[-1,1]`.

The weights obey

`n^-1 sum_ik |P_ik|=O_L(1)`,

`n^-1 sum_ik |P_ik Q_ik|`
` <= [n sqrt(n-1)]^-1 ||B R||_F ||Q||_F = O_L(n^-1/2)`.

Indeed `||R||op<=L^2` by the correlation Schur bound, so
`||BR||_F=O_L(sqrt(n))`. Equations (5)--(6) therefore factor their local
weight. Finally

`sum_ik P_ik = tr(Q R) = sum_ik Q_ik^4`,

`E F_ggg(gamma)/sqrt(6) = E[F(gamma)h3(gamma_edge)]`.

This is precisely the right side of (1).

## 6. All three derivatives on M

This allocation is distinct from Section 5 and must not be silently
identified with it. Put `V=B F(Y)`. Its contribution, apart from the fixed
factor `1/sqrt(6)`, is

`n^-1 E sum_ij V_i B_ij D_j^3 M(Y_i)`.                         (7)

The matrix `diag(V)B` has Frobenius norm `||V||_2`, with fixed moments
`O(sqrt(n))`, since `||B||op<=L` and the local responses have bounded
moments. Thus the remainder `R_M` in (4) contributes `o(1)`.

The all-edge term also vanishes. Uniformly in `i`,

`|sum_j B_ij Q_ij^3| <= [sqrt(n-1)]^-1 sum_j Q_ij^2`
`                                      <= L^2/sqrt(n-1)`.

Cauchy--Schwarz using `E||V||_2^2=O(n)` and bounded moments of `M_ggg`
therefore bounds that term by `O(n^-1/2)`.

It remains to treat a single third derivative of a nonedge tree. Write
`U_i=M_T(Y_i)V_i`. The root-hit part of (2) has Frobenius fixed moments
`O(polylog(n))`; pairing with `diag(U)B`, whose Frobenius moments are
`O(sqrt(n)polylog(n))`, again gives `o(1)`.

The root-not-hit part is exactly

`n^-1 sum_l E N_l [B diag(U) B K_3^T]_ll`.

Integrate `N_l` by parts. Differentiating `K_3` gives zero, because only
its own-free row `l` occurs in this diagonal entry. If
`J_U(i,l)=partial_(N_l) U_i`, the remaining expression is exactly

`n^-1 E <B circ J_U, B K_3^T>_F`.                             (8)

We claim `||J_U||_F=O_Lp(sqrt(n)polylog(n))`. Its two parts are

`J_U(i,l)=V_i partial_l M_T(Y_i)`
`                     + M_T(Y_i) [B J_Fcoord]_il`.

For the first part, fixed-degree hypercontractivity gives
`||V_i||_(2p)<=C_p||V_i||_2`, while each coordinate derivative has
`L^(2p)` size `O(n^-1/2)`. Consequently its squared Frobenius `L^p`
bound is at most `C sum_i ||V_i||_2^2=O(n)`.
For the second part, diagonal polynomial maxima have polylogarithmic
moments and `||J_Fcoord||_F=O_Lp(sqrt(n))` by its pointwise coordinate
influence bounds. This proves the claim.

Flatness now gives `||B circ J_U||_F=O_Lp(polylog(n))`, whereas
`||B K_3^T||_F=O_Lp(sqrt(n)polylog(n))`. Expression (8) is `o(1)`.
This finishes every term in (7), and every derivative allocation in (3).

## 7. Rademacher input and bounded response extension

The whole-functional comparison in
`fresh_weighted_unmarked_rademacher_bridge_2026_09_05.md` transfers the
polynomial theorem from Gaussian to Rademacher inputs, at error `O(n^-1)`.
Its input derivatives are taken before freezing a replacement coordinate.
The important exact third derivative is

`partial_k^3 Z_i = sqrt(6) Q_ik/(n-1)`;

its vector norm is `O(n^-1)`, although its diagonal entry is not
`O(n^-3/2)`. This precise vector estimate controls the fourth derivative
of the complete energy and justifies the comparison.

Assuming the already established uniform single-root joint convergence
of `(X_i,Z_i)` to `(gamma,sqrt(v_i)N_0)`, with `N_0` independent of the
old family and `0<=v_i<=L^2`, (1) also holds for fixed bounded measurable
`M,F` continuous outside Gaussian null sets. Use polynomial approximation
and bounded-operator Cauchy--Schwarz; independence gives the limiting
weighted error `E[(M-P)^2 Z^2]=v_i E(M-P)^2`. Fixed polynomial moment
bounds give uniform integrability. First take `n` to infinity at fixed
approximants, and only then remove approximation error.

The theorem therefore covers a fixed finite hard threshold/cutoff
construction, but does not itself supply uniform estimates as the number
of old coordinates, polynomial degree, or cutoff parameters vary.

## 8. Consequence and scope

If `M>=0`, `E M>0`, and the cubic edge coefficient
`b3=E[F(gamma)h3(gamma_edge)]` is nonzero, choosing the sign of `Z`
according to `b3` gives a strictly positive weighted direction, uniformly
over the bounded-operator signings:

`liminf n^-1 E [M(X)circ sign(b3)Z]^T B F(X) >= E M * |b3|`.

Root's exact paired means `mu_+=F_t+S H+epsilon D` and
`mu_-=-F_t+S H+epsilon D`, with `D=M clip_R(sign(b3)Z)`, cancel every
quadratic perturbation term in their half energy difference. Localizing
`M` to a fixed slack region of `F_t` therefore supplies a cube-feasible
gain after tail control. This is a finite-construction, fixed-`L` result.
A gain for every finite approximation is not automatically a gain over
an infinite-dimensional certified supremum: its approximation error
must still be compared with the dimension-dependent gain constants.

The variance-normalized version in (1a), together with the separate
inequality `n^-1 sum_i sqrt(v_i)>=1`, provides constants controlled by a
compact neighborhood of the limiting Gaussian covariance triple. This
additional argument, developed in the normalized-gain artifact, does
permit the finite approximation error to be chosen below the gain before
taking the matrix-size limit. It is not an unproved interchange of limits.
