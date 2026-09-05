# Surviving cubic channel with a local old-field weight

Date: 2026-09-05. Independent calculation of the surviving all-edge term
in the Gaussian weighted unmarked-energy expansion.

Let `B=A/sqrt(m)`, `m=n-1`, be an actual hollow symmetric signing with
`||B||op <= L`. Put `Q=B^2` and `R=Q^(circ3)`. Let `X_i` be a fixed
finite old marked-tree family, including `G_i`, and let `M,F` be fixed
polynomials in those fields. The old Gaussian limit is denoted by `X`.

## 1. The required two-root fact

For any two output roots `i,k`, the joint old-field moments are uniformly
approximated by two standard Gaussian vectors `(U,V)` with

`Cov(U_T,V_T') = 1_(T=T') Q_ik`.

Each individual vector has independent standard coordinates. The leading
tree-pairing proof gives cross-root covariance `Q_ik`: an isomorphic
whole-copy pairing has a common first free vertex, and summing that
vertex gives the row inner product of `B`. A free-free odd-parity edge
is negligible by the Boolean bilinear bound; root labels appearing as
spin vertices reduce the free-label count. This two-root fact is proved
and independently checked in the campaign's tree-moment audit. Here only
uniform convergence of fixed moments is needed, not a quantitative rate.

For any fixed polynomials `a,b`, their Gaussian covariance obeys

`|E a(U)b(V)-Ea(U)Eb(V)| <= C_(a,b) |Q_ik|`.

Indeed expand both polynomials in the multivariate orthonormal Hermite
basis. The matching multi-index of degree `d` contributes `Q_ik^d`; the
constant term is the product of the means, and Cauchy--Schwarz bounds the
remaining finite sum by `|Q_ik|` times the product of the centered `L2`
norms. No parity assumption is needed. If both polynomials are even,
one can replace `|Q_ik|` by `Q_ik^2`, but that improvement is unnecessary.

Consequently, uniformly in `i,k`,

`E[M(X_i) F_GGG(X_k)]`

`= EM(X) E F_GGG(X) + O_(M,F)(|Q_ik|) + o(1)`.

The same assertion holds with the one-level fields `Y` in place of `X`,
since their uniform fixed-moment difference tends to zero and the edge
coordinate is unchanged.

## 2. Exact coefficient matrix of the surviving term

Gaussian integration by parts in the transported `h3(G_j)` produces the
all-three-derivatives-on-the-edge-coordinate-of-`F` term

`T_main = (1/(sqrt(6)n)) sum_(i,j,k)`

`         B_ik B_ij Q_kj^3 E[M(X_i) F_GGG(X_k)]`.

Sum over `j` first, and define

`P_ik = B_ik (BR)_ik`.

Then `T_main=(sqrt(6)n)^-1 sum_ik P_ik E[M_i F_GGG,k]`.

The uniform `o(1)` in Section 1 can be inserted in this sum safely,
because

`n^-1 sum_ik |P_ik|`

`<= (n sqrt(m))^-1 sum_ik |(BR)_ik|`

`<= ||BR||F/sqrt(m) = O_L(1)`.

Here `||BR||F <= L ||R||F = O_L(sqrt(n))`, since `R` is a correlation
matrix with `||R||op <= L^2`.

The nonconstant Gaussian covariance contribution is also negligible:

`n^-1 sum_ik |P_ik| |Q_ik|`

`<= (n sqrt(m))^-1 ||BR||F ||Q||F`

`= O_L(n^-1/2)`.

Thus the response coefficient factorizes despite possible strongly
correlated pairs of roots. No coherence assumption is required.

## 3. The positive normalized trace

The remaining deterministic sum is exact:

`sum_ik P_ik = <B,BR>F = tr(QR) = sum_ik Q_ik^4`.

Therefore

`T_main = EM(X) [E F_GGG(X)/sqrt(6)]`

`         * [n^-1 sum_ik Q_ik^4] + o(1)`.

Gaussian integration by parts in the limiting edge coordinate gives

`E F_GGG(X)/sqrt(6) = E[F(X) h3(G)]`.

Equivalently the trace factor is the average local transported variance
`n^-1 sum_i (B R B)_ii`, which is at least one. This completes the exact
main-term calculation. The separate mixed-derivative and single-third-
chain estimates are still needed for the full weighted identity.

## 4. Final independent audit of the assembled theorem

The subsequently completed
`fresh_full_weighted_unmarked_projection_2026_09_05.md` passes this
independent audit, including all of its Sections 3, 4, and 6.

The mixed allocation with `a+b=3`, `a,b>=1`, is exactly the Frobenius
pairing `n^-1 <B circ A_a, B C_b>`. Flatness bounds the first factor by
a polylogarithm and the directional matrix-operator theorem bounds the
second by `sqrt(n)` times a polylogarithm. This treats the full mixed
chain rules, not merely their all-edge subterms.

The nonsingle-third-chain remainder has bounded, or polylogarithmic,
Frobenius moments by the stated squared-entry counts. In particular the
previously problematic second derivative of a non-edge field times an
edge derivative uses `sum Q_ij^2=O(n)`, not an absolute three-label sum.

For the allocation with all derivatives on the local weight, the last
own-coordinate integration by parts gives exactly

`n^-1 <B circ J_U, B K_3^T>`.

The first part of `J_U` is bounded using
`sum_i ||(BF)_i||_2^2=O(n)`, together with coordinate influences
`O(n^-1/2)`. It does not assume that each coordinate of `BF` has uniformly
bounded moments. The second part is controlled by the ordinary-coordinate
Jacobian Frobenius norm and polynomial diagonal maxima. Thus the claimed
`sqrt(n)`-scale Frobenius bound is valid, and flatness supplies the final
vanishing factor.

The complete theorem also accepts bounded deterministic root multipliers
`d_i`: these enter with the local weight, are never differentiated, and
multiply the norm estimates by at most `max_i |d_i|`. In Section 2 above,
replace `P_ik` by `d_i P_ik`; the exact surviving sum is then

`sum_ik d_i P_ik = sum_i d_i (B R B)_ii`.

Together with the separately audited whole-functional Rademacher bridge,
this proves the actual Boolean-input weighted identity and its bounded-
response extension. The dimension limit remains at fixed tree family,
responses, operator bound, and deterministic multiplier bound.
