# Whole-functional Boolean input comparison with a local unmarked weight

Date: 2026-09-05. This proves the input-comparison part of the weighted
unmarked-energy identity. The matrix remains the actual deterministic
signing throughout.

## 1. Statement

Let `B=A/sqrt(m)`, `m=n-1`, where `A` is a hollow symmetric signing and
`||B||op <= L` for a fixed constant. Put `Q=B^2`, `G(s)=Bs`, and

`Z(s)=B h3(G(s))`, where `h3(g)=(g^3-3g)/sqrt(6)`.

Let `X_i(s)` be a fixed finite vector of actual injective marked-tree
polynomials, and let `M,F` be fixed polynomial functions of that vector.
Define the full normalized functional

`T(s) = n^-1 [M(X(s)) circ Z(s)]^T B F(X(s))`.

For independent standard Gaussian inputs `N` and independent Rademacher
inputs `S`,

`E T(S) = E T(N) + O_(L,M,F,X)(n^-1)`.

Three derivatives and first-two-moment matching already give a sufficient
`O(n^-1/2)` comparison. The displayed stronger rate uses four derivatives
and matching first three moments. Neither result needs coherence of `Q`
or higher powers, nor any joint approximation over all pairs of roots.

## 2. Uniform hybrid input bounds

Use independent hybrid coordinates, each Gaussian or Rademacher. All
injective tree polynomials remain multilinear. For every fixed finite `p`,
coefficient counting and finite-degree hypercontractivity give

`sup_i ||X_i||_p <= C_p`,

`sup_i ||partial_k X_i||_p <= C_p n^-1/2`.

All higher derivatives of `X_i` in coordinate `k` are zero. These bounds
hold in hybrid product spaces: their multilinear monomials are orthogonal,
and Gaussian and Rademacher hypercontractivity can be applied successively.
It follows for any fixed polynomial local response `a` that

`sup_i ||partial_k^r a(X_i)||_p <= C_p n^-r/2`, for fixed `r`.

When coordinate `k` is fixed at `t`, the same bounds hold with a fixed
polynomial factor in `1+|t|`. This follows by expanding the affine tree
fields in that coordinate. Such factors have finite moments under both
replacement laws and justify the integral Taylor remainder below.

The transported cubic has uniformly bounded fixed moments. One direct
verification expands its cubic kernel

`K_i = sum_j B_ij b_j^(tensor 3)`.

Its squared Hilbert norm is `v_i=(B Q^(circ3)B)_ii <= L^2`. Repeated slots
have squared Hilbert norm `O_L(1/n)`, using

`K_i(a,a,b)=(Q_ib-B_ia B_ab)/m`.

Expansion in the hybrid Gaussian-Hermite/Rademacher orthogonal basis
therefore bounds the second moment, and finite-degree hypercontractivity
bounds every fixed moment. Equivalently the same hybrid moment bound
follows by expanding the fixed cubic and using that the coordinate
moments through order six are uniformly bounded.

## 3. Derivatives of the transported cubic

All derivatives in this section are ordinary derivatives with respect to
the coordinate currently being replaced, before it is frozen at `t`.
The first derivative is

`Z_i' = sqrt(3/2) sum_j B_ij B_jk (G_j^2-1)`.

For weights `w_j=B_ij B_jk`, the centered quadratic form has coefficient
matrix `B diag(w) B`. Its Frobenius norm is at most
`L^2 ||w||_2 = O_L(n^-1/2)`. In every hybrid model the quadratic form's
mean is zero and its variance is at most a fixed multiple of this squared
Frobenius norm. Thus

`sup_i ||Z_i'||_p = O_(L,p)(n^-1/2)`.

The next two derivatives have useful exact forms:

`Z_i'' = (sqrt(6)/m) [(Q s)_i-B_ik G_k]`,

`Z_i''' = (sqrt(6)/m) Q_ik`,

and `Z_i''''=0`. In particular

`sup_i ||Z_i''||_p = O_(L,p)(n^-1)`,

`||Z'''||_2 = O_L(n^-1)` deterministically.

The third derivative need not be `O(n^-3/2)` pointwise: `Q_kk=1`.
Its vector norm is the quantity needed in the comparison, and the exact
formula supplies that norm. Freezing the replacement coordinate at `t`
adds only polynomial-in-`t` factors to the preceding estimates.

## 4. Vector derivative scales for the complete local weight

Set `U=M(X) circ Z` and `V=F(X)`. For `0 <= r <= 4` and each fixed
`p >= 2`,

`|| ||partial_k^r U||_2 ||_p <= C_p n^((1-r)/2)`,

`|| ||partial_k^r V||_2 ||_p <= C_p n^((1-r)/2)`.

For derivatives of `V`, and for product-rule terms in `U` involving at
most two derivatives on `Z`, combine the pointwise estimates with

`|| ||W||_2 ||_p <= (sum_i ||W_i||_p^2)^(1/2)`.

The remaining product-rule terms have three derivatives on `Z`. Its
exact deterministic vector has coordinates `sqrt(6) Q_ik/m`. Therefore

`|| ||(partial_k^s M(X)) circ Z'''||_2 ||_p`

`<= (sqrt(6)/m) [sum_i Q_ik^2 ||partial_k^s M(X_i)||_p^2]^(1/2)`

`= O_(L,p)(n^-1-s/2)`.

For `s=0,1`, these are exactly the required third- and fourth-derivative
scales. This check prevents an incorrect pointwise estimate of `Z'''`
from entering the proof.

## 5. The whole weighted energy comparison

The fourth derivative of `T=n^-1 U^T B V` satisfies

`E|partial_k^4 T|`

`<= (L/n) sum_(r=0)^4 binom(4,r)`

`   * || ||partial_k^r U||_2 ||_2`

`   * || ||partial_k^(4-r) V||_2 ||_2`

`<= C_L n^-2`.

Every product before the outside factor `1/n` has order `n^-1`.
The polynomial-in-`t` versions of these bounds control the fourth-order
Taylor remainder uniformly through all hybrid replacement stages.
Gaussian and Rademacher variables have identical moments through degree
three. Replacing all `n` inputs consequently costs `O_L(n^-1)`.

If one uses the third derivative instead, each product before the factor
`1/n` has order `n^-1/2`, giving total comparison `O_L(n^-1/2)` by
matching only the first two moments.

## 6. Extension beyond polynomial responses

Suppose the single-root joint approximation of `(X_i,Z_i)` has already
been proved, with limit `(X,sqrt(v_i) N_0)`, where `X` is the old
independent Gaussian family, `N_0` is independent, and `0 <= v_i <= L^2`.
The bounded-operator-norm cubic module supplies exactly this approximation
uniformly in the output root, for both auxiliary input laws.

Then fixed bounded responses `M,F` that are continuous outside Gaussian
null sets follow by polynomial approximation. For example, for fixed
polynomial approximants `P,R`, operator-norm Cauchy--Schwarz bounds

`n^-1 |E [(M-P) circ Z]^T B F|`

by `L` times the product of their averaged Euclidean second-moment norms.
The first squared factor converges, uniformly along roots, to at most
`L^2 E(M(X)-P(X))^2`, by the joint approximation and its fixed polynomial
moment bounds. The term replacing `F` is treated the same way. Apply the
polynomial identity first, then take the approximation error to zero.

Thus the same limiting equality is valid for the local cutoff responses
used in a fixed finite construction. No cutoff, polynomial degree, or
smoothing parameter may grow during a dimension limit. This extension
uses local joint convergence only to control approximation errors; the
energy identity itself comes from the whole-functional comparison above.
