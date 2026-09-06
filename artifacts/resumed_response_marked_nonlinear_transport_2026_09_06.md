# Marked nonlinear unmarked-channel tests: polynomial transport

Date: 2026-09-06. Author: resumed response track.
Status: complete proposed proof, submitted for independent reconstruction.
This is the marked counterpart to the director/bound-audit nonlinear
unmarked projection module. It does not assume that a linear covariance
identity determines a conditional expectation.

## 1. Statement and fixed parameters

Let `B=A/sqrt(n-1)` be symmetric hollow with flat sign entries and
`||B||op<=L`, where `L` is fixed. Let `X_i` be a fixed finite
ancestor-closed family of exact injective marked-tree fields. Let `F` be
a fixed odd polynomial of the old coordinates, and let `M` be any fixed
polynomial. Set

```
K(x)=sum_T E[partial_T F(gamma)] h_T(x)=U*F(x).
```

Let `Z_i` be a fixed finite real linear combination of odd unmarked
transports `B h_r(BS)`, with fixed odd degrees `r>=3`; denote its Gaussian
input variance by `v_i`. Let `psi` be a fixed real polynomial. Then

```
(1/n) E [S circ M(X) circ psi(Z)]^T B F(X)
 = E[M(gamma)K(gamma)] (1/n)sum_i E psi(sqrt(v_i)N)+o(1).     (1)
```

Bounded deterministic weights or bounded root-dependent coefficients of
fixed-degree `psi_i` may be included. All constants can depend on the
finite family, degrees, polynomial coefficients, and fixed `L`. No
construction parameter grows with `n`.

The same statement holds with independent Gaussian input `N_i` in place
of `S_i`, with the explicit own input on the left also Gaussian.

## 2. Coordinate Jacobian estimates

The global coefficient-flattening proof in
`fresh_tree_matrix_derivative_and_weighted_third_chain_2026_09_05.md`
also gives ordinary-coordinate derivative matrices

```
J_T(i,k)=partial_k X_T,i
```

polylogarithmic operator-norm moments. Indeed regard the differentiated
input index `k` as the second matrix index; every oriented coefficient
flattening is then an original global tree flattening. No contraction
by a row of `B` is necessary. The remaining stochastic indices are
square-free, so the same matrix-chaos theorem applies.

For a fixed polynomial old response `P`, its ordinary Jacobian is

```
J_P=sum_T diag(partial_T P(X)) J_T.
```

The diagonal coefficients are maxima of fixed-degree polynomials with
uniform moments, hence have polylogarithmic fixed moments. Thus

```
|| ||J_P||op ||_p <= C_p log(n+2)^C.                       (2)
```

Also, for each fixed `a>=1`, injectivity makes `X_T,i` affine in an
individual input coordinate, and coefficient counting gives

```
||partial_k^a P(X_i)||_p <= C n^(-a/2).                    (3)
```

The same first-derivative operator estimate holds for the one-level
fields `Y_T=B[N h_T(X)]`: their Jacobians are

```
J_YT=B diag(h_T(X))+B diag(N) J_hT.
```

In particular (2) applies to fixed polynomial responses of `Y` as well.

### Uniformity for Gaussian/Rademacher hybrids

The operator estimate for the square-free matrices `J_T` holds uniformly
for arbitrary mixtures of the two input laws. To check this without an
unstated hybrid matrix-chaos theorem, represent every Rademacher coordinate
as `sign(G_k)` and replace it in the Gaussian polynomial by
`sqrt(pi/2) G_k`. Its conditional mean given the sign is exactly that
sign. Multilinearity implies that the hybrid matrix is the conditional
expectation of the scaled Gaussian matrix. Jensen transfers every
operator-norm moment bound. Scaling a stochastic tensor index by the
fixed factor `sqrt(pi/2)` changes a fixed-order flattening bound by only
a fixed factor. Polynomial chain rules and hypercontractive maxima
then give (2) for arbitrary fixed polynomial responses under hybrids.

The pointwise estimates (3) likewise hold under hybrids, including every
fixed moment and the fixed-coordinate Taylor segments used below.

## 3. Gaussian proof: own-input integration by parts

Work with Gaussian input `N`. Replace only `F(X)` by `F(Y)` using the
one-level `O(n^-1/2)` fixed-root Lp approximation. Its error in the
normalized functional is `o(1)` by the fixed operator bound and uniform
moments. Keep `M` evaluated on the exact own-free fields `X`.

Write `A_i=M(X_i)psi(Z_i)`. Gaussian integration by parts in `N_i` gives

```
(1/n)sum_i E N_i A_i (BF(Y))_i
 =(1/n)sum_i E A_i partial_i(BF(Y))_i +o(1).               (4)
```

The derivative of `M(X_i)` with respect to its own input is exactly
zero. The derivative on `psi(Z_i)` has fixed moments `O(n^-1/2)`, from
the exact first-derivative formula for the odd unmarked transports.
Its normalized pairing with `BF(Y)` is therefore `o(1)` by Cauchy--Schwarz
and the global bound `E||BF(Y)||_2^2=O_L(n)`.

For each nonedge tree, the exact first derivative is

```
partial_i Y_T,j
 =B_ji h_T(X_i)+sum_l B_jl N_l partial_i h_T(X_l).
```

The edge has only the first term, with `h_edge=1`. Therefore the
root-hit part of the first term in (4) is

```
sum_T h_T(X_i) [(1/(n-1))sum_(j!=i) partial_T F(Y_j)].      (5)
```

Each empirical average in (5) tends in every fixed Lp to
`E partial_T F(gamma)`. The old two-root polynomial moment theorem
gives variance `o(1)`, and fixed-degree hypercontractivity upgrades the
convergence to all fixed moments. The excluded own summand costs `o(1)`.
The proved local joint limit of `(X_i,Z_i)` is the old Gaussian family
and an independent Gaussian of variance `v_i`, uniformly over roots.
Consequently (5), paired with `A_i`, gives the right side of (1).

## 4. Every deeper derivative term vanishes

For each tree `T`, the remaining term in (4) has the form

```
R=(1/n)sum_(i,l) E A_i C_il N_l D_li,
C=B diag(partial_T F(Y)) B,
D_li=partial_i h_T(X_l).
```

Exact own-coordinate exclusion of the entire child row gives
`partial_l D_li=0`. Integrating by parts in `N_l` therefore leaves
only the following two classes.

1. Derivative on `A_i`. The matrix with entries
   `(partial_l A_i)D_li` has bounded Frobenius Lp norm: both derivatives
   have pointwise moments `O(n^-1/2)`, so the squared entry sum is
   `O(1)`. Since `||C||F` has moments `O(sqrt(n)polylog(n))`, the normalized
   Frobenius pairing is `o(1)`.

2. Derivative on `C_il`. If `J(j,l)=partial_l partial_T F(Y_j)`, then

   ```
   [partial_l C_il]_(i,l)=B(B circ J).
   ```

   Flatness and (2) give this matrix Frobenius moments
   `O(polylog(n))`, because
   `||B(B circ J)||F <= L ||J||F/sqrt(n-1)`.
   The matrix with entries `A_i D_li` has Frobenius moments
   `O(sqrt(n)polylog(n))`. Their normalized pairing again is `o(1)`.

No derivative is omitted, and no pointwise bound is substituted for a
correlated Gaussian-direction estimate. This finishes the Gaussian proof.

## 5. A root-field moment bound for the marked input comparison

For a deterministic unit vector `w`, put `T_w=w^T F(X)`. Under any
Gaussian/Rademacher hybrid, oddness of `F` gives `E T_w=0` exactly.
The product Poincare/Efron--Stein inequality and (2) yield

```
Var(T_w)<=C polylog(n).                                  (6)
```

Here is the Rademacher-coordinate detail. Comparing its central difference
with its ordinary derivative at the actual input, (3) and Taylor's
formula give

```
||D_k T_w-partial_k T_w||_p
 <= C ||w||_1/n <= C/sqrt(n).
```

The sum of these error squares over `k` is bounded. The ordinary gradient
has squared norm at most `||J_F||op^2 ||w||_2^2`, establishing (6).
Fixed-degree product hypercontractivity upgrades (6) to every fixed
moment. In particular, taking `w` to be a row of `B`,

```
||(BF(X))_i||_p <= C_p polylog(n)                          (7)
```

uniformly in `i` and the hybrid stage. A frozen-coordinate version has
an additional fixed polynomial factor in its frozen value: first make
that coordinate Gaussian, project onto its finite one-variable Hermite
expansion, and evaluate that expansion at the frozen value. This does
not infer coefficients from values at only two sign endpoints.

## 6. Whole-functional sign input replacement, including the own-hit term

Consider the polynomial functional

```
T(s)=(1/n)[s circ A(s)]^T B F(X(s)),
A(s)=M(X(s))psi(Z(s)).
```

The existing odd-transport derivative formulas give, through derivative
order four,

```
|| ||partial_k^a A||_2 ||_p <= C n^((1-a)/2),
|| ||partial_k^b F(X)||_2 ||_p <= C n^((1-b)/2).
```

These bounds hold uniformly under hybrids and on Taylor segments, with
a fixed polynomial factor in the frozen coordinate. Products with the
explicit vector `s` preserve the first scale whenever no derivative
hits that vector, by pointwise moments and Holder. Thus all allocations
of four derivatives avoiding `s` give `O(n^-2)` normalized fourth
derivatives by the operator bound.

If one derivative hits the explicit `s`, the remaining factor is
supported on the single root `k`; no second derivative of `s` occurs.
Write `a+b=3` for the derivatives on `A` and `BF`. For `a<=2`,

```
(1/n) ||partial_k^a A_k||_p
        ||(B partial_k^b F)_k||_p
 <= C n^-1 n^(-a/2) n^((1-b)/2)=C n^-2.
```

For `a=3,b=0`, the cubic transport has the genuine exceptional diagonal
term `partial_k^3 Z_k=O(n^-1)`, not `O(n^-3/2)`. Using the root-field
bound (7) gives `O(n^-2 polylog(n))` even for this exception. For every
higher odd degree, the ordinary pointwise third-derivative scale suffices.

Consequently the complete fourth derivative has L1 norm
`O(n^-2 polylog(n))`, including its frozen-coordinate polynomial factor.
Gaussian and Rademacher laws match their first three moments. Replacing
all input coordinates therefore changes `E T` by
`O(n^-1 polylog(n))=o(1)`. This proves (1) for sign input.

## 7. Scope and handoff

Equation (1) is a fixed-polynomial theorem, including the even tests needed
to preserve the marked part of a nonlinear center response. The separate
unmarked nonlinear module must establish the complementary odd tests.
The bounded soft-sign/absolute-value application still requires the ordered
polynomial approximation and uniform-integrability handoff; those steps
are not replaced by claiming that linear covariances determine a law.

The same proof permits bounded deterministic weights and fixed-degree
polynomials whose root-dependent coefficients are uniformly bounded,
because no input derivative hits those deterministic quantities.
