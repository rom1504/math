# Switching-natural two-eigenvalue products reduce to the fixed tensor channel

Status: proved structural obstruction for entry-local Cartesian products. It
does not rule out products with a growing, nonlocal bridge state.

## 1. Derive the quadratic channel first

Let `A` and `B` be symmetric zero-diagonal sign matrices of orders `m` and
`n`. Consider a signing on pairs `(i,a)` whose sign is determined locally
from the relevant parent signs. Requiring covariance under arbitrary Seidel
switching of either parent forces the following rule. On pairs differing in
both coordinates, a Boolean function `F(A_ij,B_ab)` must obey

```math
F(ta,b)=tF(a,b),\qquad F(a,tb)=tF(a,b)\qquad(t=\pm1).
```

Hence `F(a,b)=pab` for one fixed `p in {+1,-1}`. On the two axial types the
same argument gives fixed signs `qA_ij` and `rB_ab`. Every complete,
entry-local, switching-natural product is therefore

```math
S_{pqr}(A,B)=pA\mathbin\otimes B+qA\mathbin\otimes I_n
                         +rI_m\mathbin\otimes B,       \tag{1}
```

with `p,q,r in {+1,-1}`. This includes the signed analogue of the strong
graph product; changing standard graph-product vocabulary only changes the
three signs or drops an axial channel, which would leave zero entries rather
than a signing of the complete graph.

For a Boolean state written as an `m` by `n` matrix `X`, let `u_a` be column
`a` and `v_i` row `i`. Its exact quadratic energy is

```math
x^{\mathsf T}S_{pqr}x
=p\sum_{a,b}B_{ab}\,u_a^{\mathsf T}Au_b
+q\sum_a u_a^{\mathsf T}Au_a
+r\sum_i v_i^{\mathsf T}Bv_i.                        \tag{2}
```

The last two sums are controlled by the parent quadratic caps:

```math
\left|q\sum_a u_a^{\mathsf T}Au_a+r\sum_i v_i^{\mathsf T}Bv_i\right|
\le n\kappa(A)+m\kappa(B),                            \tag{3}
```

where `kappa(A)=max_z |z^T A z|=2 cap(A)`. For balanced factor orders this
is power-saving relative to `(mn)^(3/2)`. The first term in (2), however, is
the full entangled tensor channel: its columns `u_a` need not be equal or
opposite, so it is not determined by either parent cap.

On product states `x=z tensor y`, (2) specializes exactly to

```math
(z^{\mathsf T}Az)(y^{\mathsf T}By)
+qn(z^{\mathsf T}Az)+rm(y^{\mathsf T}By),             \tag{4}
```

with the first coefficient multiplied by `p`. Thus even the product-state
part retains the fixed tensor channels, while a uniform upper bound must
also control all entangled Boolean matrices `X` in (2). The conference
operator norm gives only the leading-scale estimate

```math
|x^{\mathsf T}(A\otimes B)x|
\le mn\sqrt{(m-1)(n-1)}=(1-o(1))(mn)^{3/2},          \tag{5}
```

not a summable defect.

## 2. No nontrivial two-eigenvalue output in this class

Suppose now that both factors are symmetric conference matrices:

```math
A^2=(m-1)I_m,\qquad B^2=(n-1)I_n.
```

Squaring (1) gives

```math
S_{pqr}^2=(mn-1)I
+2pr(n-1)A\otimes I
+2pq(m-1)I\otimes B
+2qrA\otimes B.                                      \tag{6}
```

The four displayed matrices are Frobenius-orthogonal. In particular, (6)
is never `(mn-1)I` for nonzero factors, so the product is not conference.

More generally a real symmetric matrix has two eigenvalues only if
`S^2=tS+uI` for scalars `t,u`. Comparing the three nonconstant coefficients
in (1) and (6) would require

```math
t=2pqr(n-1)=2pqr(m-1)=2pqr.                           \tag{7}
```

Thus `m=n=2`. For all conference orders `m,n>2`, no complete entry-local
switching-natural product of the factors has two eigenvalues.

This explains why standard Kronecker/strong-product spectral theorems do not
give an all-order conference operation: filling the zero axial pairs creates
the two extra channels in (1), and their cross terms in (6) destroy the
two-eigenvalue identity.

## 3. Scope and surviving theorem obligation

Dropping switching covariance permits an arbitrary Boolean function on the
four input sign pairs. Its Fourier expansion adds the already familiar
`J-I`, row-sum, column-sum, and `A tensor B` channels; it does not remove the
entangled term in (2). Homogeneous association-scheme products add a second
problem: regular or asymptotically conference outputs expose the all-ones
Boolean channel and hence attain normalized cap `1/2-o(1)`.

Consequently an all-order product with a summable cap defect cannot come
from a fixed local graph-product table or from the two-eigenvalue spectrum
alone. A surviving operation must use a bridge depending nonlocally on the
parent matrices (or a state whose complexity grows with order) and prove a
uniform inequality controlling the first term of (2). That is precisely the
fixed tensor-channel obligation already isolated in the ledger, rather than
a new composition theorem supplied by graph-product terminology.
