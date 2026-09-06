# A literal vector seed weave and its unitary-edge obstruction

Date: 2026-09-06. Status: exact construction and exact limitation of its
local Hilbert-norm contraction. This tests a vector-source extension after
the director's precision-alignment breakthrough. It does not rule out all
vector/block constructions or prove nonconvergence.

## 1. Actual full sign matrices, with the seed retained jointly

Let A be a symmetric full sign matrix of order k (complete a hollow seed
with any sign diagonal). Let m=kd be a Hadamard order. Index outer fibres
by (i,alpha), with i in [k], alpha in [d]. For each alpha choose a real
Hadamard H_alpha of order m, and share it across the k fibres of its group.
For alpha<beta choose one independent fair sign tau_(alpha,beta), and set

```math
S_{(i,\alpha),(j,\beta)}=A_{ij}\tau_{\alpha\beta}
\quad(\alpha\ne\beta).
```

Choose the within-group entries of S as any symmetric full signing.
The literal weave

```math
K_{((i,\alpha),a),((j,\beta),b)}
=S_{(i,\alpha),(j,\beta)}
  H_\alpha(a,(j,\beta))H_\beta(b,(i,\alpha))
```

is a symmetric full sign Hadamard of order m² by the same exact weave
identity. Reusing bases within a group does not invalidate that identity.
Retain a common row selector T_alpha of size ell in each of its k fibres,
obtaining a genuine principal full-sign matrix of order N=m ell.

For a retained spin assignment let
`h_(i,alpha)=H_alpha[T_alpha,:]^T x_(i,alpha)`. The k spins within a
group may be arbitrary and correlated; their physical-row source is a
JOINT k-dimensional finite law, not k independent scalar laws.

For each unordered group edge alpha,beta form k-by-k blocks

```math
U^{\alpha\beta}_{ij}=h_{(i,\alpha)}((j,\beta)),\qquad
V^{\beta\alpha}_{ji}=h_{(j,\beta)}((i,\alpha)).
```

Define the seed-dependent signed transpose

```math
(R_A V)_{ij}=A_{ij}V_{ji}.
```

Because A is symmetric and its entries are signs, R_A is an orthogonal
involution on R^(k²).

## 2. Exact vector Gaussian kernel after shared-sign averaging

The usual squared defect sums both orientations of every fibre edge.
For a fixed group edge its exponential weight, after averaging the SINGLE
shared sign tau_(alpha,beta), is exactly

```math
K_A(U,V)
={1\over2}\left[e^{-s\|U-R_AV\|_F^2}
                   +e^{-s\|U+R_AV\|_F^2}\right]
=e^{-s(\|U\|_F^2+\|V\|_F^2)}
  \cosh(2s\langle U,R_AV\rangle_F),\qquad s=t/\ell.
```

Thus the seed has not been replaced by |A| or A² at this stage. It remains
inside a full joint vector interaction. This is a concrete exact sign
construction on which a proposed vector source/precision argument can be
tested.

## 3. Local Hilbert-norm contraction nevertheless erases the seed

Let phi_even be the even Gaussian-Fock feature map, so that

```math
\langle\phi_{\rm even}(U),\phi_{\rm even}(V)\rangle
=e^{-s(\|U\|_F^2+\|V\|_F^2)}\cosh(2s\langle U,V\rangle_F).
```

Every orthogonal map R acts unitarily on this feature space, via its even
tensor powers. Therefore

```math
K_A(U,V)=\langle\phi_{\rm even}(U),\mathcal U_{R_A}
                         \phi_{\rm even}(V)\rangle.
```

After averaging any A-independent row/base/permutation randomness, each
group is represented by a tensor on the product of its incident edge
feature spaces. Its source coordinates may retain their entire joint law.
Orient each group edge arbitrarily and absorb the unitary U_(R_A) into
one endpoint tensor's edge leg. This leaves that endpoint tensor's Hilbert
norm EXACTLY unchanged, even when all its edge legs are entangled.

Consequently the generalized graph Cauchy--Schwarz/Finner estimate by the
PRODUCT OF VERTEX HILBERT NORMS has exactly the same right side for every
seed A. There is no need to assume independent coordinates or scalar
channels to obtain this conclusion. Increasing the vector dimension or
using a matrix-valued Gaussian precision cannot restore seed dependence
after this particular norm bound has erased all edge unitaries.

Within-group defects can be discarded as nonnegative, just as diagonal
defects are discarded in the scalar weave proof. If instead retained,
they involve only O(k²d) fibre pairs versus O(k²d²) between groups. For a
fixed seed and d growing, they cannot by their count alone be treated as
a macroscopic seed-transport mechanism; any such use needs a new estimate.

## 4. Surviving exact target

A seed-sensitive bound for this literal construction must control the
NONLOCAL tensor contraction before replacing it by the product of its
local norms. One possible target would relate a strict contraction loss
to the Boolean quadratic support value Q(A) of the seed while remaining
uniform over every admissible joint group source. No such inequality has
been established here. The unitary-edge observation does not show that
the actual contraction or actual signing cap is independent of A; only
this entire local-Hilbert-norm upper-bound route is seed blind.
