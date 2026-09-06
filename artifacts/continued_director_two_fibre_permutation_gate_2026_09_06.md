# A uniform two-fibre signed-permutation gate

Date: 2026-09-06. Status: exact director derivation, submitted for independent
audit. This is a signing-preserving operation, not a cap or landing theorem.

## 1. Exact gate on an arbitrary seed

Let S be any full symmetric sign matrix with diagonal +1, and H any full
symmetric sign matrix of order m. H need not be Hadamard for this statement.
Select two seed indices a,b and set e=S_ab. On these two m-point fibres,
use the orthonormal seed coordinates p_+=(1,1)/sqrt(2),
p_-=(1,-1)/sqrt(2). For any TWO independently chosen signed permutation
matrices U_+,U_- of order m, let R act by U_+ in the p_+ channel and by
U_- in the p_- channel, and be the identity on all other seed fibres.

Then

```math
\boxed{K'=R(S\otimes H)R^T
\text{ is an exact full symmetric sign matrix.}} \tag{1}
```

It has exactly the same spectrum as S tensor H. This is valid for every
seed S and every m, without a relation between different seed rows.
An arbitrary hollow seed can be completed as S=A+I first.

**Proof.** The selected principal seed block is
[[1,e],[e,1]]=2p_e p_e^T. Thus its transformed two-fibre block is

```math
\begin{pmatrix}H_e&eH_e\\eH_e&H_e\end{pmatrix},
\qquad H_e=U_e H U_e^T,
```

which is a sign matrix. For another seed coordinate c, the link vector
(S_ac,S_bc) lies wholly in channel t=S_ac S_bc. Its transformed two
blocks are S_ac U_t H and S_bc U_t H, both sign matrices. Transposed
blocks follow by symmetry, and all other blocks are unchanged. This
checks every entry; orthogonality of R checks spectral preservation.

Hollowing K' changes the full quadratic half-energy by at most km/2,
where k is the seed order. This error is o((km)^(3/2)) when km grows.
No claim that hollowing preserves exact eigenvalues is intended.

## 2. Exact response constraint: where the interaction changes

For Boolean spins x_a,x_b on the selected fibres define
u_+=(x_a+x_b)/2 and u_-=(x_a-x_b)/2. They take values0,+1,-1 and obey

```math
|u_+|+|u_-|=\mathbf1 \quad\text{coordinatewise}.
```

Set v_t=U_t^T u_t. The selected-fibre internal part of the FULL quadratic
form is 4v_e^T H v_e. Its cross part with other seed fibres is

```math
4\sum_{c\notin\{a,b\}}S_{ac}\,v_{S_{ac}S_{bc}}^T Hx_c.
```

All other terms are unchanged. The altered compatibility constraint is

```math
|U_+v_+|+|U_-v_-|=\mathbf1. \tag{2}
```

Thus the relative permutation changes the matching between complementary
supports in the two channels BEFORE their joint energy is maximized.
Their signs remain independently free on their supports. Paying the two
channels independently erases precisely this constraint, and cannot by
itself analyze a benefit from the gate.

## 3. Relation to local rectangle switches and remaining obligation

For two outer indices u,v, choose a four-coordinate product sign vector
w=(1,sigma) tensor(1,tau). The reflection I-ww^T/2 is one instance of
this gate: in one seed channel it is a signed transposition of u,v and
in the other it is the identity. Products of such transpositions generate
the allowed U_t. This proves the entire two-fibre family directly, without
assuming that generic one-sided Hadamard switching preserves symmetry.

After this operation, applying a gate on a DIFFERENT overlapping pair of
seed fibres need not satisfy the initial tensor-block hypotheses. Closure
under arbitrary networks is not asserted. Likewise, preserving spectrum
does not prove preservation or improvement of the same-spin Boolean cap.

The concrete next question is whether the matched-support constraint (2)
admits a quantitative cap reduction for arbitrary growing minimizing seeds,
or whether an explicit surviving witness family prevents such a reduction.
A finite cap change alone would not settle that question. Any convergence
application still requires a proved all-order or summable-loss landing
inequality; none is assumed in (1).

## 4. A converse for preserving every pure Boolean tensor

The gate above has an exact witness-preservation property: if x=u tensor v
with u and v Boolean, each selected seed pair occupies one channel only,
and its outer vector is changed by a signed permutation. Thus Rx is again
Boolean. This preserves all scalar-tensor witnesses, not every Boolean
vector and not necessarily a pure tensor.

There is a sharp elementary restriction on ANY linear map with this
property. Let a row of a real linear map be indexed by seed/outer pairs,
with coefficients c_ij. Suppose

```math
\sum_{ij}c_{ij}u_i v_j\in\{-1,+1\}
\quad\text{for EVERY Boolean }u,v.
```

Then precisely one of the following holds:

1. One coefficient is +1 or -1 and all others vanish.
2. Exactly four coefficients are nonzero. They lie on a two-by-two
   rectangle, each has magnitude1/2, and their four signs multiply to -1.

**Proof.** Fix u. A homogeneous linear function of v taking only sign
values has exactly one coefficient equal to a sign and all others zero:
expand its square, or maximize and minimize the linear form. Consequently
each g_j(u)=sum_i c_ij u_i belongs to {0,+1,-1}. Flipping u_i shows
2c_ij belongs to {0,+1,-1,+2,-2}. Thus coefficients are half-integral.
Parseval for the original Boolean-valued bilinear function gives
sum_ij c_ij²=1. A coefficient of magnitude1 uses all this mass. Otherwise
there are exactly four coefficients of magnitude1/2.

View their support as four edges of a simple bipartite graph. If two are
disjoint, their contribution to the coefficient u_i u_k v_j v_l in the
square must be canceled by the only other possible pairing, using the
other two corners of the same rectangle. Therefore all four edges are
that rectangle. If there are no disjoint edges, all four edges form a
star; a pair of star edges creates an uncanceled degree-two coefficient
in the square, impossible. On the rectangle, cancellation in the square
is exactly the negative sign-product condition. Conversely a two-by-two
Hadamard sign pattern divided by2 is Boolean on pure sign inputs.

This theorem needs no orthogonality assumption on the full map. In
particular an arbitrarily long product of gates that STILL maps every
pure Boolean tensor to a Boolean vector has rows supported on at most
four original coordinates. Products do not inherit the property merely
because their factors each have it: the intermediate images need not be
pure tensors. Any claim of unlimited growth while preserving that exact
invariant must confront this restriction.

The classification is not a cap lower bound and does not kill a one-layer
landing construction. It identifies an exact algebraic resource that an
overlapping-network argument cannot silently keep. Both the director and
the independent convergence researcher reconstructed the classification.
