# Odd Walsh labels: exact rooted amalgamation and composition-created relations

Status: rigorous task-local draft.  It uses the odd-dimensional orbit
classification in
[`linear_label_walsh_gram_obstruction.md`](linear_label_walsh_gram_obstruction.md)
and is accompanied by an exact finite verifier.  The information lower
bounds below concern an **orbit-complete universal Walsh carrier**.  They are
not claimed to be lower bounds for the scalar maximum of one fixed graph.

## 1. Intrinsic rooted presentations

Work over `F_2`.  Let

```math
V=F_2^m,qquad B(x,y)=x\mathbin\cdot y,qquad
\omega=(1,\ldots,1),                                    \tag{RA.1}
```

and suppose throughout that `m` is odd.  Thus `omega` is the unique
characteristic vector:

```math
B(x,x)=B(\omega,x)\quad(x\in V),
\qquad B(\omega,\omega)=1.                               \tag{RA.2}
```

For an ordered label tuple

```math
\mathbf a=(a_1,\ldots,a_k)
```

let `E_a=F_2^k` and let

```math
\alpha:E_a\longrightarrow V,
\qquad \alpha(c)=\sum_i c_i a_i.                         \tag{RA.3}
```

Its intrinsic rooted presentation is

```math
\mathcal S(\mathbf a)=(G_a,R_a,Z_a),                    \tag{RA.4}
```

where

```math
G_a(c,c')=B(\alpha c,\alpha c'),\qquad
R_a=\ker\alpha,\qquad
Z_a=\alpha^{-1}(\omega).                                \tag{RA.5}
```

Here `Z_a` is either empty or one coset of `R_a`.  Put

```math
U_a=E_a/R_a.                                             \tag{RA.6}
```

The map induced by `alpha` embeds `U_a` as the label span.  The rooted orbit
classification from the preceding draft says that `S(a)` classifies the
ordered tuple under `O(m,2)`.  A simultaneous orthogonal change of every
label is a coordinate permutation of the Walsh cube, fixes the common Walsh
bridge, and conjugates every linear-label child.  Hence `S(a)` is an exact
carrier for all Walsh-graph Boolean landscapes on this tuple.

This carrier is algebraic: it is computed by binary dot products and Gaussian
elimination.  It contains no Boolean maximization or response table.

## 2. The relative datum of two pieces

Let `mathbf a` and `mathbf b=(b_1,...,b_l)` have coefficient maps `alpha`
and `beta`.  Their isolated states do not determine their relative position
inside `V`.  Define the following three pieces of cross-information.

First, the **cross form** is

```math
\kappa:U_a\times U_b\longrightarrow F_2,
\qquad
\kappa([c],[d])=B(\alpha c,\beta d).                    \tag{RA.7}
```

It is well defined because internal relations map to zero.  Second, let

```math
\delta:U_a\oplus U_b\longrightarrow V,
\qquad \delta(u,v)=\bar\alpha u+\bar\beta v,             \tag{RA.8}
```

and define the **coincidence correspondence**

```math
J_{ab}=\ker\delta.                                      \tag{RA.9}
```

Equivalently, `(u,v)` lies in `J_ab` precisely when the vector represented
by `u` in the first span equals the vector represented by `v` in the second.
Its two coordinate projections are injective, so it is the graph of an
isomorphism between the actual intersection subspaces as presented by the
two pieces.

Finally, define the **combined root fibre**

```math
Z_{ab}^{\times}=\delta^{-1}(\omega).                    \tag{RA.10}
```

It is empty or a coset of `J_ab`.  It can be nonempty even when both isolated
root fibres are empty.

Call

```math
\mathcal D(\mathbf a,\mathbf b)
   =(\kappa,J_{ab},Z_{ab}^{\times})                     \tag{RA.11}
```

the rooted amalgamation datum.

## 3. Exact gluing theorem

### Theorem RA.1 (rooted bilinear amalgamation)

For two actual odd-dimensional Walsh label tuples, their isolated states
and the datum (RA.11) determine the intrinsic state of the concatenated
tuple exactly.  Explicitly, on `E_a direct-sum E_b`,

```math
\begin{aligned}
G_{a\sqcup b}((c,d),(c',d'))
 &=G_a(c,c')+G_b(d,d')\\
 &\quad+\kappa([c],[d'])+\kappa([c'],[d]),              \tag{RA.12}\\
R_{a\sqcup b}
 &=\{(c,d):([c],[d])\in J_{ab}\},                      \tag{RA.13}\\
Z_{a\sqcup b}
 &=\{(c,d):([c],[d])\in Z_{ab}^{\times}\}.            \tag{RA.14}
\end{aligned}
```

Consequently this algebra gives the exact Boolean extremal response of the
joined tuple for every graph of common Walsh bridges, up to the same
coordinate relabelling of exposed Boolean queries.

Conversely, with the isolated states fixed, the combined rooted orbit state
determines all three objects in (RA.11).  Thus (RA.11) is the minimal
relative datum for an orbit-complete carrier: every cross-summary from which
the combined rooted orbit state can be recovered must determine
`kappa`, `J_ab`, and `Z_ab^times`.

#### Proof

The coefficient map of the concatenated tuple is

```math
\gamma(c,d)=\alpha c+\beta d.                           \tag{RA.15}
```

Expanding `B(gamma(c,d),gamma(c',d'))` gives (RA.12).
Moreover `gamma(c,d)=0` exactly when the two quotient-span vectors coincide,
which is (RA.13), and `gamma(c,d)=omega` is exactly (RA.14).  The odd-
dimensional rooted orbit theorem then gives the Walsh-landscape conclusion.

For the converse, `kappa` is the off-diagonal block of the combined Gram
form after passing to the two internal quotients.  The images of the combined
relation kernel and root fibre in `U_a direct-sum U_b` are respectively
`J_ab` and `Z_ab^times`.  Hence no orbit-complete relative carrier can forget
any of the three.

The assertion is deliberately about the universal rooted orbit carrier.
A particular graph or a scalar maximum can identify several rooted orbits;
no semantic lower bound for every individual cross bit is being smuggled
into the theorem.

### Corollary RA.2 (associative presented composition)

Rooted amalgamation is associative on actual tuples.  After amalgamating
`a` and `b`, its quotient coefficient space is canonically

```math
(U_a\oplus U_b)/J_{ab}
\simeq \operatorname{span}(\mathbf a,\mathbf b).         \tag{RA.16}
```

Amalgamating this presented span with a third tuple `c` produces the same
Gram form, relation kernel, and root fibre as a one-step amalgamation of all
three tuples, independently of parenthesization.

#### Proof

Both parenthesizations are the pullback presentation of the same map

```math
U_a\oplus U_b\oplus U_c\longrightarrow V,
\qquad (u,v,w)\longmapsto u+v+w.                        \tag{RA.17}
```

Its pulled-back bilinear form, kernel, and inverse image of `omega` are
independent of the order in which quotients by intermediate kernels are
taken.

This is a genuine composition law, but not an autonomous monoidal product of
the isolated states.  At each join the relative datum must be supplied.
That datum is exactly the information created by placing the two pieces in a
common Walsh coordinate system.

## 4. The three cross resources are independent

The following constructions show that the three entries in (RA.11) cannot
be collapsed to one another.  They also quantify the amount of
composition-created information in the orbit-complete theory.

Write

```math
H=\omega^\perp.
```

For odd `m`, `H` is a nondegenerate alternating space of dimension `m-1`.
Choose a symplectic basis `(p_i,q_i)` when needed.

### Proposition RA.3 (arbitrary cross-form information)

If `m>=2(r+s)+1`, fix

```math
a_i=p_i\quad(1\le i\le r).
```

For every binary `r by s` matrix `K`, put

```math
b_j^{K}=p_{r+j}+\sum_{i=1}^rK_{ij}q_i
       \quad(1\le j\le s).                              \tag{RA.18}
```

All first-piece states are identical, all second-piece states are identical,
`J_ab=0`, and the combined root fibre is empty.  Their cross forms are
exactly the `2^(rs)` matrices `K`.

Therefore an orbit-complete universal gluing carrier needs at least `rs`
bits of cross information in this family.

#### Proof

Both label families are linearly independent and internally totally
isotropic.  The private coordinate `p_{r+j}` makes the second family
independent and makes its span disjoint from the first.  Every label lies in
`H`, whereas `omega` does not, so all root fibres are empty.  Symplectic
orthogonality gives

```math
B(a_i,b_j^K)=K_{ij}.                                    \tag{RA.19}
```

Different matrices are different off-diagonal blocks of the combined Gram
state and hence different rooted orbits.

### Proposition RA.4 (quadratically many intersection bits)

If `m>=2r+1`, let `a_i=p_i`.  For every `P in GL(r,2)`, put

```math
b_j^P=\sum_iP_{ij}p_i.                                  \tag{RA.20}
```

The two isolated states, the cross form `kappa=0`, and the combined root
fibre are the same for every `P`, but

```math
J_{ab}^P=\{(Pv,v):v\in F_2^r\}.                         \tag{RA.21}
```

Hence an orbit-complete gluing carrier requires at least

```math
\log_2|GL(r,2)|
=\sum_{i=0}^{r-1}\log_2(2^r-2^i)
=r^2+O(1)                                               \tag{RA.22}
```

bits in the worst case even when the cross form vanishes.

#### Proof

Both ordered families are bases of the same totally isotropic space, so all
internal and cross Gram forms vanish and no span contains `omega`.  The
coincidence equation is `u=Pv`, which proves (RA.21).  Different `P` give
different combined relation kernels.  Finally

```math
\log_2|GL(r,2)|
=r^2+\sum_{j=1}^r\log_2(1-2^{-j}),                     \tag{RA.23}
```

and the last sum stays between two absolute constants.

The coincidence resource is not merely orbit bookkeeping: one bit of it is
already visible in the scalar Boolean maximum at the full leading scale.
For odd `m>=5`, take

```math
a=(1,1,1,0,\ldots,0),\qquad
b=(0,0,1,0,\ldots,0).                                   \tag{RA.23a}
```

Regard the two endpoint labels `(a,a)` as one two-port piece and the middle
singleton as a second piece.  The singleton choices `(a)` and `(b)` have the
same isolated rooted state, and their cross form against the endpoint span is
the same bit `1`.  All combined root fibres are empty.  Only the coincidence
correspondence differs: the middle label lies in the endpoint span in the
first gluing and not in the second.  Joining the middle port to both endpoints
gives the two path words `(a,a,a)` and `(a,b,a)`.  The exact bent witness and
anticommutator calculation in the preceding draft give

```math
\max E_{(a,b,a)}-\max E_{(a,a,a)}
\ge {7-3\sqrt3\over2}\,n^{3/2}.                         \tag{RA.23b}
```

Thus at least the equality/intersection bit detected in this example is
semantically necessary, not only necessary for orbit reconstruction.

### Proposition RA.5 (one irreducible root bit)

For odd `m>=5`, let

```math
a=e_1+e_2,qquad
b^+=\omega+a,qquad b^-=e_3.                            \tag{RA.24}
```

The isolated singleton states, cross form, and coincidence correspondence
are identical for `(a,b^+)` and `(a,b^-)`: the norms are `(0,1)`, the cross
pairing is zero, and both pairs are independent.  Nevertheless

```math
Z_{a,b^+}^{\times}=\{(1,1)\},
\qquad Z_{a,b^-}^{\times}=\varnothing.                  \tag{RA.25}
```

Thus the root fibre is not determined by the unrooted amalgam.

#### Proof

Equation `a+b^+=omega` proves the first claim.  The four elements of
`span(a,b^-)` are `0,a,b^-`, and `a+b^-`; none is `omega` when `m>=5`.
All other assertions follow by direct parity calculation.

## 5. Pairwise gluing data does not close repeated composition

### Proposition RA.6 (a ternary relation invisible on every pair)

Let odd `m>=7` and choose three symplectic basis vectors `p_1,p_2,p_3` in
`H`.  Compare

```math
(a,b,c^+)=(p_1,p_2,p_1+p_2),
\qquad
(a,b,c^-)=(p_1,p_2,p_3).                                \tag{RA.26}
```

Every singleton state agrees between the two triples.  For every pair of
positions, the cross form is zero, the pairwise coincidence correspondence
is zero, and the pairwise root fibre is empty.  Yet

```math
a+b+c^+=0,
\qquad a,b,c^-\text{ are linearly independent}.        \tag{RA.27}
```

Thus the complete collection of edge-local two-piece amalgamation data does
not determine the state of a multi-piece composite.

The associative law in Corollary RA.2 avoids the obstruction precisely by
presenting the accumulated span.  When `c` is attached, the new
coincidence correspondence is computed against
`span(a,b)`, not separately against the two old singleton states.  In this
sense higher mixed relations are composition-created memory, not pairwise
edge labels.

#### Proof

All displayed vectors are nonzero, distinct, totally isotropic, and lie in
`H`, proving the pairwise assertions.  Equation (RA.27) distinguishes the
global relation kernels.

## 6. Complexity and strictness

For a tuple of length `t`, the rooted presentation `(G,R,Z)` has an
`O(t^2)`-bit representation: store the binary Gram matrix, a row-reduced
basis of `R`, and either an absent marker or one representative of `Z`.
This is independent of the Walsh label dimension `m`.  Raw labels require
`mt` bits.  Consequently, for `t=o(m)` the carrier is strictly
sub-landscape and even strictly smaller than the raw linear-label
presentation.

If two accumulated spans have ranks `r` and `s`, their new cross form costs
`rs` bits, while the intersection correspondence and root fibre admit an
`O((r+s)^2)`-bit presentation.  Propositions RA.3 and RA.4 show that
quadratic cross-memory is sometimes genuinely necessary for orbit-complete
composition; it is not an artifact of the encoding.

This conclusion has two boundaries.

1. The state is an exact **symmetry quotient**, not an efficient formula for
   the scalar Boolean maximum.  A representative can be reconstructed by
   finite bilinear algebra, after which optimization may remain hard.
2. For `t` comparable with or larger than `m`, `O(t^2)` need not beat the
   raw `mt` label description.  The theorem identifies a real low-information
   regime; it does not assert uniform compression at arbitrary interface
   rank.

## 7. Noncircularity audit

The construction passes the information-theoretic checks relevant to the
project.

- `G`, `R`, `Z`, `kappa`, `J`, and `Z^times` are computed by polynomial-time
  linear algebra on the labels.  None uses a Boolean maximum, an optimizer,
  or an unknown target-order value.
- The exact response conclusion comes from an explicit global coordinate
  conjugacy, not from assuming equality of responses.
- The amalgamation datum does not reconstruct the raw labels.  It reconstructs
  only their orbit under a symmetry group which acts exactly on every Walsh
  graph landscape.
- The theorem does not pretend that independently summarized systems compose
  for free.  Propositions RA.3--RA.6 quantify the missing relative and
  higher-order information.
- No claim is made that every abstract triple `(kappa,J,Z^times)` is
  realizable in a prescribed ambient dimension.  The theorem classifies and
  composes actual realizations; adding a complete abstract realizability
  criterion would be a separate Witt-embedding problem.

This is therefore a positive strict quotient for a structured dense bridge
family and, simultaneously, a precise obstruction to autonomous local
composition.  The mathematical resource that survives is not pairwise Gram
data alone but a rooted presented span whose intersection relations are
updated at every join.

## 8. Exact verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_odd_walsh_amalgamation.py
```

The verifier checks (RA.12)--(RA.14) on exhaustive small tuples and random
larger tuples; realizes every cross matrix in Proposition RA.3; enumerates
the `GL(2,2)` intersection family; verifies the root-bit example; verifies
the pairwise-indistinguishable ternary relation; and checks equality of
direct and sequential relation/root presentations.
