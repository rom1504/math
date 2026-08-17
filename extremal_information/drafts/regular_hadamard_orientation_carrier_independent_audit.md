# Independent audit: regular-Hadamard orientation carrier

**Audited files.**
`regular_hadamard_orientation_carrier.md` and
`verify_regular_hadamard_orientation_carrier.py`.

**Verdict.**  The orbit count, canonical invariants, gluing-fibre cardinality,
and cap separations are mathematically correct.  The main qualifications are
about language, not the algebra: the gluing coordinates are not canonical
until marginal sections/frames are chosen, the canonical verifier tests only
singleton gluing pieces, and the claim that the carrier is exponentially
smaller than the spin landscape needs a fixed-graph or `e=o(kn)` qualifier.
The draft otherwise states its response-minimality and WC.1 scope honestly.

An independent finite verifier is
`verify_regular_hadamard_orientation_carrier_independent_audit.py`.

## 1. Projective switching, including disconnected graphs

Encode signs by bits.  If `alpha` is the global-antipode bit and `d_i` are
vertex-switch bits, the action is

```math
s_i\longmapsto s_i+\alpha,
\qquad
t_{ij}\longmapsto t_{ij}+\alpha+d_i+d_j.             \tag{AOC.1}
```

Here `s_i` encodes `sigma_i` and `t_ij` encodes `b_ij`.  Because every
diagonal sign is nonzero, an element in the action kernel must have
`alpha=0`.  Then `d_i+d_j=0` on every edge, so `d` is constant on each of
the `c` connected components.  The kernel therefore has dimension `c`, and
the effective action has dimension

```math
1+k-c.                                                \tag{AOC.2}
```

It is free after division by that kernel.  The quotient of the `k+e`
coefficient bits consequently has dimension

```math
k+e-(1+k-c)=e+c-1,                                   \tag{AOC.3}
```

including isolated vertices and fully disconnected graphs.  Thus OC.8 is
correct.  In particular, the global antipode contributes exactly one
effective action bit even when `c>1`; there is not one antipode per component
in the joined carrier.

The energy statement is also exact.  If `T'=epsilon DTD`, then

```math
X^T(T'\otimes H)X
=\epsilon,((D\otimes I)X)^T(T\otimes H)((D\otimes I)X), \tag{AOC.4}
```

and `D tensor I` permutes the Boolean cube.  Hence the whole energy multiset,
with multiplicity, is preserved up to the sign `epsilon`.

## 2. Canonical coordinates OC.9

The displayed coordinates are invariant.  For every `v ne r_1`,
`sigma_v sigma_(r_1)` is unaffected by either part of the action.  On a
fundamental cycle `C` of length `ell`, vertex switches occur twice and
cancel, while the global antipode multiplies the bridge product by
`epsilon^ell`.  Therefore

```math
\sigma_{r_1}^{\ell}\prod_{f\in C}b_f                 \tag{AOC.5}
```

is invariant as well.  There are `k-1` orientation ratios and
`e-k+c` fundamental-cycle coordinates.

They are complete: normalize `sigma_(r_1)` first with the global antipode,
then recursively switch along each forest tree to make all forest edges
positive.  In that representative, the orientation ratios recover every
diagonal sign and (AOC.5) recovers every chord sign.  This gives a unique
representative modulo component-constant switches, which act trivially.
The independent verifier confirms that these coordinates and orbit keys are
in bijection for all 75 labelled graph supports through four vertices,
including disconnected supports.

There is one harmless order-of-operations blemish in the proof prose.  It
first says that the forest is normalized and *then* applies the antipode;
an antipode flips the forest edges too.  The canonical-coordinate paragraph
states the safe order (antipode first, forest switching second), and that is
the order the constructive proof should use.  Re-normalizing the forest
afterward would also repair it.

## 3. OC.2: surjectivity and the exact fibre

Let `Q(G)` be the projective-switching quotient.  Restriction gives a
well-defined map

```math
\rho:Q(G)\longrightarrow\prod_{a=1}^s Q(G_a),          \tag{AOC.6}
```

because the common joined antipode becomes an allowed antipode on every
marginal.  It is onto: choose an arbitrary coefficient representative of
each prescribed marginal class and assign arbitrary signs to the cross
edges.  Since all pieces and the joined graph are connected,

```math
\dim Q(G)=\sum_a e_a+r,
\qquad
\dim\prod_a Q(G_a)=\sum_a e_a.                         \tag{AOC.7}
```

The map is an affine linear quotient over `F_2`, so every nonempty fibre has
the same cardinality, namely `2^r`.  This proves both the unspoken
surjectivity step and OC.11 for arbitrary connected pieces.

After choosing representatives (sections) for the marginal classes, one
sees the split directly.  Independent marginal antipodes give `s-1` bits
modulo the common joined antipode.  Component-constant switches have rank
`s-1` on the connected cross-edge graph and normalize a cross spanning
tree.  The remaining `r-s+1` edge signs are fundamental-cycle transports.
The sum is `r`.

The **dimensions** of the two resources are intrinsic, but the coordinates
are not literally canonical after choosing only a cross-edge spanning tree.
One must also choose a representative/gauge section for each marginal
carrier (and a reference piece).  Different sections relabel the relative
antipodes and cycle transports without changing the fibre.  Thus “splits
canonically” should be read as “admits this coordinate split after the stated
forest and marginal-frame choices.”

The canonical verifier's function `check_gluing_fibres` tests singleton
pieces only, despite Section 7 saying it checks small connected pieces.  The
independent verifier adds three non-singleton cases: an edge plus singleton,
two edge pieces joined at two ports, and a triangle plus singleton joined at
two ports.  Every marginal tuple occurs and every fibre has exactly `2^r`
classes.

## 4. Spectra and cap normalizations

Because `H^2=nI`, its operator norm is `q=sqrt(n)`.  Every `k`-block Boolean
word has squared Euclidean norm `kn`, hence

```math
Q_H(T)\le {1\over2}\|T\|qkn.                          \tag{AOC.8}
```

For the two-block matrices,

```math
\operatorname{spec}(T_+)=\{2,0\},
\qquad
\operatorname{spec}(T_-)=\{\sqrt2,-\sqrt2\}.          \tag{AOC.9}
```

The all-regular pole attains `Q_H(T_+)=2qn`; (AOC.8) gives
`Q_H(T_-)<=sqrt(2)qn`.  OC.14 and its
`(2-sqrt(2))n^(3/2)` lower gap are correct.

For the triangle,

```math
\operatorname{spec}(J_3)=\{3,0,0\},
\qquad
\operatorname{spec}(T_{unbal})=\{2,2,-1\}.             \tag{AOC.10}
```

The regular pole attains `Q_H(J_3)=9qn/2`; (AOC.8) gives
`Q_H(T_unbal)<=3qn`.  The claimed universal gap is therefore at least
`3qn/2=3n^(3/2)/2`.  At `H_4`, exhaustive values are respectively 36 and
20 (the upper bound for the unbalanced word is 24), so the finite check is
consistent and actually has a larger gap.  The independent verifier checks
all eight bridge-sign words: positive bridge product always has spectrum
`(0,0,3)` and cap 36, while negative product always has spectrum
`(-1,2,2)` and cap 20.

## 5. Minimality and relation to WC.1

The scope language is mostly careful and should be retained:

1. `T/~ps` is minimal only as an injective encoding of the declared
   coefficient-conjugacy orbits.  It is not proved minimal for the scalar
   cap or for the full energy multiset; the draft explicitly says so.
2. OC.3 and OC.4 show that the two **types** of compatibility resource can
   each be macroscopically visible.  They do not show that all `e+c-1` bits
   are simultaneously exposed, nor a packing lower bound of that size.  The
   draft's final qualification is consistent with this limitation.
3. WC.1's Kronecker factorization remains the archived source of the basic
   algebraic closure.  OC restricts each onsite block to two orientations of
   one common Hadamard factor and each bridge to its signed scalar copy.  It
   does not compress WC.1's general truth-table labels or arbitrary bridges,
   and it makes no claim to do so.

Two phrases deserve tightening:

- “canonical minimal coordinate count” should include “for the declared
  projective-switching action”; otherwise it sounds like response
  minimality.
- Section 6 says `2^e` carrier states are “exponentially smaller than the
  `kn`-spin landscape.”  This is true for fixed `G` as `n` grows, or whenever
  `e=o(kn)`, but false uniformly over dense block graphs with `k` much larger
  than `n`.  The unconditional statement is only that the carrier uses `e`
  coefficient bits independent of Hadamard order.

With these qualifications, the note is a valid exact signed-graph gauge
model of composition-created orientation and cycle information, not an
overclaim about the original signing problem.

## 6. Reproduction

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_orientation_carrier.py

./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_orientation_carrier_independent_audit.py
```

Both pass.  The independent run reports 75 disconnected/connected support
classifications, 14 non-singleton marginal fibres across the three gluing
cases, and all eight triangle bridge words.
