# Self-audit: connected Walsh flux packing

Audit target:
[`walsh_connected_flux_packing.md`](walsh_connected_flux_packing.md).

Verdict: **PASS, with the normalization and scope qualifications already
stated in the draft.**  The construction gives an `h`-bit scalar-response
packing at fixed one-Walsh-block accuracy.  The bounded-degree and dense
queries really have one common connected support, and their unit-strength
connector edges do not consume any of the local triangle gap.  The only
gap loss comes from deliberately activating every non-target child, whose
total absolute mass is `9/200`.

## 1. Frozen state and independent flux coordinates

In chunk `i`, write the local triple as `(u,v,u+v)`.  In both states `u`
and `v` are independent and even.  Therefore

```math
c_1u+c_2v+c_3(u+v)=0
```

if and only if

```math
c_1=c_2=c_3.
```

Different chunks have disjoint coordinate supports, so a sum cannot cancel
between chunks.  Hence the complete relation kernel is exactly

```math
\mathcal R_h
=\{c:c_{i1}=c_{i2}=c_{i3}\text{ for every }i\}
=\langle111\text{ on chunk }i:i\le h\rangle,           \tag{ACF.1}
```

independently of `sigma`.  This is a dimension-`h` subspace with `2^h`
elements; it is not merely a list of `h` known relations.

Every label is supported on the first `4h` coordinates of
`V=F_2^(4h+1)`.  Every linear combination therefore has last coordinate
zero, while `omega=(1,...,1)` has last coordinate one.  Thus

```math
\mathcal R_omega(a^sigma)=\varnothing                  \tag{ACF.2}
```

for the **entire** coefficient space, not only for masks already in
`mathcal R_h`.

All labels have even weight, so every diagonal Gram entry is zero.
Different chunks are orthogonal.  In a zero-flux chunk all three
off-diagonal entries are zero; in a unit-flux chunk all three are one.
Because `u+v` is the third label and both diagonal values vanish,
bilinearity forces the three off-diagonal values to agree.  Thus each chunk
contains exactly one allowable Gram coordinate.  Choices in disjoint
chunks are independent, yielding the full cube `{0,1}^h` after every other
displayed state coordinate is frozen.

This is not a collision under the full sufficient state `(G,R)` from
Theorem 21.18: the point is precisely that `G` has `h` varying bits.  It is
a collision after deleting those `h` relation-cycle flux coordinates while
retaining all self-pairings, relations, root data, and cross-gadget Gram
entries.

## 2. Local theorem transport

The two triples in each chunk are coordinate permutations of the audited
tuples in Theorem 21.21.  Since `m=4h+1>=5`, that theorem applies without a
small-dimension exception.  A coordinate permutation of `V`, extended to
the corresponding coordinates of `E=V direct-sum V`, preserves the Walsh
matrix and the Boolean cube.  Therefore the local values are exactly

```math
M_0={9\over2},\qquad
M_1\le {3(1+\sqrt {17})\over4},                         \tag{ACF.3}
```

in units of `n^(3/2)`, and the zero-flux triple has one Boolean vector `x`
fixed by `F_E` and all three target children.  The proof needs no common
witness for a unit-flux triple or for any non-target child.

## 3. Connector saturation and non-destruction of separation

Let `H` have nonnegative weights and total weight `B_H`.  For every Boolean
pair,

```math
|q x^TF_Ey|\le q\|x\|_2\|y\|_2=qn=n^{3/2}.            \tag{ACF.4}
```

Thus every global assignment has connector value at most
`B_Hn^(3/2)`.  In the zero-flux target state, assign its common self-dual
witness `x` to all `3h` blocks.  Then every connector, regardless of cycles
or distance from the target, equals its positive individual ceiling because

```math
x^TF_Ex=x^Tx=n.                                        \tag{ACF.5}
```

The target triangle also reaches all six individual ceilings.  This proves
the exact favorable value `(M_0+B_H)n^(3/2)`.

In the unit-flux target state, restrict any global assignment to its target
three blocks.  That part is at most `M_1n^(3/2)` by Theorem 21.21, while the
connector part is separately at most `B_Hn^(3/2)`.  Sharing target variables
with the connector does not invalidate the sum of these pointwise upper
bounds.  Therefore the unfavorable value is at most
`(M_1+B_H)n^(3/2)`.

Consequently arbitrary nonnegative bridge padding preserves at least the
full local gap.  No assumption that the connected optimum decomposes is
being made.  Equality is proved only on the favorable side; a connector may
make the unfavorable state still worse, which can only increase the
separation.

Nonnegativity is essential to the displayed saturation argument.  With
mixed connector signs, one common self-dual vector need not attain every
signed edge ceiling around a cycle.  The theorem does not claim otherwise.

## 4. Active-label perturbation and the factor of two

One scaled non-target triangle has three child terms and three bridge
terms.  Since `C_a` and `F_E` are normalized orthogonal matrices,

```math
\left|{q\over2}x^TC_ax\right|\le {1\over2}n^{3/2},
\qquad
|q x^TF_Ey|\le n^{3/2}.                                \tag{ACF.6}
```

The absolute bound for one whole triangle is therefore

```math
3\left({1\over2}\right)+3(1)={9\over2}.               \tag{ACF.7}
```

With `h-1` non-target triangles and
`gamma_h=1/[100(h-1)]`, their total uniform perturbation is

```math
\epsilon_0
={9\over2}(h-1)\gamma_h={9\over200}.                   \tag{ACF.8}
```

A uniform perturbation of size `epsilon_0n^(3/2)` changes either maximum by
at most that amount.  Comparing two different states uses one lower bound
and one upper bound, so the local gap loses **twice** (ACF.8):

```math
\Delta_*
={3(5-\sqrt {17})\over4}-2{9\over200}
={3(5-\sqrt {17})\over4}-{9\over100}
=0.567670780786\ldots.                                 \tag{ACF.9}
```

The information-decoding threshold loses another factor two: uniform
error strictly below `Delta_*n^(3/2)/2` cannot assign one codeword to two
response states.

For `h=1` there are no non-target triangles and no perturbation.  The full
local gap is available, so the weaker uniform constant (ACF.9) remains
valid.

## 5. Query topology and activity

For `h>=2`, the bounded-degree connector path orders vertices by port index
first and gadget index second.  Consecutive vertices belong to different
gadgets, including the two transitions between port groups.  Hence its
`3h-1` unit edges are disjoint from all local triangle edges.  A path vertex
has degree at most two and receives two local triangle neighbors, so the
union has maximum degree at most four.  Every edge and onsite coefficient
lies in `(0,1]`.

For the dense version, `H` contains the nine unit cross edges for each pair
of gadgets, and the local triangles fill in all within-gadget pairs.  The
support graph is exactly `K_(3h)`.  Its connector weight is

```math
B_H=9{h\choose2}={9h(h-1)\over2}.                      \tag{ACF.10}
```

For both versions, all child coefficients are positive when `h>=2`: one
target gadget has coefficient one and every other gadget has coefficient
`gamma_h`.  Thus every marked label genuinely enters a child, in addition
to every vertex lying in the one connected bridge support.  At `h=1`, all
three labels are target labels and have coefficient one.

All `h` queries use the same support set and the same connector graph.  Only
which local triangle receives coefficient one changes with the public query
index.  No query coefficient depends on the hidden state `sigma`.

## 6. Packing and normalization

For distinct `sigma,tau`, choose a differing coordinate `i`; query `i`
orients the zero-flux state at least `Delta_*n^(3/2)` above the unit-flux
state, regardless of all other bits.  Therefore the response vectors of the
`2^h` states are pairwise separated in `l_infinity`.  This proves an
`Omega(h)`-bit lower bound for one decoder answering the fixed family of
`h` scalar queries.

Here

```math
n=2^{2m}
```

is one Walsh-block order and `n^(3/2)` is the normalization used in
Theorems 21.19 and 21.21.  The connected query has `3hn` Boolean variables.
If it is instead normalized by the total-variable scale `(3hn)^(3/2)`, the
gap coefficient is `Delta_*/(3h)^(3/2)`.  Thus this is a fixed one-port
accuracy packing, not a constant relative free-energy-density gap for the
whole growing support.

The path connector adds a public `Theta(h)n^(3/2)` ceiling and the dense
connector a public `Theta(h^2)n^(3/2)` ceiling.  The proof compares states
before any division by those totals.  Subtracting a public baseline from
all decoded responses is harmless, but no claim is made that the bad-state
connector contribution equals that baseline.

## 7. Verifier coverage and residual scope

The exact verifier reports

```text
connected Walsh flux packing checks passed: 152
```

and checks:

- all `2^h` tuples through `h=5`, including the full relation and root
  fibres;
- the local Gram cube and vanishing cross-gadget Gram blocks;
- the exact symbolic triangle sector polynomials;
- path connectivity, degree four, dense completeness, and all coefficient
  ranges;
- the exact `9/200` perturbation and (ACF.9).

It does not re-optimize the exponentially large Boolean landscapes.  The
only non-structural input is the independently audited and separately
verified Theorem 21.21.  The scalable proof uses termwise orthogonal bounds
and the theorem's explicit Boolean witness, so finite enumeration is not
being extrapolated to general `h`.

The result exposes one Gram bit per disjoint relation cycle.  It does not
show entrywise scalar minimality of a general Gram matrix, does not use a
single scalar query to encode `2^h` levels, and does not justify unit weight
on all unrelated label-dependent triangles.  Those limitations do not
affect the claimed connected `h`-query fixed-gap packing.
