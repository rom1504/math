# Independent audit: geodesic-fibre hard-core compression

**Verdict.** Accept GF.1--GF.3 as stated.  I reconstructed the proofs without
using the primary verifier and found no counterexample or hidden conversion
from state count to bit count.  The result is a genuine arbitrary-depth
quotient for spanning binary syndrome supports.  Its all-context extension
also covers raw nonspanning future supports.  It does not claim that a
nonspanning fragment can itself be encoded as an input state of the closed
algebra.

The independent finite checks are in
[`verify_phase3_geodesic_second_audit.py`](../experiments/verify_phase3_geodesic_second_audit.py)
and
[`phase3_geodesic_second_audit_results.json`](../experiments/phase3_geodesic_second_audit_results.json).

## 1. Reconstruction of the structural theorem

Let `t=sum_(b in B)b` have word length `D=|B|=D(S)`.  A nonempty linear
dependence inside `B` could be deleted from this representation, so `B` is
independent.  If `s in S cap span(B)` has `B`-coordinate support `I` with
`|I|>=2`, then

```math
t=s+\sum_{i\notin I}b_i
```

uses at most `1+D-|I|<D` generators.  Hence the zero fibre is exactly `B`.
Likewise, if distinct `s,s'` lie in one nonzero `span(B)`-coset and their
coordinate difference has weight at least three, then

```math
t=s+s'+\sum_{i\notin I}b_i
```

uses at most `2+D-|I|<D` generators.  This proves coordinate diameter at
most two for every nonzero fibre.  The argument applies to every diametral
vertex and every one of its shortest representations; no canonical choice
is required.

After translating a binary diameter-two anticode to contain zero, its
remaining points have weight one or two.  The weight-two supports are a
pairwise-intersecting edge family, hence are a star or a triangle.  Checking
the compatible singletons gives maximum size `D+1` for `D>=3`.  A full
two-cube has size four when `D=2`, so the stated exception is necessary.
Exhaustive enumeration independently gives maximum sizes `2,4,4,5` for
`D=1,2,3,4`, respectively.

## 2. Enumeration and units

For fixed ordered `B`, the zero fibre is forced.  Each of the
`2^(w-D)-1` other affine fibres is one of `a_D` diameter-two anticodes.  The
star/triangle classification gives

```math
a_D\le 1+2^D\left(2^D+D2^{D+1}+{D\choose3}\right)
    \le 4D2^{2D}\qquad(D\ge3).
```

Choosing an ordered geodesic basis in at most `2^(wD)` ways and summing over
`R<=D<=w` therefore gives exactly the logarithmic estimate GF.4.  This is an
upper count of pairs `(S,B)`, so multiple geodesics only overcount and cannot
invalidate it.

The complexity statement uses the correct units.  There are at most
`1+N_(w,R)` quotient states, so an enumerative index costs

```math
\left\lceil\log_2(1+N_{w,R})\right\rceil
```

bits.  GF.4 bounds this logarithm, not the number of states.  At
`R=floor(2 epsilon w+2)`, the fibre term is
`2^{(1-2epsilon)w+O(log w)}` bits and the basis term is `O(w^2)` bits.
The independent exhaustive count at `w=4,R=4` has 840 retained supports,
requiring 10 bits, while GF.4 gives the valid loose logarithmic bound 30.

## 3. Closed quotient and arbitrary futures

Diameter is antitone under support union.  Consequently the class

```math
I_R=\{S:D(S)<R\}
```

is an absorbing union ideal.  Retaining every support outside `I_R` exactly
and collapsing `I_R` therefore gives an associative, commutative, idempotent
Rees quotient.  Once a product enters the collapsed state, no later union
can leave it; there is no accumulated error.

For a collapsed spanning support every future union has integer diameter in
`[1,R-1]`.  Decoding by `R/2` thus has the sharp interval error
`R/2-1`.  With `R=floor(2 epsilon w+2)`, this is at most `epsilon w`, and the
threshold hypotheses in GF.10 indeed imply `4<=R<=w`.

The same reasoning handles every raw appended support `U`, spanning or not:
`S` already spans, and `D(S union U)<=D(S)`.  If `S` was collapsed, the
union remains collapsed; if `S` was retained, its exact support is available
and the union radius is computed exactly.  Exhaustive testing in
`F_2^3` covered all 92 spanning sources, all 128 raw future supports (36 of
them nonspanning), and both nontrivial thresholds: 23,552
source-threshold-future triples, including 8,320 collapsed cases.

## 4. Scope and residual gap

The theorem supplies a sub-landscape, depth-stable feature algebra at every
fixed positive relative error.  It does not determine the optimal response
metric entropy: the current lower bound is only quadratic in `w`, while this
upper bound is exponential in `(1-2epsilon)w`.  It also deliberately treats
future supports as raw queries when they are nonspanning; only spanning
fragments form the advertised closed source algebra.  Neither limitation
weakens GF.1--GF.3.
