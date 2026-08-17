# Audit of scale-sensitive spectral bridge compression

Verdict: **REPAIR**.

**Re-audit after repairs: PASS.**  The repaired draft cleanly confines
arbitrary-future uniformity to pointwise bridge replacement, declares the
bounded retained-feature query family for the roof table, and disclaims any
bounded global quotient for arbitrary bounded-degree graphs.  One harmless
wording correction remains: a fixed number `d` of radius-`P` future ports has
total field radius `O_d(P)` (at most `dP`, or `sqrt(d)P` in an orthogonal
concatenation), not literally `P`; the constants should be said to depend on
the declared port bound.  The graph section already uses that convention.

The pointwise spectral replacement theorem is correct, with its stated sharp
constant, and the omitted-edge errors really are charged once.  The draft's
main defect is a scope conflation: pointwise bridge replacement is uniform
over arbitrary common futures, but the *compressed roof table* is sufficient
only for a declared family of continuations that factor through the retained
features with a controlled field radius.  It is false for arbitrary later
futures on the Boolean configurations.  The graph claim also needs to
distinguish a bounded local factor presentation from a bounded global
compositional quotient.

## 1. Constant and scaling audit

For `D=R-S`, every Boolean pair satisfies

```math
|x^TDy|\le \|x\|_2\|D\|_{2\to2}\|y\|_2
=\sqrt{pq}\|D\|_{2\to2}.
```

Taking a maximum after adding the same arbitrary function of `(x,y)` is
one-Lipschitz in this pointwise norm.  Thus (ST.2), (ST.3), and the sum over
physical edges in (ST.13) are correct.  The constant is existentially sharp:
take singleton `X={x_0}`, `Y={y_0}` (or landscapes that pin this pair) and
`D=sigma (x_0/sqrt p)(y_0^T/sqrt q)`.  Merely substituting this `D` does not
force equality for arbitrary fixed `X,Y,H,K`, so the sharpness sentence should
say explicitly that the examples are pinned or singleton.

At `p=q=n`, deleting all singular values at most `epsilon sqrt n` incurs at
most `epsilon n^(3/2)`, so (ST.6) is correct, including at the threshold.
For `\|R\|_{op}\le C\sqrt n`,

```math
\|\Sigma^{1/2}U^Tx\|^2\le \sigma_1(R)\|x\|^2
\le Cn^{3/2},
```

and hence (ST.9) is correct.  With the proposed `eta`,

```math
2P\eta\le {\epsilon\over4}n^{3/2},\qquad
\eta^2\le {\epsilon\over16}n^{3/2},
```

so (ST.9b) is also correct.  The height mesh contributes a further
`epsilon n^(3/2)/4`; these are *additional* to the spectral-tail error, so a
single advertised `epsilon` budget must rescale the separate allocations.

For fixed `C,epsilon`, the covering count is
`exp(O_C(r log(1/epsilon)))`; it is polynomial for `r=O(log n)` and
`exp(o(n))` for `r=o(n)`.  This counts response-table entries (and, after
height normalization, the same order of bits up to logarithmic factors), not
the precision needed to transmit the singular vectors or the bridge itself.
The bridge/factorization must be treated as shared query data.  Since `H` is
an arbitrary real landscape, one must also retain one unrestricted additive
baseline (for example `max H`) separately; only the relative roof heights
have the asserted finite alphabet.  Without either quotienting additive
constants or accounting for this scalar, no finite absolute-description
claim is possible.

## 2. Fatal reading and counterexamples

The sentence that (ST.2) is uniform over arbitrary common max-type futures is
correct **for replacement of the interaction function**.  It does not pass
to the roof quotient.

Take two child states with the same retained feature, heights `0` and `-L`.
The roof/bucket state discards the lower state.  A later common future that
adds `2L` only to that lower state changes the true optimum to `L`, while no
function of the retained feature can even express this query.  This works
already at rank zero.  Likewise, if two feature values are retained but roof
heights more than a fixed range below the top are clamped, an arbitrarily
large later linear field can expose the clamped point.

Therefore the finite table theorem must declare its query language:

* a single opposite child whose feature radius is at most `P`; or
* continuations depending on the child only through the retained feature and
  having total dual field norm at most `Q`; or
* a bounded-degree graph with a fixed bound on all present and future ports.

Under such a declaration, points more than `2PQ` below the top can be removed
and only `O(PQ/(epsilon n^(3/2)))` height levels are relevant.  In the
single-edge setting `P=Q=sqrt(C)n^(3/4)`, this is indeed
`O_C(1/epsilon)`.  Under repeated complete-graph composition or unbounded
future degree, `Q` grows and the claimed constant number of levels does not
survive.

The phrase "no full response landscape is stored" is also too strong.  The
upper roof *is* the complete response landscape for all retained linear-field
queries.  What is avoided is the full Boolean/state-specific response
landscape.

## 3. Lower-bound scope

Theorem 18.3 proves that general `r`-featured landscapes can require
`2^(Omega(r))` bits at fixed error.  Its construction allows an arbitrary
feature map and arbitrary queried fields.  The draft does not show that this
packing is realizable specifically by SVD features of Boolean bridge sides,
or that every exposing field is realized by the opposite Boolean child.
Consequently it supports the scoped sentence

> exponential rank dependence is unavoidable for general finite-feature
> roofs,

but not yet the stronger claim that the displayed spectral-bridge table is
optimal for Boolean SVD ports.  That stronger claim needs an embedding or a
separate lower bound.

## 4. Graph-composition audit

The exact truncated Hamiltonian can be represented as a factor graph.  At
vertex `v`, concatenate the endpoint features for its incident edges; its
local dimension is exactly at most `sum_(e ni v)r_e`.  The total objective is
the sum of local roofs and edgewise bilinear pairings, and multi-affinity lets
one replace each local state set by its upper roof.  Hence the *local factor
presentation* and the physical-edge error accounting are valid.

Bounded degree alone does not give a bounded global reusable state.  If one
eliminates a region, it generally creates a joint factor on every boundary
port; its dimension/metric entropy is governed by the cut or treewidth, not
the maximum vertex degree.  Associativity is valid if gluing retains all
exposed port variables.  It must not be read as saying that arbitrary graph
pieces can be summarized by the same bounded local table.  A bounded-degree
expander is the simplest stress test: all vertex factors are small while a
balanced region has linearly many exposed ports.

For a fixed bounded-degree graph and no undeclared future ports, local
quantization is sound after using the concatenated radii; constants acquire a
degree dependence.  If `d` is fixed, this remains
`exp(O_(C,d)(d_v log(1/epsilon)))` entries and total error is proportional to
the number of blocks/edges.  The draft should write this accumulated budget
explicitly.

## 5. Prediction audit

1. For `alpha I+beta J`, the `n-1` residual singular values equal `|alpha|`,
   so only the all-ones direction survives an `epsilon sqrt n` threshold for
   fixed `alpha`.  This numerical-rank claim is correct.  If `beta` is
   constant, however, `\|R\|` is order `n`, so the constant-size table bound
   under (ST.8) does not apply; the text currently claims only rank one and
   should keep that distinction explicit.
2. A bipartite matrix of uniformly bounded row and column degree has operator
   norm `O(1)` and contributes only `O(n)` at balanced size.  "Bounded-degree
   bridge" should mean both endpoint degrees are bounded.
3. A dense iid sign matrix has a positive fraction of singular values of
   order `sqrt n` (for thresholds in the bulk), so the extensive numerical-
   rank prediction is correct.  It should not be stated for arbitrary fixed
   `epsilon` above the spectral edge.

## 6. Classical versus new content

The operator-norm Lipschitz estimate, Eckart--Young truncation, volumetric
nets, and low-rank feature factorization are classical.  Exact roof algebra
and its general rank packing are already established earlier in this project.
The useful new contribution here is their scale-sensitive synthesis and the
clean physical-edge error accounting.  It is a valid theorem-level bridge
hierarchy after the query-family and graph-interface scopes are repaired; it
is not an all-purpose compositional compression theorem.

## Required repairs before promotion

1. Separate arbitrary-future uniformity of pointwise bridge replacement from
   the bounded feature-query sufficiency of the roof table.
2. State a field-radius/port-degree budget and use it in the height-clamping
   and quantization claims.
3. Describe ST.2 as an exact local factor presentation; disclaim bounded
   global quotient without bounded separators/treewidth.
4. Scope the `2^(Omega(r))` lower bound to general featured landscapes unless
   a Boolean-SVD realization is proved.
5. Clarify sharpness, cumulative epsilon budgets, bridge-as-shared-data, and
   the additive baseline and two minor prediction qualifications above.
