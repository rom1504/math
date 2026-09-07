# Size-bounded pure-Max-Cut response profiles

## Scope and resource model

This note separates three promises which should not be conflated.

* The **interface promise** is unit outer load: for each labelled boundary
  vertex `i`, the total weight of edges of the component incident with `i`
  is at most one.
* The **semantic size promise** is that the component has at most `p` private
  vertices and at most `m` (nonzero, nonparallel) edges.  Edge weights are
  arbitrary nonnegative real numbers; in particular, this promise makes no
  bounded-precision assumption.
* A **description-length promise** says only that a chosen decoder receives
  `B` bits.  In the absence of restrictions on that decoder this has no
  nontrivial consequence beyond the cardinality bound `2^B`; a precise
  tightness statement is given below.

Let

```math
X_w=\{+1,-1\}^w/\{s\sim-s\}
```

with projective Hamming distance `d_proj`.  For a weighted graph `G` with
boundary `[w]`, let

```math
h_G([s])=\max_z \operatorname {Cut}_G(s,z).
```

Offsets are quotiented, because boundary-disconnected private components
can contribute an arbitrarily large constant.  Thus the metric is

```math
d_{sh}([f],[g])={1\over2}\operatorname {osc}(f-g).
```

Write `C_(w,p,m)` for the resulting response shapes under the first two
promises.  The proof below constructs internal covers (their centers are
actual response shapes), so the bounds also hold under the looser convention
allowing ambient centers.

## The main structural theorem

The proof is an instance of a general max-affine presentation lemma which is
useful beyond Max-Cut.

### Lemma (shared-parameter tropical presentation entropy)

Let `X` have `q` elements.  For each `x in X`, let `A_x` be a set of at most
`r` vectors in `R^d`, and consider the arbitrary-real-parameter family

```math
f_\theta(x)=\max_{a\in A_x}\langle a,\theta\rangle,
\qquad \theta\in\Theta\subseteq\mathbb R^d.       \tag{T.1}
```

Suppose every `[f_theta]` lies in the radius-`R` ball of
`R^X/R1` under the shape norm.  Then

```math
\boxed{
\log_2\operatorname {Cov}_\delta\{[f_\theta]:\theta\in\Theta\}
\le d\log_2\!\bigl(4(qr^2+1)\bigr)
 +d\log_2\!\left(1+{2R\over\delta}\right).}      \tag{T.2}
```

No boundedness, discretization, or finite precision is assumed for
`theta`.  The proof partitions parameter space by the at most
`q binom(r,2)` comparison hyperplanes, observes that `f_theta` is one
linear map on each face, and applies the normed-space volumetric cover on
that face.  Thus (T.2) is a semantic presentation bound: `d` shared real
parameters and `r` competing affine witnesses per query cannot carry an
arbitrary table of independently precise values.

If every witness vector belongs to `{0,1}^d`, there are at most
`(3^d-1)/2` distinct comparison hyperplanes, irrespective of `q` and `r`.
Accordingly (T.2) improves to

```math
\log_2\operatorname {Cov}_\delta
=O\left(d^2+d\log\left(1+{R\over\delta}\right)\right). \tag{T.3}
```

More generally a fixed finite coefficient alphabet gives the same
quadratic-in-parameter-dimension law, with an alphabet-dependent constant.

The same proof works with affine forms after adjoining a constant parameter,
and with minima after changing signs.  It therefore applies directly to
finite-state graphical models, finite-horizon robust control, and tropical
circuits whenever the number of shared parameters and witnesses is known.

### Theorem (finite-parametric response entropy)

For `w>=2`, `p,m>=1`, and `delta>0`,

```math
H_*:=\min\left\{2^{w+2p-2},{3^m-1\over2}\right\}. \tag{0}
```

Then

```math
\begin{aligned}
\log_2\operatorname {Cov}_{\delta}(C_{w,p,m})
\le{}&\log_2(m+1)+2m\log_2(w+p)\\
&+m\log_2(4(H_*+1))
+m\log_2\left(1+{w\over2\delta}\right).       \tag{1}
\end{aligned}
```

In particular, at macroscopic distortion `delta=epsilon w`,

```math
\log_2\operatorname {Cov}_{\epsilon w}(C_{w,p,m})
=O_\epsilon\bigl(m^2+m\log(w+p)\bigr).           \tag{2}
```

The bound is independent of the magnitudes and precision of all internal
edge weights.

Modulo constants, every private connected component which does not meet the
boundary can be deleted.  Every remaining connected component containing
private vertices has at least as many edges as private vertices.  Hence one
may take `p<=m`, and (2) becomes

```math
\log_2\operatorname {Cov}_{\epsilon w}(C_{w,m})
=O_\epsilon\bigl(m^2+m\log(w+m)\bigr).           \tag{3}
```

#### Proof

Fix first a labelled simple graph topology with `e<=m` edges and `p`
private vertices.  For a representative `s` of a projective boundary word
and a private spin assignment `z`, let

```math
a_{s,z}\in\{0,1\}^e
```

be the cut-incidence vector.  If `c in R_+^e` is the edge-weight vector,

```math
h_c(s)=\max_{z\in\{+1,-1\}^p}\langle a_{s,z},c\rangle. \tag{4}
```

This is where component size, rather than bit precision, enters.  Consider
all hyperplanes

```math
\langle a_{s,z}-a_{s,z'},c\rangle=0.             \tag{5}
```

There are at most

```math
H\le |X_w|\binom{2^p}{2}<2^{w+2p-2}              \tag{6}
```

listed comparisons.  There is a stronger bound which is important here:
both incidence vectors are binary, so every normal in (5) belongs to
`{-1,0,1}^e`, and opposite normals define the same hyperplane.  Consequently

```math
H\le {3^e-1\over2}\le {3^m-1\over2}.             \tag{6a}
```

This removes both the number of boundary queries and the number of private
spin assignments from the optimizer-cell count.  On each relatively open
face of the resulting arrangement, all
comparisons in (4) have fixed signs.  After a fixed lexicographic tie rule,
one private maximizer `z_s` is therefore fixed for every `s`; on that face
the entire response vector is a single linear map

```math
c\longmapsto (\langle a_{s,z_s},c\rangle)_{s\in X_w}. \tag{7}
```

An arrangement of `H` hyperplanes in `R^e` has at most

```math
\sum_{j=0}^e2^j\binom Hj\le [4(H+1)]^e           \tag{8}
```

relatively open faces.  (The first expression is the standard induction
bound for all faces, not only the full-dimensional chambers.)

Unit outer load implies that every response is one-Lipschitz for
`d_proj`: compare a fixed private cut before and after flipping boundary
coordinates, then compare maxima in both directions.  Since
`diam(X_w,d_proj)<=w/2`,

```math
\|[h_c]\|_{sh}={1\over2}\operatorname {osc}(h_c)\le {w\over4}. \tag{9}
```

For one arrangement face, (7) lies in a quotient linear subspace of
dimension at most `e` and in the radius-`w/4` ball of its induced norm.  The
usual volumetric argument in an arbitrary `e`-dimensional normed space gives
an internal `delta`-cover of size at most

```math
\left(1+{w\over2\delta}\right)^e.                \tag{10}
```

Indeed, take a maximal `delta`-separated set and compare the disjoint
`delta/2` balls with the ball of radius `w/4+delta/2`.

It remains to union over topologies.  Padding with isolated labelled
private vertices, every topology with at most `m` edges is encoded by at
most `m` slots, each containing either no edge or one of fewer than
`(w+p)^2` unordered pairs.  There are at most

```math
(m+1)(w+p)^{2m}                                   \tag{11}
```

such choices (a deliberately loose overcount).  Multiplying (8), (10), and
(11), using the minimum of (6) and (6a), and replacing `e` by `m` proves
(1).  In particular,

```math
\log_2(4(H_*+1))\le 2+m\log_2 3,
```

which proves (2).  The deletion and
spanning-forest observation proves `p<=m` in the projective response
problem. `square`

### Audit of the four delicate steps

* **Lower-dimensional arrangement faces.**  Counting only chambers would
  miss weight vectors with exact ties.  The total-face recurrence obtained
  when the `H`th hyperplane is inserted is bounded by
  `F(H,e)<=F(H-1,e)+2F(H-1,e-1)`.  With the evident base cases this gives
  `F(H,e)<=sum_(j<=e)2^j binom(H,j)`, the quantity used in (8).  A fixed
  tie-breaking optimizer is constant on every relatively open face,
  including zero-dimensional ones.
* **Nonconvex images and the quotient.**  No convexity of an arrangement
  face or of its response image is used.  In the quotient of the linear
  span of (7) by its constant directions, take a maximal separated subset
  of the image itself.  The translated half-radius norm balls are disjoint
  and all lie in one enlarged ball, so volume gives (10).  This also makes
  the cover internal to the response image if desired.
* **Unbounded internal weights.**  Parameter space may be unbounded and
  need not have finite volume.  Only its response image after quotienting
  constants is covered.  Unit outer load bounds that image by (9), even if
  internal weights tend to infinity or cancel at arbitrarily fine scales.
  Directions which change only the private optimum offset disappear in the
  quotient; other unbounded directions cannot occur because of (9).
* **Topologies and private labels.**  Parallel edges are first combined and
  zero edges deleted.  Labelling private vertices only overcounts.  Private
  components disjoint from the boundary are constant responses and are
  deleted before using `p<=m`.  Equation (11) then unions over every
  remaining topology, while the arrangement argument unions over every
  real weight vector on that topology.

### Corollary (robust landmark dimension)

Call `s_1,...,s_k` **gamma-independent response landmarks** for this class
if there are thresholds `t_i` such that, for every subset `U` of `[k]`, a
component `G_U` exists with its anchored response satisfying

```math
h_{G_U}(s_i)-h_{G_U}(s_0)\begin{cases}
\ge t_i+\gamma,&i\in U,\\
\le t_i-\gamma,&i\notin U.
\end{cases}                                       \tag{12}
```

Then

```math
k=O\left(m^2+m\log(w+p)+
m\log\left(1+{w\over\gamma}\right)\right).       \tag{13}
```

This is a more operational answer to “how many independent landmarks can
an `m`-edge graph expose?”  Every individual boundary coordinate is
query-exposable, but only the number in (13) can vary robustly and
independently across an `m`-edge family.

To prove it, the `2^k` response shapes in (12) form a `gamma`-packing: two
different label words differ by at least `2gamma` at one landmark while
their anchored difference is zero at `s_0`, so their shape distance is at
least `gamma`.  Apply (1) at any radius below `gamma/2`.

## Consequence for a universal unit-load compiler

The normalized distance-shell theorem in the current repository says that
unit outer load realizes the full class

```math
L_w=\operatorname {Lip}_1(X_w,d_{proj})/\mathbb R.
```

Its Hamming-code packing gives, for fixed `0<eta<1/4`,

```math
\log_2\operatorname {Cov}_{\eta w}(L_w)
\ge 2^{(1-H_2(2\eta)+o(1))w}.                    \tag{14}
```

Suppose `C_(w,m)` were an `epsilon w`-net of `L_w`.  For every fixed
`tau>0` with `epsilon+tau<1/4`, a `tau w`-cover of `C_(w,m)` would be an
`(epsilon+tau)w`-cover of `L_w`.  Combining (3) and (14) gives

```math
m^2+m\log(w+m)
\ge 2^{(1-H_2(2\epsilon+2\tau)+o(1))w}.          \tag{15}
```

In particular, first letting `w` grow and then `tau` decrease to zero,

```math
\liminf_{w\to\infty}{\log_2m\over w}
\ge {1-H_2(2\epsilon)\over2},
\qquad 0<\epsilon<1/4.                           \tag{16}
```

Thus polynomial-size components cannot approximate all unit-load responses
at a fixed macroscopic accuracy.  This conclusion still holds when every
edge weight is an arbitrary real carrying unlimited formal precision.  It
is strictly stronger than counting finite graph encodings.

For polynomial `p,m`, the response-bit entropy is polynomial while that of
the full unit-load ball is exponential in `w`.  Notice that arbitrary
numbers of private maximizers do not appear in (2): the finite cut-incidence
alphabet forces all their comparisons through only `3^m` possible normals.

## Why a bare `B`-bit theorem collapses to counting

There can be no stronger uniform conclusion from `B` alone unless the
decoder/grammar is restricted.  For every

```math
B\le \left\lfloor\log_2\operatorname {Pack}_{\epsilon w}(L_w)\right\rfloor,
```

choose any `2^B` members of an `epsilon w`-packing of `L_w`.  The normalized
distance-shell compiler realizes every one as a unit-load pure-Max-Cut
component (possibly of exponential size).  A hardwired decoder from
`B`-bit strings to these components therefore realizes exactly `B`
operational response bits.  The cardinality upper bound is attained.

For a conventional bounded-precision graph encoding with `b` bits per
weight, one of course also has

```math
\log_2|\{\hbox{encoded profiles}\}|
\le O(m\log(w+p)+mb),                             \tag{17}
```

but (17) is only representation counting.  The useful nontrivial resource
is the **max-affine presentation complexity** `(m,p)`: (1) bounds its
coarse response information even when `b=infinity`.  Equivalently, a useful
grammar must constrain the number of continuous edge parameters and private
optimization alternatives, rather than merely naming a bit budget.

## Interpretation and limitations

1. Unit boundary load controls regularity but not global complexity.  The
   private compiler plus one distance bridge can hide an exponential
   landscape behind a low-load interface.
2. Finite max-affine presentation controls global complexity.  A topology
   with `p` private spins presents each coordinate as a maximum of `2^p`
   linear forms in only `m` shared weights.  More strongly, the forms use
   binary cut-incidence coefficients, so every comparison belongs to a
   fixed ternary normal alphabet of size `3^m`.  Hyperplane cells plus
   finite-dimensional volume turn that fact into (1).
3. The estimate is probably not sharp in `m`.  Its `m^2` term is the generic
   arrangement ceiling for `3^m` ternary-normal hyperplanes.  Improving it
   would require a genuinely cut-specific bound on which of those cells can
   be robustly distinguished in response shape, not more description
   counting.
4. No computation was used: the distinction is theorem-level.  A useful
   next falsifier is to seek graph families whose weight space has
   `Omega(m^2)` robust optimizer-cell entropy.  Such a family would show that
   the arrangement term is real; its failure could motivate an improved
   near-linear law.

## Recommended next theorem

Determine the macroscopic number of *robustly separated* joint boundary
optimizer patterns of a fixed `m`-edge topology.  Concretely, prove one of

```math
\log N_{\rm robust\ cells}(\epsilon w)=O_\epsilon(m\log m)
```

or construct a unit-load family with

```math
\log N_{\rm robust\ cells}(\epsilon w)=\Omega(m^2).
```

The theorem above already shows that private-spin count is not an
independent information resource once the number of edges is fixed.  The
remaining question is whether the ternary oriented-matroid cell count is a
genuine robust response resource or only a zero-margin combinatorial
artifact.
