# Paired-selector switching: kernel partitions, cross differences, and observable cycles

**Status.** Independent proof-level report. No repository files are edited.

## Executive conclusions

Three points sharpen the current target.

1. The robust reset automaton can be much smaller than the current
   suffix-product power set.  The exact state is the **kernel partition** of
   the whole composite selector.  It has at most
   `Bell(r)-1` accepting values, and this is worst-case minimal for the full
   selector alphabet.
2. When the two compared branches use different selectors, diagonal error
   is not closed.  The forced carrier is the cross-difference array
   `D_(ij)=y_i-x_j`, on which the paired selector is
   `(i,j)->(tau(i),sigma(j))`.
3. Given an exact regular affine-selector presentation, a finite observable
   ordered-pair lift gives a necessary-and-sufficient coherent bounded-error
   theorem.  But a local face-adjacency graph need not be path-realizing:
   the translated clamp has a spurious identity-selector self-loop.

Thus a finite face/tangent skew product is complete only under an exact
Markov/path-lifting hypothesis.  Establishing such a finite quotient, not
adding more local face labels, is the unresolved step.

## 1. Mismatched selectors force cross differences

On a paired tie-free branch write

```math
x'_j=x_(sigma(j))+s_j,
qquad
y'_i=y_(tau(i))+t_i.                               \tag{1.1}
```

The diagonal error `e_i=y_i-x_i` does not close if `sigma != tau`:

```math
e'_i=y_(tau(i))-x_(sigma(i))+t_i-s_i.              \tag{1.2}
```

For an exact two-coordinate counterexample take zero translations,
`sigma=id`, `tau=(12)`, and `x=y=(0,z)`.  The old diagonal error is zero for
every `z`, while the new error is `(z,-z)`.  No state consisting only of the
old diagonal error and the discrete selector pair can determine the next
error.

Define instead

```math
D_(ij)=y_i-x_j,
qquad (i,j) in [r]^2.                              \tag{1.3}
```

Then

```math
\boxed{
D'_(ij)=D_(tau(i),sigma(j))+t_i-s_j.}              \tag{1.4}
```

This is an affine coordinate-selector update on `[r]^2`, with selector

```math
kappa(i,j)=(tau(i),sigma(j)).                      \tag{1.5}
```

The observed projective error is the half-oscillation of the diagonal
`D_(11),...,D_(rr)`.  Cross-difference arrays modulo constants are
equivalent to the pair `([x],[y])`, so their dimension is `2r-2`.  Indeed,
row differences recover `[y]` and column differences recover `[x]`.
Equation (1.2) shows why the `r-1` diagonal coordinates are generically
insufficient.

At a tie, each compatible pair `(tau,sigma)` gives a possible product
selector (1.5).  Which sequences of such choices are jointly feasible is a
separate symbolic-realizability problem.

## 2. The exact robust reset state is a kernel partition

Let a finite graph describe an exact regular selector language, and label an
edge by `sigma:[r]->[r]`.  If the old whole-path composite is `rho`, appending
`sigma` changes it to `rho circ sigma` under the convention

```math
P_sigma P_rho=P_(rho circ sigma).                  \tag{2.1}
```

For a transformation `rho`, let `ker rho` be its fibre partition.  If `Pi`
is a partition, define its pullback by

```math
i ~_(sigma^(-1)Pi) j
iff sigma(i) ~_Pi sigma(j).                        \tag{2.2}
```

### Theorem 2.1 (kernel-partition reset automaton)

Initialize every permitted path-start vertex with the discrete partition.
Store `(v,Pi)`, where `Pi` has at least two blocks, and update

```math
(v,Pi) -> (v',sigma^(-1)Pi).                      \tag{2.3}
```

Reject when the pullback has one block.  This automaton recognizes exactly
the reset-free paths and has at most

```math
\boxed{|V|(Bell(r)-1)+1}                           \tag{2.4}
```

states including one rejecting sink.

There are arbitrarily long reset-free paths exactly when the accepting lift
has a reachable cycle.  If its maximum accepting edge-height is `H`, every
legal path of length `H+1` contains a tangent reset, so the arbitrary-
residual bound is `(H+1)epsilon`.  A reachable cycle gives the linear lower
bound from Theorem 16.14.

#### Proof

Constant transformations form a two-sided ideal.  Hence a word contains a
constant-product factor if and only if its whole product is constant: one
direction uses the ideal property, and the other uses the whole word as the
factor.  Moreover

```math
ker(rho circ sigma)=sigma^(-1)(ker rho),           \tag{2.5}
```

so (2.3) stores exactly the only information needed to decide whether the
whole product is constant.  The cycle and height claims are ordinary finite
automaton pumping followed by Theorem 16.14. `square`

For the full selector alphabet this partition state is worst-case minimal.
If `Pi != Pi'`, choose `a,b` joined in exactly one partition and choose a
selector with image `{a,b}`.  Its pullback is the one-block rejecting
partition in one case and has two blocks in the other, so the states are
future-distinguishable.  A restricted graph alphabet may admit a further
Myhill--Nerode quotient.

This replaces the `2^(r^r-r)` suffix-set ceiling by a Bell-number ceiling.
The suffix set is correct but not minimal: the two-sided ideal property
already remembers every earlier reset.

For paired dynamics (1.4), the same theorem applies on `[r]^2` to the
restricted product-selector alphabet `tau times sigma`.  Its kernel
partition is not arbitrary:

```math
ker(tau times sigma)=ker(tau) times ker(sigma).     \tag{2.6}
```

Hence it is enough to store the pair of channel kernel partitions, at most
`Bell(r)^2-1` accepting values, rather than a general partition of `r^2`
coordinates.  This is the full-cross-carrier reset state; a smaller
observation-specific reset quotient would require using the declared
diagonal response relation.

## 3. A local face cycle need not pump

Consider

```math
T_delta(z)=clip(z+delta,0,1),
qquad0<delta<1/2.                                  \tag{3.1}
```

This is induced by the all-finite translated-clamp matrix already verified
in the project.  Its middle face `I=(0,1-delta)` has identity selector, and

```math
T_delta(I) cap I=(delta,1-delta) != empty.         \tag{3.2}
```

The usual local face graph therefore contains an identity-selector
self-loop.  Yet

```math
T_delta^k(z)=min(z+k delta,1),                     \tag{3.3}
```

so no orbit remains in `I` indefinitely; the number of consecutive middle-
face steps is at most `ceil(1/delta)`.

Thus

```math
reachable local face cycle => pumpable selector cycle
```

is false.  An exact presentation for this example can subdivide the middle
face into `Theta(1/delta)` progress layers.  The ordinary max-affine normal
fan does not contain that quantitative state.  Listing every locally
possible tie resolution has the same over-approximation problem.

## 4. Exact coherent theorem under a regular presentation

Let a finite directed graph have edges `e:u->v` labelled by affine selector
maps

```math
A_e z=P_(sigma_e)z+b_e,
qquad z in R^q/R1.                                 \tag{4.1}
```

At a terminal vertex `v`, let `O_v` be the ordered coordinate pairs whose
differences are observed.  For full Hilbert response take all pairs.  For
the paired carrier (1.4), take pairs of diagonal indices

```math
((i,i),(j,j)),
qquad i != j.                                     \tag{4.2}
```

Build the backward observable pair lift.  For a base edge `e:u->v`, add

```math
(v,i,j) -> (u,sigma_e(i),sigma_e(j))              \tag{4.3}
```

with weight

```math
w_e(i,j)=b_(e,i)-b_(e,j),                         \tag{4.4}
```

provided the two target coordinates differ.  If they agree, the backward
path terminates: all earlier input differences have been reset.  Retain only
vertices reachable backward from a declared observation pair.

### Theorem 4.1 (observable zero-cycle criterion)

Assume every finite graph path is an exact feasible affine-selector
itinerary and initial coordinate oscillations are uniformly bounded.

* All declared projective response differences are uniformly bounded over
  path length if and only if every directed cycle of the reachable backward
  pair lift has total weight zero.
* If the reachable lift has `N` vertices and every edge weight has absolute
  value at most `W`, the accumulated translation difference is at most
  `NW`, plus the initial oscillation.
* A reachable cycle of weight `c != 0` yields a legal repeated word with
  response drift at least `k|c|/2-O(1)` after `k` repetitions.

#### Proof

One edge obeys

```math
(A_ez)_i-(A_ez)_j
=z_(sigma_e(i))-z_(sigma_e(j))+w_e(i,j).          \tag{4.5}
```

Iterating backward expresses a terminal observed difference as the weight
of one lifted path plus an initial coordinate difference if the pair never
collapses.  If all lifted cycles have weight zero, delete closed subwalks
without changing the weight.  What remains is simple, has at most `N-1`
edges, and a possible terminating collapse edge raises the safe bound to
`NW`.  The initial term is bounded by hypothesis.

Conversely, traverse a reachable nonzero lifted cycle `k` times, reverse its
base-edge sequence, and append the fixed access word to the terminal
observation.  Exact path realization makes this a legal forward itinerary.
Equation (4.5) gives `kc+O(1)`, and half-oscillation gives the displayed
projective drift. `square`

Equivalently, the lifted weight is a scalar coboundary on every reachable
strongly connected component.  For one base vertex and one selector edge,
the theorem is exactly the current common-cycle-mean/twisted-coboundary
criterion.

## 5. Mixed switching cycles are essential

Let `R` swap two coordinates and define

```math
A(z)=P_Rz,
qquad B(z)=P_Rz+(0,1).                             \tag{5.1}
```

For the projective coordinate `t=z_2-z_1`,

```math
A(t)=-t,
qquad B(t)=-t+1.                                  \tag{5.2}
```

Individually,

```math
A^2=B^2=id,                                       \tag{5.3}
```

so each generator passes its one-map cycle-mean test.  But

```math
(B circ A)(t)=t+1,                                \tag{5.4}
```

and the mixed word drifts linearly.  The backward pair lift has a mixed
two-edge cycle of weight one, while both monochromatic two-edge cycles have
weight zero.

These are sparse monomial max-plus maps.  For any fixed finite horizon they
can be embedded on the visited bounded region into all-finite maps by making
inactive entries sufficiently negative, but this does **not** produce one
all-finite system with infinite drift: a fixed all-finite max-plus map has
bounded projective image.  The example is therefore an exact weighted-
automaton theorem and a local/tangent switching falsifier, not a global
all-finite coherent-drift example.  It still proves that checking each
continuation separately cannot establish switching stability; all
transported mixed cycles must be tested whenever the exact affine-selector
language permits them.

## 6. Ties and the exact boundary of the finite theorem

On a tie face, represent every compatible argmax selector by a parallel
edge.  Theorems 2.1 and 4.1 are necessary and sufficient only if this graph
has an exact path-lifting property:

1. every actual paired/directional selector itinerary maps to a graph path;
2. every graph path used as a cycle certificate is realized by an actual
trajectory or directional perturbation.

Soundness alone makes the zero-cycle condition a valid sufficient
certificate, but a graph nonzero cycle can be spurious without the second
condition.  The translated clamp is the finite counterexample.

For paired max-plus systems the discrete label is

```math
(paired face, tau times sigma),                   \tag{6.1}
```

acting on the continuous cross-difference carrier.  A finite exact Markov
partition of that carrier would yield a necessary-and-sufficient bounded-
error theorem by Theorem 4.1.  No such partition follows merely from having
finitely many original max-affine faces.

The zero-cycle test itself is weighted-automata / graph-potential folklore.
The response-theoretic content is the forced cross-difference carrier, the
restriction to observable diagonal pairs, reset as termination of a
backward witness, the conversion of a nonzero accessible cycle into Hilbert
response drift, and the Bell-number kernel-partition reduction for robust
residuals.

## 7. Judgment

A finite face/tangent skew product does give a complete theorem **when an
exact regular presentation is supplied**.  The naive face graph does not
supply one.  The strongest next theorem is therefore a structural Markov-
refinement result for a natural all-finite max-plus subclass, not a more
elaborate local cycle criterion.
