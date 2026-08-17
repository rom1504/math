# Exact selector presentations and observable witness cycles

**Status.** Proof source for Theorem 16.19 and Corollary 16.16a. Tiny exact
counterexamples are checked by
[`../experiments/verify_paired_selector_cycles.py`](../experiments/verify_paired_selector_cycles.py).

## 1. The minimal robust reset state

For selectors `rho,sigma:[r]->[r]`, use the convention that appending
`sigma` changes a path product `rho` to `rho circ sigma`. Constants form a
two-sided ideal. Hence a word contains a constant-product factor if and only
if its whole product is constant.

It follows that the future reset state is not the set of all suffix products.
It is the fibre partition `ker rho`. If `Pi` is a partition, update it by

```math
i\sim_{\sigma^{-1}\Pi}j
\quad\Longleftrightarrow\quad
\sigma(i)\sim_\Pi\sigma(j).                        \tag{1.1}
```

The one-block partition is the reset sink. Beyond a fixed control vertex this
gives `Bell(r)-1` accepting states, and is worst-case minimal for the full
selector alphabet: if two partitions differ, choose `a,b` joined in exactly
one and a continuation whose image is `{a,b}`. One state resets and the other
does not. This does not assert that control vertices of a restricted regular
language cannot merge.

For a paired selector `(tau,sigma)` on cross coordinates `(i,j)`,

```math
ker(tau times sigma)=ker(tau) times ker(sigma),     \tag{1.2}
```

so a pair of kernel partitions suffices; an arbitrary partition of `r^2`
points is unnecessary. For the full independent product-selector alphabet
there are `Bell(r)^2-1` accepting pairs plus one sink; diagonal-observation
or restricted languages may quotient further.

## 2. Exact regular affine-selector presentations

Let a finite directed multigraph have permitted starts and edges

```math
e:q->q',
\qquad A_ez=P_(sigma_e)z+b_e.                      \tag{2.1}
```

At a terminal control state `q`, declare ordered coordinate observations
`O_q`. The presentation is **sound** if every actual resolved orbit gives a
graph path with the displayed affine updates. It is **path-realizing** if
every graph path used by the presentation, including repetitions of a cycle
with fixed access paths, is realized by an allowed finite orbit. A finite
invariant cell partition is path-realizing when every map has one
selector-affine formula on a whole cell and maps it into one declared
successor; tie faces require a common realizable tangent refinement.
Nonempty one-step intersections are not enough.

Build the reverse witness graph with vertices `(q,i,j)`. For each edge
`e:q->q'`, add

```math
(q',i,j)->(q,sigma_e(i),sigma_e(j))                \tag{2.2}
```

of weight `b_e(i)-b_e(j)`. Retain vertices on a reverse path from a declared
terminal observation to an allowed start. Assume allowed initial coordinate
oscillation is at most `R_0`.

### Theorem 2 (observable witness-cycle dichotomy)

For a sound presentation, absence of positive relevant cycles implies every
declared directed output is at most `R_0+K`, where `K>=0` is the maximum
weight of a simple relevant reverse path, including the empty path. If the
presentation is path-realizing, a
relevant positive cycle of weight `c` gives fixed access, cycle, and exit
words `u,v,w` such that, for every `k`, the supremum over allowed finite-orbit
seeds of the output on `uv^kw` is at least `kc-C`. Thus absence of positive
cycles is necessary and sufficient for uniform directed upper boundedness.
The seed may depend on `k`; one infinite pumping orbit requires
nested-cylinder realization.

Under the same path-realization assumption, if observations and relevance
are closed under pair reversal, two-sided projective boundedness is
equivalent to zero weight on every relevant cycle. Equivalently, edge weights
are a vertex coboundary on each relevant strongly connected component.

#### Proof

One edge satisfies

```math
(A_ez)_i-(A_ez)_j
=z_(sigma_e(i))-z_(sigma_e(j))+b_e(i)-b_e(j).       \tag{2.3}
```

Iteration expresses a terminal difference as one bounded initial difference
plus the weight of its reverse witness path. If cycles are nonpositive,
deleting them cannot lower the path weight; a maximizing path is simple.
A positive relevant cycle can be repeated under path realization and (2.3)
then grows by `c` each time. Pair reversal negates all weights, so two-sided
boundedness forces every cycle to have weight zero. The usual path-independence
argument identifies zero cycle weights with a coboundary. `square`

## 3. Paired channels and the information boundary

For different selectors

```math
x'_j=x_(sigma(j))+s_j,
\qquad y'_i=y_(tau(i))+t_i,
```

the diagonal error does not close. The forced joint carrier is

```math
D_(ij)=y_i-x_j,
\qquad
D'_(ij)=D_(tau(i),sigma(j))+t_i-s_j.               \tag{3.1}
```

It is an affine selector system on `r^2` coordinates, while its physical
dimension is only `2r-2`. The observed Hilbert error uses witness pairs
`((i,i),(j,j))`. Theorem 2 therefore needs at most `|Q|r^4` witness states,
not the full trajectory or every response value. It preserves joint
cancellation before absolute values. The endpoint hypothesis concerns the
full cross carrier: `x=y` gives zero diagonal error but unbounded off-diagonal
entries when the common projective input is unrestricted. A bounded-image
prefix, reset, or bounded initial domain is needed before differing selectors
may expose those entries.

The exact-language hypothesis is real. The all-finite projective map

```math
z->clip(z-delta,0,1)
```

has a slope-one middle cell with a nonempty one-step self-intersection, but
every orbit exits after finitely many steps. Its local face graph therefore
has a spurious pumpable cycle. At ties, independently branching argmax
selectors can likewise create paths no common perturbation realizes.

The theorem completely solves the response cocycle once an exact finite
presentation is supplied. It does not prove that such a presentation exists
or is small. The remaining structural question is tropical lumpability: find
natural invariant finite partitions whose whole cells map to whole cells.
