# Wave 32 Route 1: exact completion-forest min-cuts

## 1. Conclusion

For a fixed rooted tree of partially specified projective words and a fixed
auxiliary full center, the minimum total projective-Hamming length of full
completions has an exact formula.  One first chooses **one projective
orientation per word**, shared by all its specified coordinates.  Conditional
on those orientations, the free outside bits decouple coordinate by
coordinate.  Coordinate `i` is the unit-capacity tree min-cut separating the
pinned `+` terminals from the pinned `-` terminals.  Equivalently, its cost is
the maximum number of edge-disjoint tree paths joining oppositely pinned
terminals.

Thus the exact label statistic is an orientation-consistent sum of
coordinate conflict-path packings.  Minimizing it over favorable labels,
low-row centers, trees, and a partition into at most `n^eta` groups is
**necessary and sufficient**, not merely sufficient, for the metric part of
the forest event in (10.910).  This identifies the weakest local agreement
input sought by the route.

The principal-ground exchange cocycle does not imply this statistic as a
bare algebraic consequence.  It controls signed first-order energy totals,
whereas the min-cut counts unweighted, orientation-consistent terminal
conflicts, including conflicts separated by arbitrarily many unspecified
tree nodes.  Replicated versions of the weighted zero-shore gadget from Wave
31 have arbitrarily many exact grounds with zero mutual energy regret and
linear projective separation.  This is only a scoped weighted/zero-edge wall,
not a complete-signing counterexample.  The precise missing input is a
minimizer-specific higher-order **cut-compression theorem** producing small
coordinate separator sets (or ruling out the dual conflict-path packing) on
a likely batch.

## 2. Setup and projective lifting

Let `T=(V,E)` be a tree rooted at `0`, with non-root vertices
`V^*=V\{0}`.  Fix a center `z in {+1,-1}^n`.  At vertex `v in V^*`, let
`S_v subseteq [n]` and let `[y^v]` be a projective word on `S_v`; thus
`y^v` and `-y^v` represent the same prescribed child cut.  Its projective
completion cylinder is

```math
\mathcal C_v
=\{[x]:x\in\{\pm1\}^n,\ x|_{S_v}=\pm y^v\}.
```

Define the best rooted completion-tree length

```math
D_T(z;\{(S_v,[y^v])\})
=\min_{[x^v]\in\mathcal C_v}
 \sum_{uv\in E}d_{\rm pr}([x^u],[x^v]),
\qquad [x^0]=[z].                                      \tag{R32.1}
```

Here `d_pr([x],[x'])=min(d_H(x,x'),n-d_H(x,x'))`.

**Projective lifting lemma.**  Fix the representative `x^0=z`.  For every
choice of projective completions in (R32.1), orient representatives
recursively away from the root so that every edge realizes its projective
distance as ordinary Hamming distance.  This is possible simultaneously on
all edges precisely because `T` has no cycle.  Conversely, any such oriented
representatives give feasible projective completions.  Consequently

```math
D_T
=\min_{\substack{\sigma_v\in\{\pm1\}\ (v\ne0)\\
                  x^0=z,\ x^v|_{S_v}=\sigma_vy^v}}
  \sum_{uv\in E}d_H(x^u,x^v).                          \tag{R32.2}
```

The sign `sigma_v` is a **node sign**, not a coordinate sign.  It must be the
same on every coordinate of `S_v`.  Changing the displayed representative
`y^v` to `-y^v` merely replaces `sigma_v` by `-sigma_v`, so the formula is
projectively invariant.

## 3. Exact coordinate min-cut formula

Put `S_0=[n]`, `y^0=z`, and `sigma_0=1`.  For a fixed node-orientation vector
`sigma`, coordinate `i` has terminal set and prescribed signs

```math
P_i=\{v:i\in S_v\},\qquad b_{v,i}^{\sigma}=\sigma_vy_i^v
\quad(v\in P_i).                                       \tag{R32.3}
```

The root belongs to every `P_i`.  Define

```math
\lambda_i(\sigma)
=\min_{g:V\to\{\pm1\}}
 \left\{\sum_{uv\in E}{\bf1}_{g(u)\ne g(v)}:
                 g(v)=b_{v,i}^{\sigma}\ (v\in P_i)\right\}.   \tag{R32.4}
```

**Theorem (exact completion-forest formula).**

```math
\boxed{
D_T(z;\{(S_v,[y^v])\})
=\min_{\sigma\in\{\pm1\}^{V^*}}\sum_{i=1}^n\lambda_i(\sigma).
}                                                       \tag{R32.5}
```

*Proof.*  Expand the Hamming sum in (R32.2) over coordinates.  Once `sigma`
is fixed, the only remaining variables are the free values `x_i^v` with
`i notin S_v`; variables of different coordinates occur in disjoint
summands and constraints.  The minimum for coordinate `i` is exactly
(R32.4).  Finally minimize over the shared node orientations.  `square`

Only the outside bits decouple.  Moving `min_sigma` inside the coordinate
sum is generally false.  With one non-root node, center `z=(+,+)`, and
specified word `y=(+,-)` on both coordinates, a separate orientation per
coordinate would report zero, while the correct projective distance is one.

## 4. Cut and dual path-packing characterizations

For fixed `sigma,i`, let `P_i^+` and `P_i^-` be the positive and negative
terminals.  Then

```math
\lambda_i(\sigma)
=\min\{|C|:C\subseteq E,\ 
        \text{no component of }T-C\text{ meets both }P_i^+\text{ and }P_i^-\}.
                                                               \tag{R32.6}
```

Indeed, the disagreement edges of any feasible binary labeling form such a
separator.  Conversely, after deleting a separating set, each component has
terminals of at most one sign; label it by that sign (and label terminal-free
components arbitrarily).  This uses at most the deleted edges.

Equivalently, attach a source to every terminal in `P_i^+` and a sink to
every terminal in `P_i^-` by capacity larger than `|E|`, while tree edges have
unit capacity.  Max-flow/min-cut and integrality give

```math
\boxed{
\lambda_i(\sigma)
=\max\{|\mathcal P|:\mathcal P\text{ is a family of edge-disjoint paths in }T,
\text{ each joining }P_i^+\text{ to }P_i^-\}.
}                                                       \tag{R32.7}
```

The corresponding path LP is

```math
\lambda_i(\sigma)
=\max_{f_{uv}\ge0}
 \left\{\sum_{u\in P_i^+,v\in P_i^-}f_{uv}:
 \sum_{u,v:e\in [u,v]_T}f_{uv}\le1\quad(e\in E)\right\},    \tag{R32.8}
```

and it has an integral optimum.  Combining the coordinates, (R32.5) is

```math
D_T=\min_{\sigma}
 \max\left\{\sum_{i,u,v}f^i_{uv}:
 f^i\text{ obeys (R32.8) for each }i\right\}.           \tag{R32.9}
```

The inner adversary packs coordinate-tagged conflict paths, with capacity
one for every pair `(coordinate,tree edge)`.  There is no justified exchange
of `min_sigma` and `max_f`.

A dynamic program gives another exact implementation.  Root `T` and let
`M_{v,i}(a)` be the minimum cost in the descendant subtree of `v` conditional
on label `a` at `v`.  Set it to infinity when a terminal pin is violated and
otherwise use

```math
M_{v,i}(a)=\sum_{u:\,u\text{ child of }v}
             \min_{c\in\{\pm1\}}
             \{M_{u,i}(c)+{\bf1}_{a\ne c}\}.            \tag{R32.10}
```

Then `lambda_i=M_{0,i}(z_i)`.

## 5. The exact local agreement statistic

For fixed `T,z` define the orientation-consistent tree-conflict statistic

```math
\mathsf{Conf}_T(z;\mathbf y)
:=\min_\sigma\sum_i\lambda_i(\sigma).                  \tag{R32.11}
```

By (R32.5), `Conf_T=D_T`.  It is therefore the weakest exact statistic of
the partial labels for the metric conclusion: `Conf_T=O(k_0)` is necessary
and sufficient for full completions of total rooted-tree length `O(k_0)`.
In primal language it asks for shared node orientations and coordinate edge
sets `C_i` of total size `O(k_0)` which separate every pair of oppositely
pinned terminals.  In dual language it says that, after one coherent choice
of node orientations, the total maximum conflict-path packing is `O(k_0)`.

For comparison, let `a_i(v)` be the nearest strict ancestor of a terminal
`v in P_i\{0}` which also belongs to `P_i`.  The compressed-terminal
alternation count

```math
\mathsf{Alt}_T(\sigma)
=\sum_i\sum_{v\in P_i\setminus\{0\}}
 {\bf1}_{b_{v,i}^\sigma\ne b_{a_i(v),i}^\sigma}        \tag{R32.12}
```

satisfies `sum_i lambda_i(sigma)<=Alt_T(sigma)`: for every disagreeing
terminal `v`, cut the last original-tree edge entering `v`.  Along the
nearest-terminal-ancestor chain between any two opposite terminals there is
a first sign change, and that entering edge is cut.  This convenient
pair/ancestor
certificate can overcount arbitrarily because many conflicts may be stopped
by one shared edge.  At the other extreme, agreement only on overlaps of
adjacent selector nodes is too weak: on a three-node path, a coordinate can
be pinned `+` at the root, unspecified at the middle node, and pinned `-` at
the leaf.  Every adjacent overlap check is vacuous, but its min-cut cost is
one.  Replication makes the gap arbitrary.

To state exactly the Wave 31 batch event, let `Y_a` be the allowed favorable
projective words for sample `a`, and let

```math
\mathcal Z_C(A)=\{[z]:R_2(z)\le Cn^2\}.
```

For an integer `G`, minimize over: a partition into `ell<=G` nonempty groups;
one `y^a in Y_a` per sample; one center in `Z_C(A)` per group; and one rooted
tree on that group plus its center.  Define

```math
\mathsf{GroupConf}_G(\mathbf S)
=\min \max_{1\le g\le\ell}\mathsf{Conf}_{T_g}(z_g;\mathbf y_g).
                                                               \tag{R32.13}
```

Then, with `G=floor(n^eta)`,

```math
\boxed{\mathsf{GroupConf}_G(\mathbf S)=O(k_0)}             \tag{R32.14}
```

is exactly necessary and sufficient for the metric and low-row-center part
of the structured forest event after (10.911).  On the adversarial-law event
of probability `exp{-O(rL_0)}`, (R32.14), favorable deficits at the recorded
tolerance, and `eta<c_0` imply the `T<=n^eta` row-good coset cover and hence
convergence by (10.906).  Thus (R32.14) is the sharp local-label successor;
it does not by itself prove that the event has the required probability.

Its exact falsification certificate is also clear: after every permitted
partition, favorable-label choice, low-row center, and tree choice, some
group has the property that **every** coherent node orientation admits more
than `Ck_0` coordinate-tagged, edge-disjoint conflict paths in total.

## 6. Audit against principal-ground exchanges

The exchange identity (10.913) controls signed linear boundary-energy
totals.  Neither its terms nor its nonnegativity distinguish how many
coordinates carry incompatible prescribed spins, and cancellation can make
all exchange totals zero.  The min-cut statistic instead needs a coherent
orientation and a small separator for every coordinate terminal system.
It is both unweighted and higher-order: a conflict can run through a chain of
nodes on which that coordinate is unspecified, and several branch conflicts
may share one cut edge.

There is a direct scoped algebraic wall.  In the five-vertex real weighted
gadget from (10.915), the all-positive spin and the spin obtained by flipping
`F={0,1,2}` are both exact grounds and their signed shore exchange is zero,
but their projective Hamming distance is two.  Take a block-diagonal sum of
`t` copies (so cross-block weights are zero).  Each block has an independent
zero-cost relative sign, giving `2^t` exact grounds and zero energy regret
between every pair.  Standard binary codes select exponentially many of
these grounds with pairwise projective distance `Omega(t)`.  Any completion
tree containing two such witnesses has total length at least their metric
distance, regardless of the auxiliary center.  Consequently, fewer groups
than selected separated witnesses force a group of length `Omega(t)`, while
all bare exchange regrets remain zero.

This does **not** falsify the project route: the block matrix has real and
zero weights, is not a complete signing, does not model the fixed-density
principal deck, and need not be an exact minimizer among complete signings.
It does prove that (10.913) alone cannot yield (R32.14) by algebraic
rearrangement or positivity.

The precise missing higher-order input can be stated without slack:

> Uniformly against every selector law, a batch with probability
> `exp{-O(rL_0)}` admits favorable child words and a partition into at most
> `n^eta` groups such that each group has a low-row center, a selector tree,
> and one coherent node orientation for which coordinate separators of total
> size `O(k_0)` exist.

Equivalently, it must rule out the dual conflict-path certificate following
(R32.14).  To derive this from signing minimality one needs a new coercive
principle which charges each disjoint terminal-conflict path to an
uncancelled resource of total `O(k_0)`.  The known cross-regret cycles provide
only signed shore sums and have no positive cost per packed path.  A
quadratic boundary-susceptibility inequality, a Hamming stability theorem for
principal grounds, or a true higher-order agreement theorem could supply the
charge; none is currently proved.

## 7. Verification

`tmp/forest_mincut_r32_check.py` independently checks on exhaustive/random
small instances that:

1. brute-force projective completion agrees with (R32.5);
2. the coordinate dynamic program agrees with the separator and integral
   conflict-path-packing formulations;
3. coordinatewise node orientations give the advertised false zero in the
   two-coordinate example; and
4. the five-vertex weighted gadget has the stated two ground classes, zero
   shore exchange, and positive projective separation.
