# Multitype synchronized dense bridges

Status: main-agent theorem draft.  The finite identity is accompanied by an
exact exhaustive verifier.  This extends the one-orbit `alpha I+beta J`
example without changing its essential common-section mechanism.

## Model

Partition a common coordinate set into nonempty types

```math
V=C_1\sqcup\cdots\sqcup C_K,\qquad |C_c|=n_c,
```

and let block `a` carry a spin vector `x^a in {-1,1}^V`.  Write

```math
k_{a,c}=|\{i\in C_c:x_i^a=1\}|,
\qquad s_{a,c}=2k_{a,c}-n_c.                                \tag{MT.1}
```

The internal energy `h_a` may be any function of the complete count vector
`k_a=(k_(a,1),...,k_(a,K))`.  Couple blocks `a<b` by

```math
R_{ab}=\sum_c alpha_(ab,c)D_c
       +\sum_(c,d) beta_(ab,cd)1_(C_c)1_(C_d)^T,             \tag{MT.2}
```

where `D_c` is the diagonal projector onto `C_c` and every
`alpha_(ab,c)>=0`.  The coefficients `beta` are arbitrary real numbers.

### Theorem MT.1 (exact multitype common-section quotient)

The microscopic optimum is exactly

```math
max_(k_(a,c))\left\{
 \sum_a h_a(k_a)
 +\sum_(a<b)\left[
  \sum_c alpha_(ab,c)(n_c-2|k_(a,c)-k_(b,c)|)
  +\sum_(c,d)beta_(ab,cd)s_(a,c)s_(b,d)
 \right]\right\}.                                         \tag{MT.3}
```

For every fixed count array, one common microscopic representative attains
all of the pairwise bounds in (MT.3) simultaneously.

#### Proof

Inside type `c`, let `P_(a,c)` be the plus set of block `a`.  Then

```math
\sum_(i\in C_c)x_i^ax_i^b
=n_c-2|P_(a,c)\triangle P_(b,c)|
\le n_c-2|k_(a,c)-k_(b,c)|.                                \tag{MT.4}
```

The nonnegative coefficients make the sum of these inequalities an upper
bound.  Give every type an ordering and choose, simultaneously for every
block,

```math
P_(a,c)=\{\text{the first }k_(a,c)\text{ coordinates of }C_c\}. \tag{MT.5}
```

Within each type these sets are nested, so every inequality (MT.4) is an
equality.  The block-constant term in (MT.2) is already
`sum_(c,d)beta_(ab,cd)s_(a,c)s_(b,d)` and depends only on the
counts.  This proves (MT.3). `square`

The exact **joint search grid** has

```math
\prod_c(n_c+1)^m
```

labels, hence at most `mK log_2(n+1)` label bits when `n=sum_c n_c`.
The matrices can be dense and full rank: block-constant terms make them
dense, while nonzero diagonal coefficients remove the large blockwise
kernels generically.  Thus algebraic rank is not the controlling resource.
For `K=o(n/log n)` and fixed `m`, this joint optimization label is
subextensive in the `mn` microscopic spin bits.  This is not, by itself, an
operational state for a fragment whose already chosen microscopic alignment
is exposed to an arbitrary future.  Such a rooted future can distinguish two
alignments with the same counts.  The quotient is reusable only in the
declared class where every old and new block is jointly reoptimized and all
identity channels share the common-section hypotheses.

### Theorem MT.2 (signed typewise balance)

Allow nonzero `alpha_(ab,c)` of either sign.  Balance of the signed block graph
is sufficient at every type size.  When `n_c>=2`, it is also necessary for
the edgewise overlap optima to admit one common representative section
**uniformly for every count assignment**.  Equivalently, the exact condition
is that the signed block graph is balanced, or equivalently

```math
sgn(alpha_(ab,c))=epsilon_(a,c)epsilon_(b,c)                 \tag{MT.6}
```

for typewise vertex gauges `epsilon_(a,c)`.  If this holds for every type,
flip block `a` on the whole cell `C_c` according to `epsilon_(a,c)` and
apply Theorem MT.1, with the induced changes to `h_a` and `beta`.

If one type has even size and contains an isolated unbalanced unit-sign cycle
whose blocks are pinned to zero type-magnetization, separately optimizing the
pair responses overestimates the true contribution by exactly `2n_c`.

#### Proof

The balance criterion and gauge construction apply independently on each
type.  After the gauge, (MT.5) supplies a common section.  Conversely, delete
one edge from an unbalanced cycle and choose signs `epsilon_a` along the
remaining path so that every path-edge sign is
`epsilon_a epsilon_b`.  Assign count one where `epsilon_a=1` and count
`n_c-1` where `epsilon_a=-1`.  Along a positive path edge, simultaneous
pairwise optimality forces the equal one-sets (or equal co-one-sets); along a
negative path edge it forces a one-set and its complement.  Hence the path
forces one fixed singleton `P` at every `+` vertex and `P^c` at every `-`
vertex.  The deleted edge has sign different from the product of its endpoint
gauges.  If it is positive, the forced endpoint sets are complements rather
than maximally aligned; if it is negative, they are equal rather than
maximally anti-aligned.  It cannot attain its pair optimum.  Thus no common
section works for all count assignments.

For the quantitative even-size witness, the product of realized edge signs
is `+1` at every coordinate while the desired product is `-1`; at least one
unit edge is missed, costing two.  Pair a maximizing coordinate pattern with
its global negative to satisfy the zero-magnetization constraint and attain
total loss `2n_c`. `square`

### Theorem MT.3 (fixed-type thermodynamic limit)

Fix `m,K`.  Suppose `n_c/n -> p_c>0` and the typewise signed graphs are
balanced.  Perform the gauges (MT.6) first.  In these gauged coordinates,
suppose the diagonal coefficients converge to `a_(ab,c)>=0`, the transformed
block-constant coefficients are `b_(ab,cd)/n`, and the transformed internal
landscapes satisfy, uniformly,

```math
n^(-1)h_(a,n)(x^a)
 -> f_a((u_(a,c))_(c<=K)),
\qquad u_(a,c)=s_(a,c)/n_c,                                 \tag{MT.7}
```

for continuous `f_a` on `[-1,1]^K`.  Here `u_(a,c)` is the **gauged**
type-magnetization; explicitly the transformation sends
`u_(a,c)` to `epsilon_(a,c)u_(a,c)` and
`b_(ab,cd)` to `epsilon_(a,c)epsilon_(b,d)b_(ab,cd)` before the displayed
assumptions are read.  Then the normalized optimum converges to the maximum
over `u in [-1,1]^(mK)` of

```math
\sum_a f_a(u_a)
+\sum_(a<b)\left[
 \sum_c a_(ab,c)p_c(1-|u_(a,c)-u_(b,c)|)
 +\sum_(c,d)b_(ab,cd)p_cp_du_(a,c)u_(b,d)
\right].                                                    \tag{MT.8}
```

#### Proof

Theorem MT.1 reduces the finite optimum to the product of the typewise
parity grids.  After division by `n`, every term converges uniformly to
(MT.8); those grids become dense in the compact cube. `square`

## Interpretation

The reduced joint optimization is predicted by orbit structure plus an
optimizer-compatible common section, not by matrix rank.  Each new coordinate
type adds one aggregate count and one independent signed-holonomy audit.  The
growth law `O(mK log n)` is therefore a concrete intermediate regime for the
declared jointly reoptimizable family.  It is not an arbitrary-future
fragment quotient.  The theorem ceases to close when edge pairs use
incompatible coordinate partitions, when their signed type graphs are
unbalanced, or when a future freezes an alignment chosen earlier; the missing
information is then joint overlap geometry rather than another marginal
count.
