# Query-local Walsh amalgamation

Status: rigorous task-local draft.  The upper theorem is an exact response
quotient for a declared family of graph futures.  The lower theorem is a
target-scale semantic packing, not merely an orbit count.  An exact finite
verifier accompanies the draft.

## 1. Declared graph futures

Work over `F_2`.  Let `V=F_2^m` have its standard bilinear form and let
`omega=(1,...,1)`.  Two marked linear-label pieces are presented by

```math
alpha:E_a=F_2^r\longrightarrow V,
\qquad
beta:E_b=F_2^s\longrightarrow V.                       \tag{QL.1}
```

Write `R_a=ker(alpha)`, `R_b=ker(beta)`,

```math
U_a=E_a/R_a,
\qquad U_b=E_b/R_b,                                    \tag{QL.2}
```

and identify these quotient spaces with their images in `V`.  The relative
rooted amalgamation datum is

```math
\begin{aligned}
kappa(u,v)&=u\mathbin\cdot v,\\
J&=\{(u,v)\in U_a\oplus U_b:u=v\},\\
Z^times&=\{(u,v)\in U_a\oplus U_b:u+v=omega\}.
\end{aligned}                                           \tag{QL.3}
```

A **Walsh graph query** is a real weighted graph on the marked blocks.  A
vertex carries a scalar multiple of its usual linear-label Walsh child, and
an edge carries a scalar multiple of the common Walsh bridge.  A block is
active iff it occurs in a nonzero child term, a nonzero incident bridge, or
an attached rooted term.  We also allow independent componentwise rooted
continuations which transform equivariantly under that component's Walsh
coordinate relabelling.  No external spin or field may be shared between two
nominal components; such a term joins their supports.

Fix a finite declared family `Theta` of such queries.  For a query `theta`,
discard only blocks absent from every term and let its **connected supports**
be the active-vertex sets in the connected components of its nonzero bridge
graph.  An active vertex with no bridge counts as a singleton support.  Let
`C_Theta` be the inclusion-maximal supports occurring in any declared query.

For `C in C_Theta`, let

```math
P_C=alpha(span\{e_i:i\in C\cap[r]\})\subseteq U_a,
\qquad
Q_C=beta(span\{e_j:j\in C\cap[s]\})\subseteq U_b.       \tag{QL.4}
```

Define the **query-local relative carrier**

```math
D_Theta(alpha,beta)=
\left(
 kappa|_(P_C\times Q_C),
 J\cap(P_C\oplus Q_C),
 Z^times\cap(P_C\oplus Q_C)
\right)_(C\in C_Theta).                                \tag{QL.5}
```

Repeated entries on overlaps are identified.  Thus (QL.5) means equality of
all local restrictions, not a particular redundant serialization.

## 2. Exact local quotient

### Theorem QL.1 (connected-support Walsh quotient)

The isolated rooted states of the two pieces, restricted to the coordinate
subspaces occurring in `C_Theta`, together with (QL.5), determine exactly
the answer to every query in `Theta`.

More explicitly, for each `C` they reconstruct the rooted orbit triple of
the concatenated subtuple on `C`.  Equality of these local triples for two
systems is sufficient even if there is no single orthogonal map carrying
one **whole** marked tuple to the other.

Among carriers required to reconstruct every local rooted orbit triple,
(QL.5) is coarsest: with the isolated local states fixed, each displayed
restriction is recoverable from the corresponding concatenated local
triple.  This minimality is orbit-complete and query-local.  The scalar
Boolean maxima of a particular family may identify further states.

#### Proof

Fix `C` and restrict the coefficient maps in (QL.1) to its coordinate
subspaces.  The coefficient map of the local concatenated tuple is

```math
gamma_C(c,d)=alpha(c)+beta(d).                           \tag{QL.6}
```

The rooted amalgamation formulas give

```math
\begin{aligned}
G_C((c,d),(c',d'))
 &=G_a(c,c')+G_b(d,d')\\
 &\quad+kappa([c],[d'])+kappa([c'],[d]),\\
R_C&=gamma_C^{-1}(0),\\
Z_C&=gamma_C^{-1}(omega).                               \tag{QL.7}
\end{aligned}
```

Only the restrictions in (QL.5) occur in these formulas.  Conversely, the
off-diagonal block of `G_C` descends to `kappa|_(P_C times Q_C)`, while the
images of `R_C` and `Z_C` in `P_C direct-sum Q_C` are respectively the two
other restrictions in (QL.5).  This proves the reconstruction and converse.

The all-dimensional rooted Walsh orbit theorem now supplies, separately for
each connected support, an orthogonal coordinate relabelling carrying one
local landscape to the other.  A query energy is a sum over its disconnected
components, and its Boolean variables are disjoint between components.
Hence its maximum is the sum of the component maxima.  The orthogonal maps
may therefore be chosen independently on different components and even
differ from one declared query to another.  No compatibility datum between
supports that never occur in one connected component is observable.  The
same argument applies to the allowed componentwise rooted continuations.
`square`

For a purely unrooted component, the ambient-orbit theorem sharpens this
state exactly: its isolated state is only `(G,R)` and its relative local
datum is only

```math
(kappa|_(P_C times Q_C),J cap(P_C direct-sum Q_C)).     \tag{QL.7a}
```

The local root fibre is retained only on supports carrying a coordinate-
rooted or external-field continuation.  This is a semantic, not approximate,
collapse.  It follows by embedding every label as `(0,a)` in the full Walsh
coordinate space and applying the ambient Witt conjugacy separately on each
unrooted component.

### Complexity consequence

Put `p_C=dim(P_C)`, `q_C=dim(Q_C)`, and `d_C=p_C+q_C`.  A direct binary
serialization of the relative part of (QL.5) uses at most

```math
B_Theta
\le\sum_(C\in C_Theta)
\{p_Cq_C+d_C^2+d_C+1\}                                 \tag{QL.8}
```

bits: store the cross form, a row-reduced basis for the local coincidence
subspace, and either an empty flag or one representative of the local root
coset.  Overlap consistency can only reduce this count.

If every maximal connected support has size at most `w` and their total
incidence is `L=sum_C |C|`, then

```math
B_Theta=O(wL).                                          \tag{QL.9}
```

In particular, a partition, bounded-occurrence family, or any declared
support family with bounded `w` and `L=O(t)` has an exact `O(t)`-bit local
carrier for `t=r+s` labels, rather than the `O(t^2)` bits of the unrestricted
rooted orbit state.  Bounded `w` alone is not enough: the family of all pairs
has `L=Theta(t^2)`.  The theorem does not claim compression when one connected
future can contain all `t` labels.

The theorem also identifies why replacing `C_Theta` by its global linear
span is generally wasteful.  Such a replacement inserts cross-Gram,
coincidence, and root compatibility between supports that no declared query
ever joins.  The union of local orbit charts, not the span of their union,
is the query-generated object.

## 3. A sharp semantic lower bound for local coincidence memory

The preceding coarseness statement concerns local rooted orbits.  The next
result shows that its linear scaling for bounded connected supports is
already forced by ordinary scalar Boolean graph maxima.

### Theorem QL.2 (one exposed coincidence bit per path component)

For every `h>=1` there are two fixed isolated piece states and `2^h`
relative realizations with all of the following properties.

1. The isolated rooted states are identical for all realizations.
2. The complete cross form and complete combined root fibre are identical.
3. The only varying resource is the coincidence correspondence `J`.
4. There are `h` declared three-block path queries, each supported on a
   different component, whose response vectors are pairwise separated in
   `l_infinity` by at least

```math
Delta n^(3/2),
\qquad Delta={7-3sqrt3\over2}>0.                         \tag{QL.10}
```

Consequently any summary with one decoder answering every declared query to
uniform error `epsilon n^(3/2)`, where `epsilon<Delta/2`, requires at least
`2^h` states, or `h` bits.  The exact query-local carrier uses `O(h)` bits,
so its order of growth is sharp.  Here `n=q^2` is the order of one active
Walsh block, not the total number of variables across all marked ports.

#### Construction and proof

Choose an odd

```math
m=3h+r,
\qquad
r=\begin{cases}1,&h\text{ even},\\2,&h\text{ odd},\end{cases}              \tag{QL.11}
```

leaving `r` coordinates unused.  Put `q=2^m`, `n=q^2`, and in the `i`th
three-coordinate chunk define

```math
a_i=(1,1,1),
\qquad c_i=(0,0,1),                                    \tag{QL.12}
```

extended by zero outside that chunk.  The first piece consists of two
marked copies of every `a_i`.  For `sigma in \{0,1\}^h`, the second piece
has labels

```math
b_i^sigma=
\begin{cases}
a_i,&sigma_i=0,\\
c_i,&sigma_i=1.
\end{cases}                                             \tag{QL.13}
```

Every second-piece tuple is an independent orthogonal family of anisotropic
vectors: its Gram matrix is the identity and its relation kernel is zero.
The unused coordinates make every isolated and combined root fibre empty.
The first piece is fixed.  Moreover

```math
a_i\mathbin\cdot b_j^sigma=1_(i=j)                     \tag{QL.14}
```

for every `sigma`, so the full cross form is fixed.  On the other hand,
`b_i^sigma` belongs to the first-piece span exactly when `sigma_i=0`, and in
that case it equals the vector represented by either endpoint copy of
`a_i`.  Thus `J` records the `h` bits independently.

For query `i`, retain only the path whose endpoint labels are the two copies
of `a_i` and whose middle label is `b_i^sigma`; give every unused marked block
zero child and bridge weight.  If
`sigma_i=1`, the local word `(a_i,c_i,a_i)` has maximum
`7n^(3/2)/2`.  If `sigma_i=0`, the word `(a_i,a_i,a_i)` has maximum at most
`3sqrt3 n^(3/2)/2`.  These are the bent-witness saturation and
anticommuting spectral bound of Theorem 21.12, transported to the `i`th
coordinate chunk by a coordinate permutation.  Their difference is
(QL.10).

If two bit strings differ in coordinate `i`, query `i` separates their
answers by at least (QL.10).  Two landscapes assigned the same summary state
would receive the same decoded answer at every query, so one of the two
errors at coordinate `i` would be at least `Delta n^(3/2)/2`.  This proves
the packing lower bound. `square`

## 4. What the theorem does and does not say

The result gives an exact and strictly smaller composable state whenever the
declared future interaction hypergraph has small connected supports.  It
also proves that local coincidence memory cannot in general be discarded or
replaced by separately stored child states and cross-Gram data.

It does **not** solve the semantic minimality of the unrestricted Walsh graph
family.  In particular:

- arbitrary cross-form bits are known to be orbit-visible, but no matching
  scalar Boolean response packing is proved here;
- the characteristic-root fibre is required by a canonical rooted query but
  is exactly invisible to the entire declared unrooted weighted-graph class;
- if future connected supports grow to size `t`, (QL.8) returns the full
  quadratic orbit memory and supplies no compression.

Thus query-locality is a genuine strict theorem, not a route around the
dense connected-interface obstruction.  It identifies the exact boundary:
global `O(t^2)` compatibility is needed only to the extent that declared
future components can make the corresponding labels coexist.

The verifier checks local Gram/fibre reconstruction and every coefficient-
level bit of the packing through `h=6`.  The path-value separation imports
the independently verified bent witness and anticommuting ceiling from
Theorem 21.12; the finite script is a structural regression test rather than
a standalone exhaustive proof of all Boolean maxima.
