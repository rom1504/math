# Switching gauge quotients for optimized bridges

**Status.** Director-critical scope theorem.  The statements below are
elementary but exact.  They distinguish a labelled fixed-bridge response
from the unlabelled bridge optimization that occurs in the original signing
problem.  No claim is made that the surviving bridge fibre is response
separated.

## 1. Exact covariance of a block parent

For a hollow complete signing `A` of order `m`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_{x\in\{+-1\}^m}|H_A(x)|.
```

For complete signings `A,C` of orders `m,n` and an exact-sign bridge
`B in {+-1}^{m by n}`, write

```math
P(A,C;B)=
\begin{pmatrix}A&B\\B^T&C\end{pmatrix}.            \tag{GQ.1}
```

For sign vectors `s,t`, let `S=diag(s)`, `T=diag(t)`, and set
`A^s=SAS`, `C^t=TCT`.

### Theorem GQ.1 (independent child switches are pure gauge after bridge transport)

For every `A,C,B,s,t`,

```math
P(A^s,C^t;SBT)
=
\begin{pmatrix}S&0\\0&T\end{pmatrix}
P(A,C;B)
\begin{pmatrix}S&0\\0&T\end{pmatrix}.             \tag{GQ.2}
```

Consequently the two parents have exactly the same Boolean energy multiset,
not only the same cap:

```math
H_{P(A^s,C^t;SBT)}(x,y)
=H_{P(A,C;B)}(s\mathbin\odot x,t\mathbin\odot y),
\qquad
Q(P(A^s,C^t;SBT))=Q(P(A,C;B)).                     \tag{GQ.3}
```

#### Proof

Block multiplication gives (GQ.2).  The displayed change of Boolean
variables is a bijection, proving (GQ.3). `square`

The bridge must be transported.  Holding `B` fixed is not a symmetry in
general; Section 4 gives the smallest useful falsifier.

## 2. Minimal covariance hypothesis for optimized bridge languages

Let `\mathcal B(A,C)` be any allowed bridge family.  It may depend on the
two children and may contain exact signs, weighted bridges, sparse bridges,
or bridges satisfying additional constraints.  Call it **switching
covariant** if

```math
\boxed{
\mathcal B(A^s,C^t)=S\mathcal B(A,C)T
 :=\{SBT:B\in\mathcal B(A,C)\}}
                                                        \tag{GQ.4}
```

for every `s,t`.  Equality, rather than one inclusion, is the useful minimal
hypothesis: it says that bridge transport is a bijection between the two
feasible fibres.

### Theorem GQ.2 (optimized composition factors through child switching classes)

If (GQ.4) holds, then the attainable parent-cap multisets agree:

```math
\{Q(P(A^s,C^t;B')):B'\in\mathcal B(A^s,C^t)\}
=
\{Q(P(A,C;B)):B\in\mathcal B(A,C)\}.               \tag{GQ.5}
```

Therefore every relabelling-invariant functional of this multiset is
switching invariant.  In particular,

```math
F_{\mathcal B}(A,C)
:=\inf_{B\in\mathcal B(A,C)}Q(P(A,C;B))             \tag{GQ.6}
```

satisfies

```math
F_{\mathcal B}(A^s,C^t)=F_{\mathcal B}(A,C).        \tag{GQ.7}
```

The same conclusion holds for a minimum, supremum, maximum, quantile, or
expectation under a bridge law transported by `B mapsto SBT`.  It also holds
with an auxiliary bridge cost, provided that cost is invariant under the
same transport.

#### Proof

By (GQ.4), `B mapsto SBT` is a bijection of feasible fibres.  Equation
(GQ.3) preserves each value under that bijection, which proves (GQ.5) and
all its consequences. `square`

The unrestricted exact-sign family is switching covariant.  So are
operator-norm, rank, fixed-support, and other constraints that themselves
survive left--right diagonal signs.  Thus the response in the question,

```math
F(A,C)=\min_{B\in\{+-1\}^{m\times n}}Q(P(A,C;B)),   \tag{GQ.8}
```

depends only on the two coordinate-labelled switching classes.  Even if
the future block `C` is kept as a fixed public query, optimizing `B`
separately for each query makes every switching of `A` have the same
response to `C`.

Covariance is a structural guarantee, not a logically necessary condition
for one accidental equality of minima.  Without it there is no universal
invariance theorem, as Section 4 shows.

## 3. Exact classification of labelled and unlabelled future languages

The distinction is most transparent for an abstract contextual language.
Let a group `G` act on systems `X`.  Let `K_X` be the admissible context
fibre and `V(X,k)` its scalar value.  Suppose there are bijections

```math
\tau_g:K_X\longrightarrow K_{gX},
\qquad
V(gX,\tau_gk)=V(X,k),                              \tag{GQ.9}
```

with the usual composition law.  Put `r_X(k)=V(X,k)`.

### Proposition GQ.3 (what quotients a gauge orbit)

Under (GQ.9):

1. The labelled response table is transported, not fixed:

   ```math
   r_{gX}(k')=r_X(\tau_g^{-1}k').                  \tag{GQ.10}
   ```

2. A scalar response `\Lambda(r_X)` factors through `X/G` whenever
   `\Lambda` is invariant under this relabelling.  Optimization over the
   whole fibre, the unordered response multiset, and integration against a
   transported measure are examples.

3. A fixed-label public language factors through `X/G` **if and only if**
   its table is pointwise invariant:

   ```math
   V(gX,k)=V(X,k)
   \quad\hbox{for every declared fixed label }k.    \tag{GQ.11}
   ```

   Equivalently, in a common invariant context fibre,
   `r_X=r_X circ \tau_g^{-1}`.  A single discrepancy

   ```math
   |V(gX,k)-V(X,k)|\ge\delta                        \tag{GQ.12}
   ```

   is an exact contextual falsifier at every error below `delta/2`.

#### Proof

Equation (GQ.10) is (GQ.9) with `k=\tau_g^{-1}k'`.  Applying a
relabelling-invariant `\Lambda` proves item 2.  Item 3 is exactly the
definition of equality of fixed-label response tables, and the triangle
inequality gives the distortion claim. `square`

For block signings, the natural context is `(C,B)` and

```math
\tau_s(C,B)=(C,SB)                                 \tag{GQ.13}
```

when only the first child is switched.  Hence the following languages have
different answers.

| Future language | Does it quotient the child switch? | Reason |
|---|---:|---|
| all bridges, minimized separately for each parent/query | yes | (GQ.4)--(GQ.7) |
| any left--right switching-covariant bridge family, minimized | yes | transported feasible-fibre bijection |
| invariant bridge distribution and any distributional statistic | yes | measure-preserving transport |
| the complete context orbit, but context labels forgotten | yes | response multiset is merely permuted |
| a fixed bridge `B_0`, with labels held fixed | not in general | the required context is `SB_0`, not `B_0` |
| a non-invariant public bridge bank | not in general | (GQ.11) is an additional substantive condition |
| one bridge shared across several independently switched alternatives | only for a common switch | one `B` cannot absorb several different left gauges |
| a gauge-covariant bank with labels co-switched | yes | comparisons are made after `\tau_s` |
| the same bank with labels pinned in the laboratory frame | only if (GQ.11) happens | relative gauge is observable |

Thus “public” by itself is not the dividing line.  The dividing line is
whether the future is transported/optimized and then quotiented, or instead
provides a fixed reference frame against which a relative switching label
is measured.

## 4. A four-vertex fixed-bridge falsifier

Let `A=J_3-I_3`, let `C` be the hollow order-one matrix, and take the fixed
bridge

```math
B_0=(-1,-1,-1)^T.
```

Then `P(A,C;B_0)` is a switching of the all-positive order-four signing, so

```math
Q(P(A,C;B_0))=6.                                   \tag{GQ.14}
```

Now switch only `A` by `s=(-1,-1,1)` but hold `B_0` fixed.  The six edge
signs are

```text
12:+, 13:-, 23:-, 14:-, 24:-, 34:-.
```

Fixing the fourth spin to `+1` (global spin inversion loses nothing), the
eight energies are

```text
2, 4, 0, -2, 0, -2, 2, -4,
```

and hence

```math
Q(P(A^s,C;B_0))=4.                                 \tag{GQ.15}
```

This is a fixed-label response gap of two.  There is no contradiction with
GQ.1: transporting the bridge to `SB_0=(1,1,-1)^T` restores cap six.
Thus noncovariance can expose switching gauge already at order four.

## 5. The residual bridge fibre: two different quotients

Quotienting child switches does **not** quotient away an arbitrary dense
bridge.  Two natural counts must be kept separate.

### Proposition GQ.4 (a bare bipartite bridge modulo row/column gauges)

The action

```math
B\longmapsto SBT,
\qquad (s,t)\in\{+-1\}^m\times\{+-1\}^n          \tag{GQ.16}
```

on a bare exact-sign `m by n` bridge has effective orbit size
`2^(m+n-1)`.  Therefore it has

```math
2^{mn-(m+n-1)}=2^{(m-1)(n-1)}                     \tag{GQ.17}
```

orbits.  A complete invariant is given by the rectangle products

```math
B_{ij}B_{i1}B_{1j}B_{11},
\qquad 2\le i\le m,\quad2\le j\le n.              \tag{GQ.18}
```

#### Proof

If `SBT=B`, then `s_it_j=1` for every `i,j`, so `s,t` are the same constant
sign.  This is the two-element kernel, proving the orbit count.  Row and
column switches can make the first row and column positive, after which the
remaining entries are exactly (GQ.18). `square`

### Proposition GQ.5 (the bridge between anchored complete child classes)

Fix coordinate-labelled switching classes `[A]` and `[C]` of complete
signings of orders `m,n>=2`.  Consider all triples `(A',C',B)` with
`A' in [A]`, `C' in [C]`, and arbitrary exact-sign `B`, modulo switching of
the complete parent.  After choosing fixed representatives `A,C`, every
class has a representative `(A,C,B)`, and the only residual identification
is

```math
B\sim-B.                                           \tag{GQ.19}
```

Consequently the anchored bridge fibre has exactly

```math
\boxed{2^{mn-1}}                                   \tag{GQ.20}
```

classes.

#### Proof

Diagonal switching acts freely and transitively on the projective switching
orbit of a complete signing: if `SAS=A`, completeness and nonzero entries
give `s_is_j=1` for every `i ne j`, hence `s` is constant.  Thus there are
unique projective switches returning `A'` and `C'` to the chosen
representatives.  Choosing the opposite sign-vector representative on
exactly one shore leaves its internal block fixed and sends `B` to `-B`;
choosing it on both shores fixes `B`.  There are no other identifications.
This proves (GQ.19)--(GQ.20). `square`

The smaller count (GQ.17) is appropriate only when row and column gauges of
the **bare bridge** remain disposable.  Once complete internal blocks anchor
those vertex gauges, the `m+n-2` row/column alignment bits are physical in
addition to the `(m-1)(n-1)` rectangle holonomies; only one global bridge
sign remains redundant.  Conversely, (GQ.20) is a coefficient-fibre count,
not a proof that all `2^(mn-1)` choices have separated optimized caps.
The exponent identity is

```math
mn-1=(m-1)(n-1)+(m-1)+(n-1).                       \tag{GQ.21}
```

## 6. Consequence for Theorem 36.28 and the original optimization

Theorem 36.28 is mathematically correct for its declared language.  It holds
one common bridge `B` fixed while the child and query fill range over matched
switching poles.  Switching only the child naturally transports `B` to
`SB`; transporting both a child pole `s` and a fill pole `q` transports it
to `S_sBS_q`.  The theorem keeps `B` fixed instead.  These transported
bridges are generally absent from its bank.  The theorem therefore measures
a **relative gauge against a public reference bridge**.  Its `Theta(n)`
response rate is real for fixed public bridges, shared-bridge compilers,
pinned laboratory-frame queries, and any architecture that must answer all
those labelled futures.

It does **not**, by itself, obstruct the original unlabelled cross-order
optimization.  In that optimization the bridge is a decision variable for
each candidate parent.  The map `B mapsto SBT` is a bijection on all exact
sign bridges, so GQ.2 collapses the entire switching orbit of each child to
one state.  In particular, Theorem 36.28 does not rule out a composition
state that stores no child switching label and chooses the bridge after a
gauge representative has been selected.

This scope correction does not make optimized composition easy.  After the
two child gauges are fixed, the arbitrary bridge still has the
`2^(mn-1)`-element anchored coefficient fibre in GQ.5, and minimizing over
that fibre can retain essentially the full parent Boolean difficulty.  The
correct conclusion is therefore:

```text
Theorem 36.28 proves labelled fixed-reference response heaviness;
it does not prove unlabelled optimized-bridge incompressibility.

Switch labels are pure gauge for the original bridge minimum;
the optimized bridge/holonomy fibre is the surviving obstruction.
```

The theorem remains relevant to a proposed recurrence only after that
recurrence is shown to require a common/noncovariant bridge across several
alternatives, or another fixed reference frame.  A recurrence allowing an
independently optimized covariant bridge cannot inherit its lower bound.

## 7. Archive comparison

- The fixed-bridge observation already appears qualitatively in the
  quadratic dense-bridge compression audit: switching changes a generic
  fixed bridge from `B` to `SB`.
- The polynomial fractional-bridge report proves exact covariance for one
  special universal-double family and explains why randomizing a covariant
  signed permutation only renames states.
- Theorem 21.41 develops projective signed-graph gauge and cycle holonomy for
  a structured common-Hadamard block family.
- Theorem 36.28 proves the new all-bounded-cap fixed-reference response bank,
  and its audit explicitly exempts narrower co-switching languages.

GQ.1--GQ.3 make the general optimized-versus-labelled distinction explicit;
GQ.4--GQ.5 audit the two easily confused bridge-fibre counts; and Section 6
records the resulting scope correction.  These are not a new lower bound or
a recurrence theorem.
