# Examples and obstruction atlas

This file records the smallest examples that forced the current definitions.
An equality of summaries is useful only after its resolution, normalization,
and permitted relabelings have been declared.

## Atlas

| Summary retained | What it really preserves | Counterexample or boundary | Missing information | Status |
|---|---|---|---|---|
| Support-sensitive energy entropy | location of every occupied energy interval | It fixes the limiting maximum for arbitrary landscapes | degeneracy and response geometry | Proved |
| Positive-rate upper-tail entropy | exponentially populated upper levels | It fixes every subsequential normalized maximum for homogeneous quadratics | labels, roots, and coupling response | Proved |
| Complete energy histogram | all one-replica energy multiplicities | order-eight signings below have the same histogram but one-vertex caps `16` and `20` | correlation of energy with the incident-sign feature | Exact computation |
| Complete energy--energy--global-overlap histogram | all two-replica queries depending only on the two energies and total overlap | left/right Curie landscapes have identical exact data but block response `1/2` versus `1` | decomposition of total overlap among labeled blocks | Proved |
| Complete code pair-distance enumerator | unrooted two-codeword geometry | `C^r` and `D^r` below have radii `2r` and `3r` | worst-root distance profile | Proved |
| Full spectrum | averaged Euclidean quadratic geometry | order-eight signing pair has caps `14` and `12` | orientation of cube vertices relative to eigenvectors | Exact computation |
| Every bounded restriction profile | all fixed local patterns and their induced moments | scalable Hadamard lift retains profiles through size six but has a leading cap gap | a growing zero-entropy resonance | Imported proved repository construction |
| Fixed traces, signed-graph densities, and finitely many replicas | bounded algebraic/traffic tests | bounded-operator Walsh planting has normalized caps at most `1/2` and at least `2/3` | a planted Boolean direction invisible at fixed arity | Imported proved repository construction |
| Strict high-temperature pressure on a bounded `beta` interval | smoothed bulk free energy | the same Walsh pair has asymptotically equal pressure there and a cap gap at least `1/6` | the zero-temperature endpoint | Imported proved consequence |
| One action profile without uniform integrability | weak profile laws | an archived scalable pair has vanishing action distance but separated caps | tail control for the energy observable | Imported proved repository construction |
| Sign-near posterior barycenter | almost every edge sign | it permits terminal rounding with negligible cap error | it forgets little: obtaining it costs `binom(n,2)-o(n^2)` bits | Imported/proved information converse |
| Full pinning response roof | every optimum after a linear spin field | it recovers `H(x)` at every cube vertex | nothing relevant; it is the full landscape | Proved |

The complete evidence cards, including normalizations and archive locations,
are in
[`drafts/obstruction_atlas_report.md`](drafts/obstruction_atlas_report.md).

## Example 1: pair-overlap data does not close a labeled coupling

Let `n=2m`, with even `m`, and split the coordinates into fixed labeled halves
`L` and `R`.  Define

```math
Q_L(x)=\frac{(\sum_{i\in L}x_i)^2-m}{m},
\qquad
Q_R(x)=\frac{(\sum_{i\in R}x_i)^2-m}{m},
```

and normalize each by its maximum `m-1`.  Swapping the two halves is a
bijection satisfying

```math
Q_R(\pi x)=Q_L(x),
\qquad
\langle\pi x,\pi y\rangle=\langle x,y\rangle.
```

Therefore the two landscapes have exactly the same finite energy histogram
and the same exact count of every triple

```math
(H(x),H(y),n^{-1}\langle x,y\rangle).
```

Keep the external apparatus fixed on `L` and impose zero left-block overlap.
Then

```math
\frac{1}{2(m-1)}
\max_{R_L(x,y)=0}\{Q_L(x)+Q_L(y)\}
=\frac{m-2}{2(m-1)}\longrightarrow\frac12,
```

whereas the same response for `Q_R` is exactly `1`.  The global overlap
summary forgot which species carries the energetic rigidity.

This pair is weighted and block-sparse, not a complete `+/-1` signing.  It is
an exact quadratic Boolean counterexample to universal sufficiency of global
pair data, not a claim about dense-sign realizability.

## Example 2: scalar entropy fails even for complete signings

In first-row-positive gauge, the two order-eight masks

```text
1466915
1068688
```

produce complete hollow sign matrices with the identical energy histogram

```text
-14:2, -12:6, -10:8, -8:12, -6:18, -4:22, -2:32,
0:40, 2:38, 4:34, 6:24, 8:12, 10:6, 12:2.
```

Both absolute caps are `14`.  Add one new vertex with the same all-negative
incident sign vector.  The extended caps are respectively `16` and `20`.
Thus the complete scalar density of states does not determine even the
simplest fixed-interface composition response.

This is an **exact finite computation**.  The matrices and exhaustive checks
are saved in
[`experiments/entropy_overlap_results.json`](experiments/entropy_overlap_results.json).

## Example 3: an exact tensor-stable coding collision

In the four-dimensional Hamming cube, set

```math
C=\{0000,0011,0101,0110\},
\qquad
D=\{0000,0011,0101,1001\}.
```

Both ordered distance enumerators are

```math
4+12z^2,
```

but their covering radii are `2` and `3`.  Cartesian powers therefore have
the same enumerator `(4+12z^2)^r` and radii `2r` and `3r`.

Let the energy of a cube point record whether it belongs to the code.  The
ambient cube is distance regular, so equality of the code--code enumerator
implies equality of every exact count

```math
#\{(x,y):(1_C(x),1_C(y),d(x,y))=(a,b,j)\}
```

for `C^r` and `D^r`.  Yet for `lambda>4r` the rooted query

```math
V_C(z;\lambda)=
\max_x\{\lambda 1_{C^r}(x)-d(x,z)\}
=\lambda-d(z,C^r)
```

has worst-root value `lambda-2r` for `C^r` and `lambda-3r` for `D^r`.
Unrooted pair geometry misses a linear worst-root separation.

This is a proof, not only a finite pattern.  It tensorizes at every `r`.

The phenomenon is not repaired by choosing a larger fixed replica count.  For
every `k`, Theorem 3.3 in [`theorems.md`](theorems.md) constructs parity-half-
cube codes whose complete ambient membership-and-distance histograms agree
through `k` points but whose covering radii differ.  Cartesian powers preserve
all those histograms and keep a positive normalized radius gap.  Thus an
unrooted finite overlap hierarchy is universally insufficient for code
covering radius; a rooted outer statistic or a structural synchronization
theorem is genuinely required.  Exact checks of the first members are in
[`experiments/code_replica_hierarchy_results.json`](experiments/code_replica_hierarchy_results.json).

## Example 4: why exact full pinning is not compression

Let `Omega={-1,+1}^n` and expose the feature `phi(x)=x`.  The feature hull is
the cube.  Every Boolean vertex is extreme, and a probability measure with
mean `x` must be the point mass at `x`.  Consequently the upper response roof
satisfies

```math
\widehat H_\phi(x)=H(x)
```

at all `2^n` vertices.  Any exact state sufficient for every linear pinning
query determines the entire landscape.  This is the clean information
boundary behind several bridge constructions that became equivalent to full
parent maximization.

## Example 5: state information versus landscape information

For every edge-sign vector `A`, shift the quadratic

```math
q_A(x)=a\sum_{i<j}A_{ij}x_ix_j
```

by its own maximum, obtaining `H_A=q_A-max q_A`.  All `2^binom(n,2)`
landscapes now have the same maximum, zero.  A strong field in direction `u`
pins the optimizer to `u`, and the degree-two Walsh coefficients of the
response recover all edge signs.  Uniform counterfactual response therefore
has quadratic information rate although the scalar maximum has zero rate.

By contrast, locating one arbitrary Boolean optimizer costs at most `n` bits.
The distinction is not “extrema are always expensive”; it is that a query
family rich enough to inspect every counterfactual direction can encode the
whole interaction skeleton.

## Exact finite census

After fixing first-row-positive gauge, the experiments exhaustively enumerated
every residual unlabeled graph representative through order eight using the
NetworkX graph atlas.  At order eight:

- `1044` rooted-gauge unlabeled graphs were checked;
- they produced `243` distinct exact
  `(H(x),H(y),x dot y)` histograms; and
- no two equal pair histograms had different **multisets** of one-vertex
  extension responses.

This is useful negative evidence about the smallest dense-sign collision; it
is not a theorem beyond order eight.  The full representative dataset is
[`experiments/quadratic_landscape_order8.json`](experiments/quadratic_landscape_order8.json).

## What distinguishes the examples

The examples point to three increasingly rich pieces of information:

1. **support resolution** distinguishes an absent extreme level from a
   subexponential cloud;
2. **rooted/interface geometry** records how energetic states correlate with
   the feature seen by a future coupling; and
3. **closure data** records correlations of features created by repeated
   composition.

The first is enough to recover a quadratic maximum from upper-tail entropy.
The second is needed by the explicit block and code queries.  Whether the
third can be compressed without rebuilding the whole landscape is the main
open structural question.
