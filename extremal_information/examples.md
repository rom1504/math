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
| Outer code distance spectrum | every root-averaged distance pressure | equal spectra can have appended-fragment radii `2` and `1` | alignment of distance layers with the syndrome group | Proved/exactly checked |
| Vanishing conditional species variance | average prediction from total overlap | a rare matching fibre keeps response error `rho/2` | uniform cross-root calibration on exposed fibres | Proved |
| Child response-separation polytopes | posterior width of each retained channel | same-space addition can cancel completely; max-plus gluing can make `Gamma={0}` | relative displacement/cross-Gram geometry | Proved |
| Convex signed-sum roof | every linear support query | identical zonotopes can have discrepancy gap `2d` | holes in a growing-dimensional reachable set | Proved |
| Exact tropical factor rank | separable channel count under uniform lattice-scale error | rank `r` can have rank-one normalized MSE `1/r` | query mass of the exposed anchors | Proved |
| Canonical code-transversal crossings | exponential exact channel count and raw sub-half-unit robustness | the graph-code table has rank-one normalized MSE `1/(16t)` | a positive-density joint witness, rather than a zero-density transversal | Proved |
| Syndrome block state | basis-versus-dense generator choices in `Theta(w)` direct-sum blocks | subset-selecting future fragments give `Omega(w)` response separation and a linear information rate | complexity beyond an explicitly supplied block decomposition | Proved |

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

## Example 6: equal outer spectra, different fibre composition

Over `G=F_2^3`, take parity-check column lists

```math
H_A=(1,2,3,4),
\qquad H_B=(1,2,4,7).
```

Their syndrome coset-leader profiles are

```math
(0,1,1,1,1,2,2,2),
\qquad(0,1,1,2,1,2,2,1).
```

Both have histogram `(1,4,3)` and their two-word kernels therefore have the
same complete outer polynomial

```math
2+8z+6z^2.
```

Append the same full-rank fragment `E=(1,3,5,6)`.  The first union of column
types omits syndrome seven and has covering radius two; the second contains
all seven nonzero types and has radius one.  The root-averaged pressure state
forgot the group label on each distance layer.  The labeled syndrome profile
repairs the composition without reconstructing the code or its full root
table.  The exhaustive check is in
[`phase2_code_syndrome_profiles_results.json`](experiments/phase2_code_syndrome_profiles_results.json).

## Example 7: average synchronization misses an exposed fibre

On `2m` states, choose a perfect matching.  Let `R_1` equal `rho` on one
matching edge and zero on all other off-diagonal pairs; let `R_2` equal
`rho` on the other matching edges and zero elsewhere.  Give both kernels
diagonal one.  They are PSD, every nonnegative mixture is exactly
ultrametric, and i.i.d. sampling gives weakly exchangeable Gram arrays.

For `q=(R_1+R_2)/2`, conditional variance on the matching fibre is

```math
{\rho^2\over m}\left(1-{1\over m}\right)\longrightarrow0.
```

Nevertheless that fibre contains both species-one values `rho` and zero, so
every scalar predictor `L(q)` has uniform error at least `rho/2`.  The two
signed zero-temperature queries expose the same gap.  The fibre is a set of
disconnected matching edges: local ultrametric no-crossing holds, while the
cross-root linkage in Theorem 9.2 fails exactly.

## Example 8: a convex roof can forget leading discrepancy

For every coordinate `j<=d`, compare vector pairs

```math
A_j=(2e_j,2e_j),
\qquad B_j=(3e_j,e_j).
```

Their concatenations have the identical zonotope `[-4,4]^d`, so every linear
support response agrees.  Yet

```math
\min_{s\in S(A)}\|s\|_1=0,
\qquad
\min_{s\in S(B)}\|s\|_1=2d.
```

The Shapley--Folkman response bound is therefore intrinsically
effective-dimension dependent.  In fixed dimension its error is
subextensive; when dimension grows with the number of summands, the holes can
carry a leading extremal gap.

## Example 9: uniform tropical robustness versus average erasure

Let `D_r` have zero diagonal and every off-diagonal entry one.  Its diagonal
crossing set proves min-plus rank `r`, and Theorem 11.1 shows that every
uniform approximation within error below `1/2` still needs at least `r`
terms.  The rank-one all-one matrix, however, differs only on the diagonal:

```math
{1\over r^2}\|D_r-\mathbf1\|_F^2={1\over r}\to0.
```

Thus neither exact rank nor a uniform fooling-set gap implies an average-loss
lower bound.  The query law must give quantitative mass to witnesses that
remain monochromatic under every small-channel factor assignment.

## Example 10: exponential code rank with vanishing diffuse error

Let

```math
C_t=\{(z,z):z\in\mathbb F_2^t\}
```

and split its `2t` coordinates into the two displayed blocks.  Its normalized
conditional distance table is

```math
M_t(x,y)={d_H(x,y)\over2t}.
```

The exact min-plus rank is `2^t`, and the canonical transversal crossing set
remains rank-obstructing under raw uniform error below `1/(4t)` after this
normalization.  Under independent uniform `x,y`, however,
`d_H(x,y)` is `Binomial(t,1/2)`.  The rank-one constant table `1/4` therefore
satisfies

```math
\mathbb E\left(M_t-{1\over4}\right)^2={1\over16t}\longrightarrow0.
```

The query-weighted four-cell certificate sees the same loss of mass: even for
one channel and the complete transversal witness graph its value is exactly
`2^{-t}/8`.  Thus exact algebraic size, uniform lattice-scale resolution, and
diffuse normalized rate--distortion are genuinely different complexity
coordinates.

## Example 11: joint code queries amplify microscopic support bits

Decompose `F_2^w` into `q=w/L` fixed `L`-dimensional blocks.  In each block a
latent bit chooses either a basis, whose Cayley covering radius is `L`, or all
nonzero vectors, whose radius is one.  A legal appended fragment indexed by
`P subset [q]` leaves the basis choice visible exactly on the selected blocks.
The resulting radius is

```math
q+(L-1)|\{j\in P:a_j=0\}|.
```

Thus a single query can coherently add any selected subset of the latent
effects.  Hamming packings and the decoder triangle inequality imply that,
for every fixed `epsilon<1/8`, uniformly answering all future-fragment radii
within `epsilon*w` needs `Omega_epsilon(w)` bits on this family.  The latent
block vector gives a matching `O(w)` exact state.  This is the positive
counterpart to Example 10: query mass erased its pairwise transversal, while
an algebraically legal joint query here sums the effects before distortion is
paid.
