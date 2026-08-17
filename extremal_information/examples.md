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
| Raw mixed-holonomy dimension | `D kappa` exact kernel-valued gluing coordinates | discrete and two-scale carriers collapse all or almost all endpoint responses | Hausdorff separation in the query metric and cheap carrier presentation | Proved |
| Carrier Hausdorff geometry | all distance-to-carrier endpoint responses | presentation cost can erase an arbitrary carrier; diffuse laws can give witness balls negligible mass | access cost and query exposure | Proved |
| Projected carrier | distance geometry after a coarse metric quotient | identity quotients and high-entropy projected hyperspaces give no compression | small fibres are insufficient without entropy reduction and descended composition | Proved |

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

## Example 12: locally trivial fragments create full-rate information

Let `W=F_q^D` and `Q=F_q^k`.  Pair a zero lift and a lift by `v_j in W` over
each quotient basis direction.  Either fragment alone is removed by a linear
kernel-fixing shear.  After union, however, the kernel-endpoint profile is

```math
F_V(u)=\min_{z\in F_q^k}
\{2\operatorname{wt}(z)+\operatorname{wt}(u+Vz)\}.
```

Inside an asymptotically good `[D,D/4,>D/8]` host code, its `k`-subspaces
give at least `q^(3Dk/16)` profiles separated by more than `D/16` whenever
`k<=D/32`.  Thus unlabelled endpoint queries retain `Theta(Dk log q)` bits:
a constant fraction of the exact compatibility information born at
composition.

The exact profile is more detailed than its span.  In the binary case it
recovers precisely the set of distinct generator columns of Hamming weight at
least three; general changes of channel basis are therefore not exact
symmetries.  At macroscopic scale it lies within `2 rank(V)` of distance to
`im V`, which is why Grassmannian packing becomes the right coarse geometry.

## Example 13: the same gauge dimension can collapse completely

Two counterexamples delimit Example 12.

1. Give `F_q^D` the discrete metric in which every nonzero kernel element has
   cost one.  Every mixed profile is simply `F_V(u)=1_(u ne 0)`, independent
   of all `Dk` holonomy coordinates.
2. Let `varpi:F_q^D->F_q^r` and set

   ```math
   d(x,y)=D\mathbf1_{\varpi x\ne\varpi y}
          +\mathbf1_{x\ne y}.
   ```

   This metric still has linear diameter.  For `k=o(D)`, every endpoint
   profile is within `2k+1` of a decoder indexed only by the subspace
   `varpi(im V)`.  With fixed `q,r`, that is a constant-size macroscopic
   quotient.  Retaining the labeled map `varpi V` gives all-future error one
   with only `q^(rk)` states, versus `q^(Dk)` exact gauges.

The missing information is not relation rank but fine geometry inside
diameter-one metric fibres.  This is the sharp negative evidence that forced
the presented-carrier and metric-quotient laws.

## Example 14: rank geometry is full-rate, Hamming has a duality gap

In `End_(F_q)(F_(q^D))` with rank distance, the Gabidulin space

```math
\left\{x\mapsto\sum_{i<r}a_ix^{q^i}:a_i\in F_{q^D}\right\}
```

has dimension `rD` and minimum nonzero rank `D-r+1`.  Its `k`-subspaces
therefore produce `q^(k(rD-k))` mixed-channel profiles separated by at least
`D-r+1-2k`.  With `r=floor(D/2)` and `k<=D/16`, this is an
`Omega(D^2k log q)` macroscopic response lower bound.  Row-supported
anticodes have dimension `Da`; the code--anticode inequality proves this is
the exact anticode dimension at rank scale `a`.

Binary Hamming geometry behaves differently.  At scale
`a=floor(delta D)`, its largest linear anticode has dimension `a`, but the
Hamming sphere-packing bound gives

```math
D-A_W(a)-s_W(a)
\ge\bigl(H_2(\delta/2)-\delta-o(1)\bigr)D.
```

Thus the lower common-host certificate and the optimal quotient certificate
have a genuine leading gap.  This does not yet prove an equally large gap in
the complete response entropy: subspace carriers can be packed without all
lying inside one common separated code.  It is the first rigorous example in
the framework where a third geometric invariant is demonstrably needed.

## Example 15: identical quotient geometry, different rooted response

In one four-bit Hamming block compare

```math
C^{(2)}=\operatorname{span}(1100),qquad
C^{(1)}=\operatorname{span}(1000).
```

Both quotient normed spaces are exactly the three-bit Hamming cube.  Hence
their complete quotient leader distributions and all sparse-flat spectra
agree.  After `r` direct sums, however, the maximum weight inside the first
carrier is `2r` and inside the second is `r`.  Their Hausdorff distances to
the zero carrier differ by `r`.

This is a scalable rooted-versus-unrooted counterexample.  Knowing every
short-syndrome flat controls directed distance *to* the stored carrier, but
does not control reverse or future-carrier queries.  The internal kernel
metric, and in general its coupling to quotient lifts, is genuine response
information.

## Example 16: one-channel carrier entropy is nonlinear coding entropy

Binary one-dimensional carriers are indexed by nonzero words `v`.  Their
Hausdorff metric is

```math
d_H(\operatorname{span}(v),\operatorname{span}(w))
=\max\{\min(wt(v),wt(v+w)),\min(wt(w),wt(v+w))\}.
```

The maximum size of a packing at distance greater than `t` lies between
`A_2(D,t+1)-1` and `A_2(D,t+1)`.  Thus its asymptotic packing exponent is the
unrestricted binary coding rate, not merely the best linear-code rate.
Puncturing gives `2^(D-t)` quotient states.  A probabilistic cover built from
the exact line-ball volume has exponent at most `1-H_2(delta)`, proving the
genuine state overcount

```math
(H_2(\delta)-\delta-o(1))D.
```

This holds at `t=delta D` and proves that the code--anticode gap is present
in the actual response entropy already at one channel.  For `k=Theta(D)`, a
systematic chart instead produces a growing alphabet of size `2^k`; ordinary
alphabet coding meets Singleton, so only coherent same-input recoupling can
improve the puncturing exponent there.

## Example 17: a finite directed alphabet amplifies before absolute values

For local presented responses `f_a`, retain only

```math
r(a,b)=\sup_x(f_a(x)-f_b(x)).
```

Under direct-product composition, this finite table closes exactly:

```math
\|F_{\boldsymbol a}-F_{\boldsymbol b}\|_\infty
=\max\left\{
 \sum_i r(a_i,b_i),
 \sum_i r(b_i,a_i)
 \right\}.
```

The seven nonzero `[7,3,4]` simplex words give Hamming line carriers with
two-sided carrier gap four and presentation radius two.  A seven-letter
outer code of relative distance `3/4` therefore gives exponentially many
responses separated by `3m/2`.  The seven nonzero multiplication maps of
`F_8` give rank-metric lines with gap three and the same presentation cost,
yielding separation `3m/4`.

This is a same-sign joint inequality: matching blocks pay nothing, and all
differing blocks are evaluated in one orientation before the final absolute
value.  It is stronger by a linear amount than applying the global carrier
bound and paying every block separately.

## Example 18: pure Grassmannian packing can escape every common host

Seven explicit two-planes in `F_2^6` have both directed Hausdorff distances
equal to three.  Concatenating them through a seven-letter outer code gives
`2^((0.0573549...-o(1))m)` carriers in `Gr_(2m)(F_2^(6m))` at relative
Hausdorff distance tending to `3/8`.  Every product carrier contains a
nonzero vector of weight at most four, so none lies in a common host whose
minimum distance grows with `m`.

This falsifies common-host control of pure carrier packing.  It does not by
itself give presented responses: its standard presentation radius is `4m`,
larger than its carrier gap.  The simplex lines in Example 17 are the
presentation-compatible replacement.

At `(D,k,t)=(5,2,2)`, four explicit two-planes already have pairwise
Hausdorff distance three, whereas a distance-three linear host has dimension
at most two and contains only one two-plane.  Conversely, the general
injection-distance construction based only on counting low words in
`C+C'` is provably bounded by the common-host Gilbert exponent.  Placement
of those low words relative to both rooted subcodes, not their scalar count,
is the missing information.

## Example 19: pure Max-Cut realizes every projective separator table

A width-`w` Max-Cut boundary has `2^(w-1)` spin assignments after quotienting
global flip.  Positive-edge pinning gadgets expose every one of these classes,
so the conditional cut profile is the exact contextual state.

More strongly, every nonnegative function `F` on the projective cube is the
response shape of a positive-weight Max-Cut graph of treewidth at most
`w+1`.  A common anchor converts boundary spins into gauge-relative spins;
one selector per oriented word implements a lookup table; and direct edges or
two-edge paths convert every resulting signed Ising term to positive cut
edges.  After a common padding constant, the whole cube `[0,W]^(2^(w-1))` is
realizable.  Thus its response description cost is

```math
Theta(2^(w-1) log(W/epsilon))
```

bits at additive error `epsilon`.  The construction has exponential size and
boundary load, so this is a worst-case contextual lower bound rather than a
unit-sensitivity approximation lower bound.

## Example 20: interface regularity changes approximate response complexity

For arbitrary width-`w` tables, additive error `epsilon w` still costs
`Theta(2^w)` bits.  If the response shape is one-Lipschitz in boundary
Hamming distance, values on a Hamming net suffice.  Midpoint McShane
extensions give

```math
log_2 Cov_(epsilon w)
 <=2^((1-H_2(epsilon)+o(1))w),
```

while balanced two-level values on an `epsilon w`-separated code give

```math
log_2 Pack_(epsilon w)
 >=2^((1-H_2(epsilon)+o(1))w).
```

Using a code separated by more than `2 epsilon w` gives the operational
covering lower bound

```math
log_2 Cov_(epsilon w)
 >=2^((1-H_2(2 epsilon)+o(1))w).
```

The same contextual quotient can therefore have different lossy state-growth
laws under different regularity promises.  Exact coordinate exposure alone
does not determine approximate complexity.

## Example 21: metric holonomy and a tropical bottleneck

For a finite metric space `(Y,d)`, isometry `g`, and strength `lambda`, the
kernel

```math
D_(lambda,g)(a,t)=lambda d(t,g(a))
```

has the exact interacting law

```math
D_(lambda,g) star D_(mu,h)
=D_(min(lambda,mu),h circ g).
```

An arbitrary chain is represented by its weakest link and accumulated
isometry holonomy, and its directed row gaps are
`min_i(lambda_i)d(a,b)`.  Uniform entrywise perturbations accumulate only
additively.  This is a nonproduct validation of the directed-response
framework on every finite metric, not only on a binary chain.

For a `q`-state edge preferring `t=pi(a)` with strength `J`, min-plus serial
composition has the exact law

```math
K_(J,pi) star K_(L,rho)
=-max(J,L) 1+K_(min(J,L),rho circ pi).
```

An arbitrary chain is represented by a baseline, the permutation product,
and the weakest link.  All off-diagonal directed row responses equal that
bottleneck.  This is a nonproduct validation of the directed-response
framework: continuation preserves a directed face exactly while every
relevant input remains exposed, and clips it when a weak edge hides those
witnesses.  The binary case is the familiar signed zero-temperature Ising
chain, but the permutation law is nonabelian for `q>=3`.

## Example 22: a unit-load distance shell hides a universal compiler

Every projective one-Lipschitz table on a width-`w` separator is the response
shape of a pure weighted Max-Cut component with boundary load one at every
vertex.  Compile the arbitrary table at a private interface, then connect
each exposed spin to its private copy by a two-edge unit path.  Maximizing the
middle spins applies the max-plus distance projector

```math
(P_df)(x)=max_y(f(y)-d(x,y))=f(x).
```

Conversely, boundary load one forces precisely this Lipschitz condition.
Thus normalized boundary sensitivity does not reduce the realizable response
class below the full Lipschitz ball.  At error `epsilon w` it still needs
exponentially many response bits:

```math
2^((1-H_2(2epsilon)+o(1))w)
 <=log_2 Cov_(epsilon w)
 <=2^((1-H_2(epsilon)+o(1))w).
```

The compiler is exponentially large, so polynomial component size is a
genuinely stronger resource promise.  The shared-parameter theorem now
quantifies it: an `m`-edge component has only
`O_epsilon(m^2+m log(w+m))` coarse response bits.

## Example 23: a vanishing idempotence defect accumulates a transition toll

On the line `0,...,q-1`, put

```math
K_\delta(i,j)=a|i-j|-\delta\quad(i\ne j),
\qquad K_\delta(i,i)=0,
```

with `a>delta`.  Its min-plus idempotence defect is exactly `delta`, and it
is within `delta` of the exact path metric.  Nevertheless,

```math
K_\delta^{star T}(i,j)
=a|i-j|-\delta\min(T,|i-j|)
```

off the diagonal.  The row-shape drift is `(T-1)delta/2` up to depth
`q-1`.  Taking `a=2c/q`, `delta=c/q`, and `T=q-1` gives bounded-diameter
kernels with vanishing local defect but fixed response drift.

The lost datum is a reward of `delta` for each useful nonzero transition.
Composition turns that microscopic toll into macroscopic information.  Thus
one-step closeness to a metric shell does not imply a depth-stable quotient;
bounded useful path length or an exact metric semilattice is a necessary
kind of additional hypothesis.

## Example 24: polynomial presentation cannot hide a universal separator

Fix a pure Max-Cut topology with `m` edges and arbitrary nonnegative real
weights.  Every boundary response is a maximum of binary-incidence linear
forms in the same `m` parameters.  Although there may be exponentially many
private cuts and boundary queries, their optimizer comparisons have only
`(3^m-1)/2` possible hyperplane normals.  Arrangement faces and
finite-dimensional volume therefore give, under unit boundary load,

```math
\log_2 Cov_(epsilon w)
=O_epsilon(m^2+m log(w+m)).
```

This does not assume bounded numerical precision.  In contrast, the full
unit-load Lipschitz response class has exponentially many response bits.
Any family approximating every such response must have

```math
liminf (log_2 m)/w
 >=(1-H_2(2epsilon))/2.
```

The distance shell can hide a universal compiler behind a cheap interface,
but the compiler itself cannot have polynomial max-affine presentation
complexity.

## Example 25: four weighted states lump to two tropical aggregates

Partition four max-plus automaton states into `I_0={1,2}` and `I_1={3,4}`.
Choose transition rows that are microscopically different inside each block
but whose maximum into every target block equals the same quotient
transition `S(a,b)`.  Then no raw state can simply be deleted: depending on
the incoming vector, either member of a block may win.  Nevertheless every
future suffix depends only on

```math
(P_0,P_1)
=(max(p_1,p_2),max(p_3,p_4)),
```

and every appended letter updates this pair by max-plus multiplication with
the `2 by 2` quotient matrix.  Suitable quotient suffixes expose both
coordinates, giving exact response complexity
`Theta(2 log(1+B/epsilon))` on a bounded box.

This is a strict composable state that is not a boundary table.  It is
predicted by contextual response plus derivative compatibility: the
residual quotient says what futures can see, while tropical lumpability says
why the same quotient remains closed after every continuation.

## Example 26: gauges telescope, while weak links reset memory

Let `D_(lambda,g)(a,b)=lambda d(b,g(a))` be a metric shell.  Perturb every
layer by compatible endpoint potentials,

```math
K_t(a,b)=D_(lambda_t,g_t)(a,b)
          +phi_(t-1)(a)-phi_t(b)+c_t.
```

All internal potentials cancel before the min-plus minimum, so an arbitrary
chain is exactly

```math
D_(min_t lambda_t,g_T circ ... circ g_1)
+phi_0(a)-phi_T(b)+sum_t c_t.
```

The directed row-table error is only
`phi_0(a)-phi_0(a')`, independent of chain length.  Vanishing rectangular
and adjacent-interface circulations are the exact finite certificate for
this cancellation.  A nonzero circulation on a pumpable cycle instead grows
linearly when the cycle is repeated.

There is a second mechanism.  A max-plus transition whose projective image
has diameter `rho` forgets its input up to `rho`.  If such a transition block
recurs every `L` steps, a local quotient defect `epsilon` stays below
`rho+2L epsilon` after the first reset.  For a binary zero-temperature Ising
bond of strength `J`, the image diameter is `2|J|`: a weak bond is a small
reset even though every nonzero bond has global projective Lipschitz
coefficient one.

## Example 27: query packing is exactly Lipschitz exposure

For a finite metric query space `X`, choose a set `C` whose points are at
least `2gamma` apart.  On every balanced split of an even subset of `C`, put
values `+gamma` and `-gamma`; McShane extension realizes all those patterns
as one-Lipschitz responses.  Conversely, if a response language realizes all
balanced margin patterns on two queries, Lipschitzness forces those queries
to be at least `2gamma` apart.

Hence the balanced exposure dimension of the full Lipschitz ball is exactly

```math
2 floor(Pack_X(2gamma)/2).
```

This recovers the normalized Max-Cut lower entropy from a Hamming packing of
the separator itself.  In a weighted automaton, robust suffix pins instead
expose a coordinate cube; upper-orthant VC dimension shows that no more than
the raw number of coordinates can be independently exposed.  The result is
static: neither exposure certificate supplies a derivative congruence.

## Example 28: facets create robust shared-parameter responses

Let `P=conv(V)` be a full-dimensional `0/1` polytope.  Use `V` as one
max-affine query and, for every facet `F`, use the witness set obtained by
deleting all vertices of `F`.  At the normalized outward facet normal, every
query except the matching deletion has one common value; the matching query
is exactly one lower.  Thus the response shape is

```math
c_F 1-e_F.
```

Different facets have shape distance one and radius one half.  Known
`0/1` polytopes with `(cm/(log m)^2)^(m/2)` facets therefore force
`Omega(m log m)` response bits at every fixed error below one half, even
with arbitrary real parameter precision.  This validates shared-parameter
complexity as a semantic resource, while leaving a genuine gap to the
`O(m^2)` common-arrangement upper bound.

## Example 29: mean-field curvature synchronizes a histogram to one number

For heterogeneous binary local fields, the score at fixed occupancy `k` is
the sum of the `k` largest fields.  Hence one scalar chemical potential
recovers the full sorted-field profile, and anonymous block union becomes
sorted multiset union.  Rounding every microscopic field once to a common
grid gives a histogram homomorphism: its `O((1+B/eta)log n)` bits incur at
most `eta n/2` response error at mass `n`, regardless of merge depth.

Now add the same pair reward `J {k choose 2}` to every block and parent.
Linear terminal fields see only the concave envelope of the conditional
profile, but bilinearity makes that envelope an exact congruence under
composition.  For

```math
A={0,0},\qquad A'={a,-a},\qquad0<a<min(B,J/2),
```

the raw middle-fibre values differ while the endpoint-chord roofs agree, so
no future in the declared same-`J` language distinguishes them.  More
strongly, at mass `n` every bounded-field block collapses to this chord
exactly when `J>=4B/n` uniformly over the class.  The surviving exact state
is only the total field.

Thus the same microscopic family displays three distinct complexities under
one explicitly declared future algebra: an exact real slope multiset, a
sublinear-bit approximate histogram, and a one-number synchronized quotient.
The collapse is not an optimizer coincidence; it is a composable contextual
equivalence.

## Example 30: coherent clamps are stable without a kernel gauge or small reset

On the two-dimensional max-plus projective line, two all-finite kernels can
induce

```math
P_0(z)=clip(z,0,1),
\qquad P_delta(z)=clip(z,delta,1).
```

They are `delta/2` apart in the half-Hilbert metric and both are idempotent,
so that distance is unchanged at every positive depth.  Their kernel
difference has nonzero rectangular circulation and every power retains
image diameter near one half.  Thus the previously tempting completeness
claim “endpoint gauge or small full-image reset” is false for fixed coherent
continuations.

The counterexample does not justify an unconstrained third category: on the
paired orbit, both clamps become stationary after one step.  Its finite
certificate is the common idempotent relation.  By contrast,
`clip(z+delta,0,1)` is equally close in one step but drifts to order-one
distance after `Theta(1/delta)` repetitions.

For fresh adversarial residuals the picture is sharp again.  On an
`r`-coordinate selector cell, absence of a tangent-reset factor for `T`
steps lets an adversary align at least `floor(T/[r(r-1)])` residuals on one
output pair.  Hence depth-uniform stability is equivalent, up to the stated
factor, to syndetic rank-one selector products after endpoint gauges are
removed.  Coherent algebraic relations and adversarial robustness must not
be conflated.

## Example 31: one additive algebra has exact degree, conditioning, and lattice scales

Let a mass-`n` system contain `c_j` copies of each of `d` atom types, and let
type `j` contribute the bounded future-response function `phi_j`.  Equal
mass histograms are contextually equivalent precisely when their difference
lies in the integer kernel of

```math
z mapsto sum_j z_j phi_j.
```

The quotient therefore has `Theta(n^r_Z)` states, where `r_Z` is the integer
rank of the atom-response differences.  Histogram addition is the exact
future congruence; no dynamic-programming table was assumed.

For equally spaced heterogeneous mean-field bins, the atom functions are
hinges `(gamma_j+lambda)_+`.  Values at their knots recover tail sums of the
histogram.  This gives real conditioning at least `Delta/4` and literal
lattice margin exactly `Delta`, so below response error `Delta/2` the exact
external cover has

```math
{n+d-1 choose d-1}
```

members.  The example also separates arithmetic and robust information: a
single real query can encode rationally independent atom values with large
exact integer rank but zero real conditioning.  Query-parameter dimension by
itself is therefore not a state-growth law.

## Example 32: a composable zonotope sketch from atom response nets

For vectors `v_i` in the Euclidean unit ball of fixed dimension `p`, query
the directional support of all signed sums:

```math
F_V(theta)
=max_(epsilon_i in {+-1})<theta,sum_i epsilon_i v_i>
=sum_i|<theta,v_i>|.
```

One vector is therefore an additive response atom
`phi_v(theta)=|<theta,v>|`.  Uniformly over unit directions,

```math
||phi_v-phi_w||_infty<=||v-w||_2.
```

Quantize every vector once to one common root-scale Euclidean net and store
the type histogram.  A net of mesh `eta` has at most
`(1+2/eta)^p` cells, composition adds histograms exactly, and every support
query changes by at most `n eta`.  Taking `eta=n^(-a)`, `0<a<1/p`, gives
both `o(n)` bits and `o(n)` uniform error.

This validates atomic type quantization outside mean-field hinges.  It is an
upper theorem only: without exposure or lattice conditioning, different
histograms can have the same zonotope support response.

The finite verifier checks 24,000 exact signed-support identities, 24,000
root-scale quantization bounds, and 1,200 histogram merge identities.

## Example 33: adaptive refinement cannot recover a collapsed child type

Let atoms be `0` and `1`, each contributing that scalar on a single query.
A child which uses the one-centre net `{1/2}` at radius `1/2` stores only its
mass: all numbers of ones have the same coarse histogram.  If a parent later
switches to the exact net `{0,1}`, no update from those child states can
recover how many ones they contained.

Thus atomic type quantization is depth-stable only when one root-scale net
and quantizer are fixed throughout the merge tree, or when explicit
refinement data is retained.  A sequence of increasingly accurate static
nets is not automatically a composition congruence.

## Example 34: coordinate atoms make the type threshold sharp

Let there be `D` future queries and `D` atom types, with type `j` responding
one to query `j` and zero to every other query.  The response of a mass-`n`
system is exactly its histogram vector.  Distinct histograms are separated
by at least one in sup norm, so below error `1/2` the response cover contains
all

```math
{n+D-1 choose D-1}
```

states.  If `D/n` stays positive along an infinite subsequence, this costs
`Omega(n)` bits on that subsequence.
The sufficient condition `D(eta_n)=o(n)` for vanishing type-histogram rate is
therefore a sharp universal threshold; improvements require algebraic
dependencies among atom responses.

## Example 35: fractional consensus without a finite tangent reset

On `R^2/R1`, let

```math
P_alpha=((1-alpha,alpha),(alpha,1-alpha)),
\qquad0<alpha<1.
```

Its Dobrushin coefficient is `rho=|1-2alpha|`, so Theorem 16.18 gives exact
depth-`T` fresh-residual gain

```math
epsilon sum_(k=0)^(T-1)rho^k.
```

For `alpha ne 1/2`, no finite power has identical rows. Convex mixing can
therefore forget old response directions asymptotically rather than at one
finite reset. Merely touching a max-plus tie face is insufficient: the
segment can spend an arbitrarily small fraction of its length on one
selector, making the mixing coefficient arbitrarily weak.

## Example 36: different selectors force a cross-difference carrier

Let `x=y=(0,z)`, let the first channel use the identity selector, and let the
second swap its two coordinates. The old diagonal error is identically zero,
but the new diagonal error is `(z,-z)`. Thus diagonal error plus the selector
labels does not close. The joint array

```math
D_(ij)=y_i-x_j
```

does close under paired selectors, by (16.142), and retains the cancellation
needed for Hilbert response.

## Example 37: a local selector cycle need not be pumpable

The all-finite max-plus projective map

```math
T_delta(z)=clip(z-delta,0,1)
```

has a slope-one middle cell with `T_delta(C) intersect C` nonempty. The usual
local face graph therefore contains a self-loop. Nevertheless every orbit
loses `delta` per visit and leaves the cell after finitely many steps. A
nonzero witness cycle in such an over-approximation is not a drift
certificate. Whole-cell invariance or another exact path-lifting theorem is
needed before Theorem 16.19 has a converse.
