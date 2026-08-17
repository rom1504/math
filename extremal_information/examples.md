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

## Example 38: affine normal directions do not define an invariant fan

For `F(x)=x+1`, the linear normal of the hyperplane `x=0` is fixed.  But the
negative sign atom is not mapped into one atom: `F(-2)<0` while
`F(-1/2)>0`.  A switching certificate must retain the full affine form,
including its offset.

Even zero additive holonomy of normal directions is insufficient under
rescaling.  For `F(x)=2x`, pullback of `x-1` produces the distinct boundaries
`x=2^(-k)`.  The finite zero-holonomy closure law is valid for unit transport
of oriented normals, as occurs for coordinate-difference forms under
selector maps; it is false for arbitrary affine scaling.

## Example 39: a five-piece compact system has exponential future memory

On `[0,1]`, use the two encoders

```math
E_0(x)=x/3,
\qquad E_1(x)=(x+2)/3,
```

and the continuous three-piece decoder which maps both outer thirds back to
the whole interval and folds the middle third.  It satisfies `R E_b=id`.
After `t` encoders, the `2^t` ternary histories are pairwise separated by at
least `1/3` under one of the future observations
`h,h circ R,...,h circ R^(t-1)`, where `h(x)=x`.

Thus error below `1/6` requires at least `2^t` predictive states, although
the static and complete one-step response images have `O(1/epsilon)` covers.
The writing modes are contractions; the permitted expanding decoder is what
turns microscopic stored digits into macroscopic future response.  Removing
the decoder removes the explosion.

## Example 40: local tie compatibility has a sharp Helly failure

For cyclic coordinates modulo constants, set

```math
phi_i(x)=max{x_(i+1),x_i+c_i},
\qquad sum_i c_i=1.
```

Every proper subfamily of the `m` tie faces intersects, because deleting one
cycle edge leaves a forest of difference equations.  All `m` cannot
intersect, and the precise robust gap is

```math
min_x max_i|x_(i+1)-x_i-c_i|=1/m.
```

Pairwise, or even `(m-1)`-wise, optimizer-face compatibility therefore does
not certify a realizable switching word.  The obstruction can first appear
at arbitrarily large dimension.

## Example 41: strict width-three Ising exposes its whole boundary table

For a strip prefix, let `f(x)` be the best energy conditioned on its exposed
width-`w` spin column.  Arbitrary legal future columns give the exact
contextual metric

```math
D_ctx(f,g)=||f-g||_infinity.
```

A single strong ferromagnetic column plus fields pins any requested old
boundary coordinate.  Hence the usual transfer table is not merely a
convenient dynamic-programming state: it is the coarsest contextual state.

This is not an artifact of allowing arbitrary tables.  A strict width-three
nearest-neighbor strip has an explicit integer point whose reachable response
map has determinant `-1024`; its affine chart realizes a closed
eight-dimensional sup-cube of radius `1/2`.  There are genuinely eight
absolute and seven projective continuous response directions.

## Example 42: a restricted strip lumps exactly, while fresh rounding drifts

A two-letter rational width-two alphabet has seven reachable normalized
messages but only two weighted residual classes.  Partition refinement
discovers the exact quotient

```math
       c_0       c_1
A   (A,3)     (B,3)
B   (A,1)     (B,3).
```

Thus a transfer table can admit a strict path-realizing quotient when the
future alphabet is restricted.

By contrast, repeatedly quantizing an antiferromagnetic one-spin submodel
with coupling `-K` and field `s` rounds the projective shape back to flat at
every step but predicts scalar toll `K+s` each time.  The exact tolls
alternate `K+s,K-s`, so the error after `n` steps is `ns` for even `n` and
`(n-1)s` for odd `n`.  Projective state stability does not control the
absolute optimum unless the scalar reward cocycle is also compatible.

## Example 43: contraction and branching have an exact response-entropy law

Put an independent coordinate block at every node of a finite `q`-ary tree.
Let each input move to one child and multiply all surviving coordinates by
`rho<1`; observe only the current root.  The depth-`T` response tree is the
box

```math
prod_(k=0)^T[0,D rho^k]^(p q^k).
```

Its exact sup-covering number is the product of the one-dimensional covering
numbers.  It attains the contraction-weighted context-tree bound factor by
factor.  At fixed accuracy, contraction truncates relevant depth at
`Theta(log(1/epsilon))`, but branching turns that into
`Theta((1/epsilon)^(log q/log(1/rho)))` response states.  Mixing time changes
the scale at which static entropy is paid; it does not simply multiply a
static state count.

## Example 44: exact automaton state can jump while response error is continuous

A four-state, two-letter max-plus automaton has a coarsest strong block
quotient with two states.  Its quotient is discovered by block-signature
refinement, and actual depth-one suffixes expose both blocks.

Increase one microscopic maximizing self-loop by `delta>0`.  The exact strong
quotient now refines to all four singleton states, but the old two-state
quotient has one-step response defect only `delta`.  Repeating that self-loop
gives error exactly

```math
n delta
```

at every depth `n`; the nonexpansive sum-of-defects bound is sharp.  This
separates exact state complexity, which is discontinuous at `delta=0`, from
finite-horizon response distortion, which changes continuously.  The
repeatable defect cycle—not the mere failure of exact lumpability—is the
obstruction to depth-uniform reuse.

## Example 45: rational compact selector dynamics have finite symbolic memory

On a compact rational projective polytope, let every branch copy coordinates
and add rational offsets. Pulling back one affine observation wall can create
new offsets, but its normal only redistributes the finitely many labelled
coefficients among coordinate bins. The offsets lie on one rational lattice,
and only finitely many lattice translates cross the compact carrier.

Thus the complete contextual refinement is finite without enumerating orbit
words. This is an algebraic path-realization mechanism distinct from
contraction: the maps may retain old differences exactly. It preserves a
finite polyhedral observation coloring, not an arbitrary accumulated real
reward.

## Example 46: irrational rotation defeats every finite autonomous predictor

Let `F_alpha(x)=x+alpha mod 1` on the circle, with irrational `alpha`, and
observe the circle coordinate in geodesic distance. Any trajectory of a
finite autonomous predictor is eventually periodic. Along one eventual
phase, the true orbit advances by `p alpha` and is dense, so its distance
from the phase's fixed decoded value has supremum `1/2`. Consequently no
finite predictor has uniform error below `1/2`.

This remains true for badly approximable `alpha`: equicontinuity and strong
Diophantine regularity do not replace a finite invariant net. A nontrivial
arc coloring likewise has no finite exact predictor. The failure is temporal
phase drift, not numerical instability.

## Example 47: a two-centre cover can require linear dynamic memory

On the invariant grid `{0,alpha,...,nalpha}`, use

```math
F(x)=min(x+alpha,nalpha),
\qquad h(x)=x.
```

At error `eta=kalpha`, the static behavioral metric has internal covering
number `ceil((n+1)/(2k+1))`, but the exact minimum number of states in an
infinite-depth `eta`-predictor is `n-2k+1`. For `k=n/4`, two static centres
coexist with order-`n` reusable memory. A repeated predictor state too early
would be both close to its transient value and, on its eventual cycle, close
to the saturated value. The two requirements differ by more than `2eta`.

## Example 48: one contracting edge per control cycle gives the sharp gain

Take an `L`-cycle of metric fibres with `L-1` identity edges and one edge of
coefficient `theta<1`. Add an aligned fresh error `eta` on every transition.
Observed just before the contracting edge, one circuit sends

```math
z |-> theta z+Leta.
```

The limiting error is exactly `Leta/(1-theta)`. At `theta=1` it grows
linearly. This simultaneously proves the sharpness of the block-contraction
bound and shows why merely knowing that some local maps are contracting is
insufficient: every reachable cycle must pay contraction.

## Example 49: periodic tests miss a noninvertible reward diamond

On `[0,1]^2`, let

```math
A(x,y)=(0,0),
\qquad B(x,y)=(0,x),
\qquad r_A(x,y)=x,
\qquad r_B=0.
```

Every nonempty word has only the zero periodic point and all periodic rewards
vanish. Yet the coterminal paths `A` and `BA` have rewards `x` and zero, so no
single state potential can telescope both. This does not cause drift: the
diamond is transient and every total reward is at most one.

Lifting by the finite realized affine germs separates the two phenomena.
Transient path mismatches contribute a bounded endpoint error; a nonconstant
label on a recurrent lifted cycle is repeatable and grows linearly. At
per-step error `epsilon`, the packing number of that cycle's mean-reward
image at separation `2epsilon` lower-bounds the simulator state count.

## Example 50: exposed cycle response detects a toll repair, not exact refinement

In the perturbed four-state weighted automaton, the old two-block quotient
underprices one maximizing `A` self-loop by `delta`, so `A^n` has error
`n delta`. Raising the corresponding quotient self-loop toll by `delta`
removes every recurrent slope discrepancy. Even representatives give exact
microscopic lifts thereafter; an odd initial state pays at most the one-time
loss `2+delta`. Thus the corrected two-state quotient has zero cycle-response
distance and uniformly bounded all-word error although the coarsest exact
strong partition still has four states.

Arbitrary syntactic cycles would give a false obstruction: an unperturbed
microscopic `0->1->0` cycle has mean discrepancy `-3/2` from its quotient
loop but is never maximizing. Cycle response must be formed from the
path-realizing **exposed reward dynamics**, not every edge of a raw
presentation.

The Ising rounding example has the same diagnosis. Its exact two-step mean
is `K`, while the freshly rounded scheme emits mean `K+s`, causing slope
error `s`. Retuning the one-state toll to `K` leaves only a bounded parity
error at most `s`; two phase states are needed only for finer absolute or
exact response. Recurrent response entropy controls asymptotic slope memory,
while anchored simple-path residual controls the bounded transient scale.

## Example 51: an invariant grid survives switching without contraction

On `[0,1]`, use the two rational selector-PWA maps

```math
F_-(x)=max(0,x-1/3),
\qquad F_+(x)=min(1,x+1/3).
```

Both have slope-one cells, so neither gives strict global forgetting. For
every `N`, however, the grid

```math
C_N={0,1/(3N),...,1}
```

is exactly invariant under both maps. Starting from the nearest grid point,
nonexpansiveness keeps the error at most `1/(6N)` under every switch word and
at every depth. This is not repeated rounding: the chosen grid point follows
one genuine trajectory. Exact enumeration verifies 120,066 grid/word
inequalities for `N<=12`.

## Example 52: rational nonexpansiveness need not preserve a finite net

The map `F(x)=x/2` on `[0,1]` is rational, continuous, and nonexpansive. If a
finite set `C` satisfies `F(C) subseteq C`, however, every positive `x in C`
would force the infinitely many distinct points `x/2^j` into `C`. Thus the
only finite forward-invariant set is `{0}`, which is not an internal
`eta`-net below the interval's covering radius. The invariant-grid theorem
uses integer-selector preservation of one common lattice, not rational PWA
data alone.

## Example 53: three inputs expose an exponential selector orbit

On the weight-`floor(r/2)` binary slice, a cycle and an adjacent
transposition act by rational max-plus permutation selectors and generate the
whole symmetric group. A repeatable identity input observes one centered
coordinate. Every pair of distinct subsets can be moved by the same suffix
so that the probe rewards differ by one. Consequently exact and every
`epsilon<1/2` cumulative-response predictor needs
`binom(r,floor(r/2))` states, while one state suffices at error `1/2`.

This example survives compactness, rationality, unique active branches,
isometry, zero translation holonomy, and a fixed three-letter alphabet. The
full germ group has `r!` elements, but its stabilizer is invisible; the
minimal response state is exactly the exponential constant-weight orbit.

## Example 54: bounded-reward congruences do not form a quotient lattice

Use states `{I,A,B}`. Every `a` resets to `A`, every `b` resets to `B`, and
only the edge `(A,b)` has reward one. Both partitions

```math
{I,A}|{B}
\qquad\hbox{and}\qquad
{A}|{I,B}
```

admit scalar tolls with at most one total error. Their join has one block.
The words `a^n` and `b^n` force both one-state tolls to zero, but `(ba)^n`
from `A` earns `n`. There is therefore no unique coarsest bounded-error
congruence. Yet two raw starting states have same-word response difference at
most one, because their trajectories synchronize after the first input.
Asymptotic pairwise response equivalence misses quotient-created cycle
incidence.

## Example 55: a switching envelope forgets aligned words

Take two one-state max-plus alphabets with

```math
(A_a,A_b)=(1,0),
\qquad
(B_a,B_b)=(0,1).
```

Their maximum-over-letters envelopes are the same one-loop system of weight
one. Nevertheless the aligned response difference is `+1` on `a` and `-1`
on `b`. A critical graph of the switching envelope answers an adaptive
maximum over letters, not a common-word query. The finite-projective-
semigroup theorem retains the synchronized word action and detects both
cycle slopes.

## Example 56: one lattice-PWA map is a binary counter

Use `m` Boolean bits in dual-rail coordinates `(u_i,v_i)=(b_i,1-b_i)`.
Prefix minima of the `u` rails compute carries, prefix maxima of the `v`
rails compute their complements, and a continuous min/max lattice circuit
swaps exactly the rails whose bits must flip when adding one. The circuit has
`O(m)` shared gates, is additively homogeneous, and every PWA branch selects
input coordinates. It preserves the projective span-one polytope and has one
actual orbit of period `2^m`.

The centered most-significant-bit probe takes values `-1/2,+1/2` in a
primitive cyclic word. An identity query freezes a chosen phase, so any two
phases can be shifted to probe gap one and that gap can be repeated. Exact
and every sub-half-error predictor needs `2^m` states; one state works at
error `1/2`.

Without the identity query, a bounded reward around the sole periodic orbit
has one cycle mean and phase totals differ only by a bounded remainder. Thus
the exponential orbit is not automatically exponential response
information. The declared future query is what exposes and pumps it.

Replacing each dual rail by a constant-weight block digit strengthens the
same circuit mechanism. With growing block length and polynomially many
shared min/max gates, the genuinely realized exposed orbit has
`2^(r-o(r))` phases. The exponential lower bound is therefore not hidden in
an exponential list of cells, clauses, or transitions.

## Example 57: approximate reward quotienting contains graph coloring

Give every vertex of a graph one fixed raw state and every edge `uv` one
reward coordinate: `+1` at `u`, `-1` at `v`, and zero elsewhere. All dynamics
are the identity, so every partition is a transition congruence. Its optimal
per-step error is half the largest reward-coordinate range inside a block.
At threshold `1/2`, its blocks are therefore exactly independent sets and
the minimum quotient size is the graph's chromatic number.

The hardness is not caused by switching, orbit growth, or holonomy. It is
the simultaneous geometry of the declared query family. A rank-one-reset
variant shows that arbitrarily strong forgetting does not remove it.

## Example 58: exact word spectra can hide all rowwise path state

Let states and letters both be `[r]`. Under letter `e`, every edge leaving
state `e` has weight zero and every edge leaving another state has weight
`-C`. For any word, following its successive letter names gives a closed
zero-weight path, while no edge is positive. Every word therefore has
spectral radius zero, exactly matching a one-state zero system.

A coherent path lift tells a different story. Two raw rows in one fibre must
have block maxima within the sum of the upper and lower local defects. Letter
`i` separates row `i` from every other row by `C`, so any sub-`C` relational
certificate retains all `r` states. A diagonal gauge cannot repair both
directions simultaneously.

The scalar optimum uses a word-dependent critical witness. A path lift asks
one continuation relation to work from every raw representative. This
quantifier gap can be macroscopic even under exact response equality.

A two-letter de Bruijn version makes it exponential. Its states are binary
windows of length `m`; a letter shifts a matching first bit and appends either
bit. Every periodic word has a closed zero window path, so one scalar state
still gives every spectral response, but each window is the unique zero row
for its own length-`m` word. Exact coherent path lifting needs all `2^m`
states.

Blocking any `L>=2m` letters makes the separation sharper, not weaker. Each
block product has one all-zero row and all other rows equal to `-C`, so it is
projectively rank one, maximally contracting, and has a uniquely exposed
critical node. That node is the first `m`-window of the block word and ranges
over all `2^m` possibilities. Thus wordwise uniqueness and contraction do not
create a reusable state: they select a different hidden witness after each
future word is revealed.

## Example 59: width-two Ising separates anticipatory and forward state

On boundary spins `(s_1,s_2) in {+-1}^2`, three width-two Ising transfer
letters can be chosen so that their optimal edges factor as

```math
D_a times K_+,\qquad D_b times K_+,\qquad D_c times K_-,
```

where `D_a={s_2=1}`, `D_b={s_1s_2=1}`, `D_c={s_1s_2=-1}`, and
`K_+` or `K_-` fixes the first target spin. Each `K_q` meets every `D_e`, so the
future letter selects a valid microscopic predecessor and resets the
two-state support to `K_(tau(e))`. One support cannot be stable under both
`a` and `c`, while the four raw source signatures are distinct. Exact
unrooted response therefore costs two anticipatory states but exact forward
lifting costs four.

Adding an antiferromagnetic horizontal bond on `a`, with baseline penalty
`C>4`, makes the response order-sensitive without changing the state count:

```math
rho(T_w)=2N_(ca)^cyc(w)-N_a(w).
```

The words `aabccb` and `abbcac` have the same letter counts and responses
`-2,+2`, proving that no scalar per-letter toll suffices. The two-state
support remembers precisely the previous reset sign, not the discarded
second boundary spin.

## Example 60: deterministic de Bruijn reverses the carrier hierarchy

On words `I=E^m`, let input `e` deterministically shift and append `e`. Give
that one edge weight zero and every other edge weight `-C`. Every word map has
a periodic fixed point, so all word spectral responses are zero. A one-state
rowwise path lift is exact because every source has its deterministic tight
successor.

Backward-surjective support lifting behaves oppositely. After a length-`m`
word, every nonempty support maps to its singleton suffix. Exact anticipatory
support therefore needs all `q^m` singletons. With `N` support states, its
optimal certificate toll is `Theta(C/(1+log_q N))` before the exact threshold.
The exact checker verified 5,379 periodic fixed-point instances, 425
length-`m` singleton images, and 981 support/potential edges for `q=2,3` and
`m<=5`.

This is not semantic response complexity. It proves that source-total and
target-surjective proof carriers are incomparable, not successive levels of
one hierarchy.

## Example 61: rank one has `2^n` exact atoms but fixed-error compression

Put

```math
p_x={sum_i2^(i-1)x_i\over2^n-1},\qquad H(x)=-p_x^2.
```

The field `2p_x` uniquely exposes `x`, so the exact rank-one upper roof has
all `2^n` Boolean atoms. The margins are exponentially small. An exact
rational verifier checked all 510 exposures through `n=8`, together with 171
two-body response identities and 513 three-body associativity identities.

At fixed error, feature bucketing instead uses only
`(1+2PQ/epsilon)^r` cells. A code construction proves that
`2^(Omega(r))` bits are nevertheless necessary when rank grows. Exact atom
count and macroscopic response rate are therefore different scales.

## Example 62: a matching bridge has an extensive live interface

Choose an exponential Hamming code `C_n` and give every codeword an
independent binary energy bonus of size `delta n`. Under the degree-one
bridge `<x,y>`, query `y=c` uniquely exposes the bonus at `c`. The response
class contains `2^|C_n|` tables separated by `delta n`, even projectively
after restricting to constant-weight bonus labels.

The finite verifier checked 560 hostile/friendly pinned queries through order
eight. Bounded edge degree and treewidth of the isolated bridge do not imply
compression when the live boundary has `n` endpoints.

## Example 63: signed overlap holonomy costs `2n`

Take three balanced spin blocks with identity couplings `+,+,-` around the
triangle. Every separately optimized edge contributes `n`, but at each
coordinate the joint reward

```math
x^1x^2+x^2x^3-x^1x^3
```

is at most one. The true optimum is `n`, an extensive gap of `2n`. More
generally, an unbalanced signed cycle of length `ell` has optimum
`(ell-2)n` instead of `ell n`.

For a balanced signed graph, a vertex gauge makes all identity couplings
nonnegative and nested plus sets realize every pair optimum simultaneously.
Enumeration verified five heterogeneous multi-block instances and unbalanced
cycles of lengths three through five.

## Example 64: zero-temperature Viterbi mixing is all or nothing

For an all-finite max-plus transfer matrix, the global projective Lipschitz
coefficient is zero exactly for additive-rank-one matrices `K_ij=u_i+v_j`
and is one otherwise. Exact enumeration of all 625 integer `2 by 2` matrices
with entries in `[-2,2]` found 85 rank-one resets and an exact ratio-one
witness for each of the other 540 matrices.

Thus ordinary positivity and irreducibility do not reduce worst-future
best-path memory at zero temperature. A uniform `eta`-perturbation of rank one
does give an `eta`-accurate shape reset, providing the correct approximate
benchmark.

## Example 65: one local orbit probe makes de Bruijn phase observable

In the deterministic de Bruijn shift, let cyclic coordinate rotations act on
windows and read only the first symbol. The orbit of this one readout lists
all symbols of a window. Under a nontight-edge gap `C`, repeating the readout
`k<C` times separates all `q^m` length-`m` phases by at least
`k/(q-1)` in rooted response.

This is an intrinsic packing generated by two navigation letters and one
local probe, not by a raw-state lookup table. With finite leakage, however,
no common future can separate endpoint vectors by more than `C`; the phase
rate disappears when the probe reward overwhelms the filter gap. The example
simultaneously proves observability and its ceiling.

## Example 66: a dense sign bridge is incompressible at `n^(3/2)`

For a random sign matrix `B` and random queries `y_c`, put
`x_c=sign(By_c)`. With exponentially many queries at a sufficiently small
rate, concentration gives

```math
x_c^TBy_c>=d_0n^(3/2),
\qquad x_d^TBy_c<=d_1n^(3/2)\quad(d!=c),
```

with `d_1<d_0` and `||B||_(2->2)=O(sqrt n)`. Independent bonuses on the
`x_c` become independent exposed response coordinates. Uniform target-scale
accuracy therefore needs exponentially many bits.

The obstruction is not solely a programmable lookup landscape. The linear
children `h_c(x)=-x^TBy_c` have responses
`||B(y-y_c)||_1`; pairwise separated queries give an `Omega(n)`-bit
projective packing at the same scale.

## Example 67: the Potts shield has a clipped tropical state

For a `q`-state ferromagnetic Potts edge of strength `K`, every outgoing
normalized max-sum message lies in, and every point realizes,

```math
\{r in [-K,0]^q:max r=0\}.
```

Arbitrary unary futures expose every coordinate, so this is the coarsest
projective state. Its optimal contextual codebook has
`Theta_q((K/epsilon)^(q-1))` elements. If factors and couplings lie on one
lattice, the clipped carrier is exactly finite and closes at arbitrary tree
depth; rounding factors once avoids message-by-message drift.

## Example 68: vertex cover, not edge degree, controls sparse bridges

If a bipartite bridge has a vertex cover of size `k`, conditioning on its
cover spins yields a universal `2^k`-entry response table. A matching of `k`
edges makes this sharp: a bridge weight larger than the child-table spread
selects `x=y` and exposes an arbitrary table coordinate at every query.
The projective response rate is
`Theta(2^k log(1+D/epsilon))` bits.

The isolated matching has degree one and treewidth one, but its minimum
vertex cover and simultaneous live interface both have size `k`. Sparse
structure compresses only when it reduces the live cut, not when each edge is
locally simple.

## Example 69: planted quadratic poles survive a dense bridge

For every sign vector `z`, the complete sign quadratic

```math
H_z(x)=((x^Tz)^2-n)/2
```

has only the pole pair `+-z` as optimizer after adding any field with
`infinity` norm below `n/2`. A random sign bridge supplies exponentially
many query fields whose own sign poles have diagonal correlation
`Omega(n^(3/2))` and pairwise cross-correlation smaller by a fixed amount.
The resulting genuine sign-quadratic responses form an
`exp(Omega(n))` projective packing at `n^(3/2)` separation. Exact checks
cover 1,400 pole-locking instances; a seeded `n=32` certificate has absolute
and projective gaps 56 and 68.

This is the first dense-bridge lower bound inside the complete
`+-1`-coefficient quadratic class, but it deliberately exposes its own
limit: every child has cap `(n^2-n)/2`. It therefore proves an extensive
information rate for all sign quadratics, not for bounded-cap or
near-minimizing signings.

## Example 70: a derivative martingale compresses branching extremes

In a boundary-case branching random walk, the finite-depth derivative state
does not close by itself:

```math
Z_(r+m)=sum_(|u|=r)e^(-V(u))
             (Z_m^(u)+V(u)W_m^(u)).
```

At the extremal limit `W_m^(u)` vanishes. The remaining derivative mass
`Z_infinity` composes by weighted addition, and the full unmarked extremal
process is conditionally a decorated Cox process of intensity
`lambda Z_infinity e^x dx`. Every Laplace-functional query therefore sees
the realized environment only through one scalar.

This is a compression produced by critical renormalization, not by a finite
max-plus quotient. It is also sharply query-relative: if futures may inspect
genealogy labels, the state expands from the total mass to a derivative-mass
measure over branches.

## Example 71: one sparse mask list serves every signed feature landscape

Fix `m` bounded public features on a finite state set. Bernoulli thinning,
importance reweighting, and one union bound over the state rows show that a
single input has a uniformly accurate sparse surrogate with positive
probability. A second probabilistic argument needs only `m+1` public masks to
cover all `2^m` coefficient signings. The encoded state is one mask index and
the retained signs.

For Boolean quadratics at error `epsilon n^(3/2)`, this discards at least an
`epsilon^2/4` fraction of the `binom(n,2)` signs into a sparse `2`-bounded
surrogate. The same theorem sparsifies Boolean Littlewood polynomials,
bounded-CSP dictionaries, and codeword-correlation landscapes. Its limitation
is equally informative: the surrogate answers every shared future once, but
need not remain inside the original model or close under repeated internal
composition.

## Example 72: a cap-`1/2` Walsh family carries `sqrt(n)` bits

At orders `n=2^(2m)`, regularize the Walsh matrix by a self-dual bent vector,
hollow it, and switch it by the `2^(sqrt(n))` Maiorana--McFarland vectors

```math
s_g(u,v)=(-1)^(u dot v+g(v)).
```

Every child has the exact spectral cap `n^(3/2)/2`. The Walsh bridge maps a
dual Boolean query to the flat field `sqrt(n)s_g`. A low-bias code of
`exp(Omega(sqrt(n)))` tables keeps every off-diagonal product in a bounded
Rayleigh sector. A one-line resolvent bound then puts every cross response at
most `11n^(3/2)/8`, while the matched response is `3n^(3/2)/2`.

Thus the family has projective response separation `n^(3/2)/8` and needs
`Omega(sqrt(n))` bits. This is a genuine sub-landscape state—its `sqrt(n)`-
bit truth table is also sufficient—but it shows that even exact
conference-scale cap does not force bounded response complexity.

## Example 73: sparse near-top entropy becomes a linear response rate

Let `H` be any Boolean landscape whose states within `d_0n^(3/2)` of its
maximum occupy at most `exp((log2-kappa)n)` points. A random sign bridge and
an exponential almost-orthogonal query code turn coordinate switches of
`H` into an `exp(Omega(n))` projective response packing. The row event is a
weighted sign disagreement of probability arbitrarily close to `1/2`; its
relative-entropy exponent defeats the entire near-top set.

This applies at two very different cap scales. Complete quadratic poles have
only `exp(O(sqrt(n)log n))` fixed-deficit near-top states. More importantly,
the regular-Walsh child has exact cap `n^(3/2)/2`, and Hanson--Wright supplies
the required entropy gap. Its full switching family therefore has
`Theta(n)` contextual response bits at fixed `n^(3/2)` error: `Omega(n)` by
packing and `O(n)` by storing the switch.

The theorem does not make extremal entropy a sufficient state. It makes an
entropy *deficit* a certificate that hidden rooting information becomes
observable under a suitable low-operator-norm interaction.

## Example 74: Walsh composition exposes a commutation cocycle

For one Walsh child, bias and a single pair-product Rayleigh coordinate
control the explicit response separation. They do not survive repeated
composition. Take a linear truth table `g_a(v)=a dot v`. The normalized child
involution `C_a` and Walsh bridge involution `F` satisfy

```math
FC_a=(-1)^(a dot a)C_aF.
```

Choose one nonzero even-weight `a_0` and one odd-weight `a_1`. Constant words
of either label have identical lists of individual biases and within-word
pair correlations. On a bipartite `k`-block graph, however, the commuting
word simultaneously saturates every child and edge, while the anticommuting
word is bounded by the spectrum of `I+A(G)^2`. On a long path their optima
differ by

```math
\left({3-\sqrt5\over2}+o(1)\right)k n^{3/2}.
```

Thus composition creates an extensive observable from a relative
commutation phase invisible to the retained scalar summaries. The exact
Kronecker presentation uses `k sqrt(n)` truth-table bits and is closed under
graph gluing, but it is only coefficient-level closure; it does not itself
compute the Boolean maximum.

## Example 75: permutation bent switches gain a logarithmic factor

Replace the Boolean truth table in the first Walsh construction by a
permutation `pi` of `F_2^m` and use

```math
s_pi(u,v)=(-1)^(u dot pi(v)).
```

All `q!` such vectors remain bent. For two permutations, the rooted Rayleigh
coordinate is the signed bias with which
`tau(v)=v+pi(v)+sigma(v)` preserves the Boolean inner product. If that bias
exceeds `1/4`, Fourier Parseval forces `m-O(1)` independent large characters.
Their joint Hoeffding cost has speed `q log q`, enough to survive conditioning
two random functions to be permutations and then applying Turan's bound.

The resulting explicit-family code has `exp(Omega(q log q))` children at the
same `n^(3/2)/8` response gap, where `n=q^2`. Thus its response information is
`Theta(sqrt(n)log n)`, matched by listing the permutation. The construction
is weaker than the probabilistic full-orbit linear-rate theorem but exposes a
useful algebraic invariant: approximate self-isometry of the Walsh bilinear
form.

## Example 76: discrete adversarial chains have a lower-spectral limit

For a finite-state nearest-neighbor chain, keep the local disorder alphabet
finite and adversarial rather than random or convexified. At inverse
temperature `beta`, every disorder letter is a positive transfer matrix.
The minimum partition function differs by only a fixed factor from the
minimum norm of a matrix product, so its pressure is the log lower spectral
radius. The uniform soft-max sandwich then proves existence of the
zero-temperature adversarial ground-state density with error at most
`log(q)/beta`.

When the transfer entries are uniformly positive, projective cavity updates
contract in Hilbert metric. A `delta`-net of the `(q-1)`-dimensional cavity
simplex has `O(delta^{-(q-1)})` states and approximates the asymptotic mean
pressure within `delta/(1-kappa)`. This is a genuine nonconvex-disorder limit
and a benchmark success for response plus contraction.

It also diagnoses its own transfer ceiling. Every dense sign bridge has
Boolean bilinear maximum at least `n^(3/2)/sqrt(3)`, so a balanced dense cut
already carries a leading interface and the standard transfer state has
`2^n` coordinates. Fixed-width lower-spectral multiplication cannot simply
be transplanted there without a new nonlocal quotient.

## Example 77: relative involution algebra is robust but conditional

For a symmetric involutory child `C` and bridge `F`, the anticommutator gives
an exact spectral ceiling on every graph composition. On a bipartite graph,
a Boolean child pole transported to another Boolean vector by `F` gives the
opposite certificate: its failure to remain a child pole is exactly

```math
{1\over2}||[C,F]s||_2^2.
```

Thus a fixed anticommutator gap below two forces an extensive separation
from an exactly commuting family on every positive-degree regular bipartite
graph. The result robustifies the Walsh parity example without pretending
that two operator norms form a complete quotient. The Boolean transported
pole is an indispensable realizability hypothesis.

## Example 78: a binary Gram collision changes a Walsh path at leading scale

Linear Walsh labels can have the same complete binary Gram matrix while
different hidden linear relations change the simultaneous maximizing
section. On three blocks, `(a,a,a)` and `(a,b,a)` can have all Gram entries
one, yet the second word exactly attains `7n^(3/2)/2` while the first is at
most `3sqrt(3)n^(3/2)/2`.

This is a deliberately narrow falsifier. Ordinary truth-table overlaps see
whether two linear labels are equal and therefore distinguish the example.
What fails is the tempting Gram/commutation quotient used alone.

## Example 79: a rooted relation form is an exact Walsh orbit state

In every label dimension, the binary Gram form and relation kernel can omit
whether a label combination equals the characteristic vector
`omega=(1,...,1)`. A single Walsh-rooted future sees that omitted bit with an
`n^(3/2)/6` response gap.

Adding the rooted relation coset closes the classification in both parities:
Gram, kernel, and `omega`-coset determine an ordered label tuple up to the full orthogonal
group. The resulting `O(k^2)`-bit state is an exact Boolean-extremal quotient
for every graph on `k` Walsh blocks, independently of the ambient label
dimension. It is nevertheless not separately reusable under gluing, because
cross-Gram values and mixed relations are newly created at the interface.

## Example 80: Walsh gluing is a rooted bilinear amalgam

Two independently summarized linear-label pieces become exactly composable
after supplying three relative objects: their cross bilinear form, the graph
of coincidences between their presented spans, and the fibre of pairs that
sum to the characteristic root. These data reconstruct the combined Gram,
relation, and root state by linear algebra and compose associatively on the
accumulated presented span.

The price is real. Fixed isolated states realize every one of `2^(rs)` cross
forms, and a separate family realizes `|GL(r,2)|=2^(r^2+O(1))` different
intersection correspondences with zero cross form. Worse, all singleton and
pair amalgams can agree while a ternary relation differs. Dynamic memory is
therefore the accumulated span, not an edge-local compatibility label. The
state remains strict--`O(t^2)` rather than `mt` bits when `t=o(m)`--but it
does not compose for free.

## Example 81: numerical rank depends on the extremal scale

At balanced size `n`, deleting a bridge singular direction of size at most
`epsilon sqrt(n)` changes every Boolean interaction by at most
`epsilon n^(3/2)`. The retained singular features therefore give an exact
low-rank roof plus a uniform spectral-tail error. For a bounded number of
feature-visible ports, a net uses
`exp(O(r_epsilon log(1/epsilon)))` cells.

This puts three bridge regimes on one scale. The identity part of
`alpha I+beta J` and every bounded-degree bipartite bridge are subleading at
`n^(3/2)`; a typical dense sign bridge has extensive numerical rank. The law
produces local factor tables, not a global state for a graph with a large
separator, and arbitrary state-specific futures can still inspect states
discarded inside one feature fibre.

## Example 82: a microcanonical hypograph retains finite-rate rare branches

At one exponential speed, store the upper-semicontinuous log-count profile
over a compact descriptor. Cartesian composition becomes supremal
convolution whenever the number of decompositions is subexponential, and
hypograph recovery preserves every finite-rate branch. The nonconcave
microcanonical entropy of the mean-field BEG model is the classical
benchmark: bounded canonical temperature data sees only its concave envelope.

Two limits are sharp. Exponentially many distinct descriptor fibres can add a
missing mass exponent unless descriptor complexity is charged, and one
maximum versus `e^(sqrt(n))` maxima is invisible at speed `n`. Thus this is a
leading-rate rare-event state, not an extremal point process or an automatic
structured realization theorem.

## Example 83: complete Walsh spectra can miss a rooted response bit

Every word in a linear-label Walsh child `J_a` reduces to a Weyl operator
whose trace is nonzero only when an even-length label sum vanishes. Its phase
is determined by binary Gram data. Consequently the Gram matrix and relation
kernel determine the complete spectrum of every weighted graph carrier on a
marked tuple; the characteristic-root fibre is spectrally invisible.

For synchronized copies of the rooted-collision children, this remains true
across every weighted graph experiment, while a canonical rooted Boolean
future has an `n^(3/2)/6` projective gap. Spectral completeness is therefore
strictly query-relative. A new appended label can add cross-Gram information,
so no claim is made for arbitrary unmatched contexts.

## Example 84: a larger ambient symmetry deletes the Walsh root exactly

Embed every linear label as `(0,a)` in the full `2m`-dimensional coordinate
space of the order-`2^(2m)` Walsh bridge. The ambient characteristic vector
is `(omega,omega)` and never belongs to the embedded label span. Therefore
two tuples with equal Gram form and relation kernel extend to one ambient
orthogonal coordinate permutation, which conjugates every child and bridge.

For odd `m`, constant words `omega^k` and `e_1^k` have different label-space
root fibres but identical entire landscapes on every unrooted real weighted
graph; the even-dimensional analogue uses `e_1+e_2`. A fixed rooted pole
breaks that larger symmetry and separates them. The root bit is not hidden
unrooted information--it is information created by declaring a coordinate
anchor.

## Example 85: local coincidence memory has a sharp linear rate

On `h` disjoint three-coordinate chunks, fix duplicated endpoint labels
`a_i=111` and choose each middle label independently as `a_i` or `c_i=001`.
All isolated states, cross-Gram data, and root fibres remain fixed. Only the
coincidence correspondence records the `h` choices. The `i`th ordinary
three-block path maximum changes by at least

```math
{7-3sqrt3\over2}n^{3/2}
```

when its middle label changes. Thus bounded connected supports admit an
`O(h)` exact local carrier, but no `o(h)`-bit fixed-error summary. Localizing
the orbit charts removes quadratic global compatibility without making the
remaining coincidence information free.

## Example 86: interaction mass prices approximate Walsh compatibility

For a public partition of a weighted linear-label Walsh graph, delete the
edges crossing its parts.  Each retained component is decoded exactly from
its local Gram and relation state, while the deleted interaction changes the
upper or absolute optimum by at most

```math
\left(\sum_{e\text{ crossing}}|w_e|\right)n^{3/2}.
```

Thus a unit path on `t` labelled blocks has an `O(t/eta)`-bit local carrier
with additive error `eta t n^(3/2)`.  The same architecture cannot compress a
complete graph with all edge magnitudes at least `c`: error at most
`eta t^2 n^(3/2)` forces `sum_C |C|^2 >= (1-2eta/c)t^2`.  Static orbit
compression and dynamic locality are separate resources, and interaction
mass is the exact exchange rate for this particular approximation.

## Example 87: an unrooted triangle sees binary Gram flux

Take three even linear Walsh labels with the sole relation `a+b+c=0`.  One
tuple can span a totally isotropic plane, while another has all three
off-diagonal pairings equal to one.  Their self-pairings, relation kernels,
and characteristic-root fibres agree.  Nevertheless the ordinary unweighted
triangle maximum is respectively

```math
{9\over2}n^{3/2}
\quad\hbox{and at most}\quad
{3(1+\sqrt{17})\over4}n^{3/2}.
```

The gap `3(5-sqrt(17))n^(3/2)/4` needs no root, field, or pinned spin.  Hence
the triangle flux in the off-diagonal Gram form is genuinely scalar-visible;
it cannot be discarded merely because the root fibre disappears in the
unrooted language.

## Example 88: regular Hadamard amplification closes a dense response state

Let `H` be a symmetric Hadamard matrix with a positive regular Boolean
eigenvector and amplify a fixed outer quadratic template by `H^(tensor r)`.
Tensoring every Boolean block witness with that eigenvector embeds its whole
normalized cross-correlation matrix at the next depth. The resulting compact
convex response sets are nested in fixed dimension and converge.

For `H=W_4`, a trace-zero sign template yields genuine hollow dense signings
of orders `d4^r`, and their normalized Boolean maxima converge. This is a
near-original thermodynamic-limit benchmark whose state is exponentially
smaller than the Boolean landscape. Its limitation is precise: exact tensor
amplification supplies the recovery map, while arbitrary orders and
non-tensor perturbations remain untreated.

## Example 89: second-speed entropy needs saddle mass

For fixed finite descriptors, a tower of separated logarithmic speeds is an
exact lexicographic tropical response algebra. At speeds `(n,sqrt(n))`, for
example, it distinguishes one extremizer from `exp(theta sqrt(n))`
extremizers and composes by lexicographic supremal convolution.

Growing decomposition fibres break the bare pointwise state. In
Vandermonde's identity, each central summand contributes second-speed
coefficient `-1` at scales `(n,log n)`, whereas their `Theta(sqrt(n))`
near-saddle multiplicity changes the exact coefficient to `-1/2`. The
missing information is tangent counting mass, not another value of the
pointwise entropy profile.

## Example 90: summable one-sided Boolean recovery tolerates non-tensor edits

A common Boolean lift that preserves every pair response within errors
`epsilon_r` makes the fixed-port convex carriers approximately nested. If
the errors are summable, compactness alone forces a Hausdorff limit; no
reverse projection and no consistent optimizer are needed.

This is realizable inside dense sign matrices. After each regular-Hadamard
amplification, flip all signs on a perfect matching. The level has
`Theta(N_r)` non-tensor edge edits, but their perturbation operator norm is
exactly two, so the normalized defects `2/sqrt(N_(r+1))` are summable.
Hollowing changes the normalized objective by only `O(N_r^(-1/2))`.

The summability scale cannot formally be weakened to square summability:
the scalar kernels `sin(log(r+2))` have square-summable step changes and two
separated subsequential limits.

## Example 91: one connected Walsh language carries many flux bits

Use `h` disjoint relation triangles of even Walsh labels.  In each triangle,
hold the relation `a+b+c=0`, all self-pairings, the root fibre, and every
cross-triangle pairing fixed, while choosing its three off-diagonal Gram
values to be either all zero or all one.  This gives an `h`-cube of states.

For each coordinate, its ordinary triangle maximum detects that flux bit.
Nonnegative Walsh bridges may connect all `3h` ports without reducing the
gap: the favorable common eigenvector saturates every connector, while the
operator norm gives the same connector ceiling in every competing state.
Activating all nontarget triangles with total small weight leaves a fixed
gap.  Hence the same bounded-degree-four connected support, or the same
complete support, carries an `h`-bit scalar response packing.

This closes scalar visibility at fixed one-block accuracy.  It does not
give an extensive whole-system rate, because the total-variable normalized
gap decays as `h^(-3/2)`.

## Example 92: an automatic dense signing has a continuous scale phase

Order the tensor powers of the regular order-four Walsh matrix so that each
is the leading block of the next.  Hollow the prefixes of the resulting
infinite sign matrix.  Every fixed base-four prefix phase is a fixed outer
template tensored with deeper Walsh powers, so its normalized Boolean
maximum converges.  Operator control of the vertices added between nearby
prefixes promotes these pointwise limits to one continuous mantissa profile.

The profile is not constant.  Full powers have ratio `1/2`, while an exact
48-spin certificate, amplified by the regular Boolean eigenvector, gives
ratio at least `89/(96sqrt(3))` on every `3*4^r` prefix.  Thus one explicit
coherent all-order dense signing sequence genuinely fails to converge even
though every fixed phase does.

This is not nonconvergence of the minimizing sequence.  It is the sharp
falsifier for the claim that cross-scale recovery alone synchronizes outer
phases.

## Example 93: residual synchronization and witness support are different

If every length-`D` max-plus product has row rank one, its normalized row is
one of finitely many suffix residuals.  One letter acts on these profiles by
a deterministic state transition plus a scalar toll.  Context-cycle tolls
then compute exact spectral response, and a violated potential budget pumps
linearly.

That residual quotient need not carry one locally thresholded witness
support.  The two rank-one matrices

```math
\begin{pmatrix}0&-1\\-2&-3\end{pmatrix},
\qquad
\begin{pmatrix}-2&-3\\1&0\end{pmatrix}
```

share the residual `(0,-1)` and give every word spectral radius zero, but
their zero-threshold relations have no common nonempty invariant support.
Thus an empty support core falsifies one proof presentation, not the scalar
response.  Arbitrary terminal pins, by contrast, expose projective residual
profiles at exact distance half their coordinate oscillation.

## Example 94: balanced Boolean pullbacks erase a scale phase

Let a finite or compact family of normalized dense quadratic operators carry
a phase label.  If convex combinations of signed coordinate replications
pull each phase into the next level with small operator defect, their phase
marginals define a Markov kernel.  A delay window that refreshes an
`alpha_r` fraction onto one full-support phase law forces the limiting
same-spin response to be constant whenever recovery plus transfer error is
`o(alpha_r)`.

Neither ingredient can be omitted.  A two-state scalar family with refresh
mass `2^(-r)` and transfer defect `2^(-r-1)` retains response phases zero and
one forever.  Applied to the nonconstant Walsh prefix profile, the theorem
forces every purported balanced-reordering certificate to pay at least
`89/(48sqrt(3))-1.01` defect per unit refreshed mass.

## Example 95: Gaussian tangent profiles form a finite response semigroup

A lattice Gaussian landscape is determined asymptotically by leading height,
mean, covariance, logarithmic power, and total-mass amplitude.  Under
convolution, heights and means add, covariances add, logarithmic powers gain
`d/2`, and total-mass amplitudes multiply.  This finite parameter tuple is
associative and recovers every future convolution without storing the
growing descriptor grid.

The same saddle-mass rule repairs Vandermonde's missing `1/2 log n` and the
`(q-1)/2` multinomial correction.  It is not universal: a quartic saddle has
`n^(3/4)` tangent mass.  Truncation, exponential rescaling, and flooring give
finite integer landscapes realizing the Gaussian semigroup at every large
order.

## Example 96: public bridges cannot amplify local hidden Walsh state

Put hidden labels in `k` Walsh blocks of order `n`, and allow an arbitrary
state-independent scalar bridge landscape.  If hidden state enters only
through bounded onsite child terms, two entire landscapes differ pointwise
by at most `Dkn^(3/2)=DN^(3/2)/sqrt(k)`.  Maximization cannot amplify this
difference.  Signed, dense, or enormous public connector baselines are
irrelevant because they cancel before optimization.

For the connected triangle-flux cube this proves the earlier one-port gap has
the correct `h^(-3/2)` full-cube scaling and that every subcode has
total-scale diameter `O(h^(-1/2))`.  Unequal disjoint cell sizes do not help:
at fixed total-scale distortion, their response packing has bounded entropy.
Any escape must make each local bit alter `Omega(k^(3/2))` unit interaction
atoms, or otherwise leave the bounded state-local scalar architecture.

## Example 97: logarithmic averaging forgets scale order

An automatic response sequence with phase profile `L(t)` has a canonical
logarithmic empirical law: push `dt/(t log h)` through `L`.  This remains
true when the pointwise sequence oscillates.  It is a genuine compression,
but not a reusable last-scale response state.  The distinct profiles
`cos(2pi log_h t)` and `cos(4pi log_h t)` have the same logarithmic law.

Every positive power-weighted mean retains a phase function `C_alpha(s)`,
and the exact inverse

```math
L(s)=C_alpha(s)+(s/alpha)C_alpha'(s)
```

recovers the whole ordered profile.  Thus a unique averaged thermodynamic
statistic can coexist with future contexts that recover all the information
it discarded.

## Example 98: perfect forgetting can hide an arbitrary reward algebra

An approximate row-rank-one suffix controls every normalized terminal query
with an error paid once, independent of the preceding depth.  It does not
control accumulated max-plus reward.  Two exact rank-one binary transfer
matrices already have contraction zero but incur an optimal one-state reward
error `delta/4` per letter because their directed endpoint compatibilities
have a nonzero cycle defect.

The obstruction is universal below the shell scale.  Multiplying any
bounded all-finite weighted automaton by `alpha` puts every product in a
single projective row shell of radius `alpha/2`, while preserving its full
spectral response algebra up to that scale.  Static residual proximity and
dynamic scalar sufficiency are therefore incomparable without a compatible
reward cocycle.

## Example 99: gauge refresh is cheap; semantic refresh is not

Signed-permutation conjugates of one dense operator can be two units apart
in fixed-coordinate operator norm while carrying identical Boolean response.
Their phase label is pure gauge, and arbitrary exact Markov refresh is
implemented by conjugating the same pullback.  One bit suffices.

For an observable response phase, a common-stationary expander has a
different cost.  If one phase exceeds the stationary average by `D`, then a
scrambling window must accumulate that excess as transfer toll.  At one-step
toll `delta`, the reusable state requires `Omega(D/delta)` bits.  Applied to
the Walsh prefix obstruction at defect `O(N^(-1/2))`, this excludes every
`exp(o(sqrt N))` half-scrambling phase quotient.

## Example 100: bounded-operator signs have linear visible SVD rank

An `n`-square sign bridge has Frobenius mass exactly `n^2`.  If its operator
norm is at most `C sqrt(n)`, then a fixed positive fraction of its singular
values must exceed every fixed threshold below `sqrt(n)`.  Thus the
scale-sensitive SVD roof cannot stay subextensive at the `n^(3/2)` response
scale.

This is not general incompressibility.  A full-rank orbit may collapse by
rearrangement, and the rank-one all-ones matrix collapses to magnetization.
The example separates spectral interface dimension from nonlinear response
congruence.

## Example 101: a closed rate roof can have an unclosed tangent law

Power costs `a|x|^p` close exactly under infimal convolution for every
`p>1`.  Their central lattice saddles have mass `n^(1-1/p)`, but the
normalized tangent density is `e^(-a|x|^p)`.  Its self-convolution stays in
the same scale family only for `p=2`: otherwise dyadic stability would
contradict the central limit theorem.

Hence the leading extremal roof can have a finite algebra while the next
composition query generates a growing function-valued state.  Gaussian
tangent closure is a rigidity phenomenon, not a generic consequence of
convex homogeneous rates.

## Example 102: compatibility broadcast restores a total-scale rate

Encode `Theta(k)` hidden bits as an alternating bilinear form on a public
`F_2` label list.  Every nonzero form changes a constant fraction of the
`k^2` edge phases.  A single random public sign dressing makes all resulting
quadratics spectrally flat, while a support-only discrepancy bound separates
every pair by `Omega(k^(3/2))` under a predeclared negative-clone context.

Thus state-dependent dense broadcast escapes the local-state ceiling and
forces `Omega(k)` response bits at total scale.  The price is explicit: a
quadratic-size shared public base, dense bit incidence, and same-support
additive contexts whose overlay is not itself an exact signing.

## Example 103: bounded fan-in has exactly `t sqrt(N)` total-scale capacity

For `O(N^2)` bounded quadratic atoms, if each atom reads at most `t` hidden
bits, the sum of all neighbouring contextual gaps is at most the total
bit--atom incidence.  Hence only `O(t sqrt N)` hidden coordinates can each
have an `N^(3/2)` gap.  This is sharp even at the information level: a
cellwise switching library followed by an outer code gives
`exp(Omega(t sqrt N))` pairwise-separated exact signings.

The example deliberately separates resources.  The sharp unrestricted
family need not be flat and may use an exponential public library.  The
flat alternating-form family is necessarily different: in every coordinate
system its average edge fan-in is linear and every hidden bit changes a
quadratic number of edge phases.  Dense broadcast, rather than a clever
repackaging of local scalar channels, is the source of its macroscopic rate.
