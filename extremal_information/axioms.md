# Candidate definitions and surviving principles

These are working definitions forced by examples.  They are not axioms of an
established theory, and they should be renamed or discarded when a sharper
object appears.

## 1. Declare the experiment before the summary

A finite **landscape experiment** consists of

```math
(\Omega,H,\phi,\Theta),
```

where `Omega` is the state space, `H:Omega -> R` is a normalized energy,
`phi:Omega -> R^d` is the feature visible to the environment, and `Theta` is
the permitted set of linear interventions.  Its response is

```math
V_H(\theta)=
\max_{x\in\Omega}\{H(x)+\langle\theta,\phi(x)\rangle\}.
```

The normalization and the query family are mathematical data.  A summary
that preserves `V_H(0)` need not preserve `V_H(theta)`, and a state sufficient
for uniform fields need not be sufficient for a labeled bridge.

Two experiments are **exactly response-equivalent on `Theta`** when their
response functions agree on `Theta`.  For approximate questions use

```math
d_\Theta(H,G)=
\sup_{\theta\in\Theta}|V_H(\theta)-V_G(\theta)|.
```

This is generally a pseudometric: the theory should identify landscapes that
no declared query can distinguish.

## 2. The exact quotient for a linear feature interface

For `u` in `conv(phi(Omega))`, define the upper response roof

```math
\widehat H_\phi(u)=
\max\left\{
\mathbb E_\lambda H(X):
\mathbb E_\lambda\phi(X)=u
\right\}.
```

Equivalently, retain the upper faces of the lifted response body

```math
K_\phi(H)=
\operatorname{conv}\{(\phi(x),H(x)):x\in\Omega\}.
```

For the full query family `Theta=R^d`, convex duality proves that the response
function and the roof determine one another.  Hence the roof is the minimal
exact quotient for this experiment.  For a restricted `Theta`, only the
support values in those directions are operationally minimal; the whole roof
may retain unnecessary faces.

This restriction caveat is important.  “The response roof is minimal” is not
valid without specifying a determining query family.

## 3. Extremal rate--distortion

Let `H_n` be a declared class and let a summary `Z=S(H)` be decoded into a
response function.  The deterministic uniform information cost is

```math
R_n^{\rm det}(\epsilon;\Theta)
=\inf_S\left\{
\log_2|\operatorname{range}S|:
\sup_{H\in\mathcal H_n}d_\Theta(H,\widehat V_{S(H)})
\le\epsilon
\right\}.
```

Under a prior `Pi_n` and loss `ell(H,Z)`, define the Shannon version

```math
R_{\Pi_n}(D)=
\inf_{P_{Z|H},\widehat V:
\mathbb E\ell(H,Z)\le D} I(H;Z).
```

Candidate asymptotic rates must specify whether information is normalized by
`n`, by the number of interactions, or by another natural model size.  A
landscape with `Theta(n^2)` independent couplings and a Boolean optimizer with
`n` labels live at different rates.

For deterministic uniform distortion, covering and packing numbers in
`d_Theta` give the exact elementary bounds.  For ensembles, source
rate--distortion, Fano inequalities, and posterior entropy give lower bounds.

## 4. Composition must act on the state

For additive composition with a common feature space,

```math
H_\oplus(x,y)=H_1(x)+H_2(y),
\qquad
\phi_\oplus(x,y)=\phi_1(x)+\phi_2(y),
```

the roofs compose by sup-convolution.  For one bilinear coupling

```math
H_B(x,y)=H_1(x)+H_2(y)+\phi_1(x)^TB\phi_2(y),
```

the two child roofs determine the exact parent maximum.

That does **not** imply closure under iteration.  A later coupling may inspect
a joint feature of `(x,y)` not determined by the two marginal feature means.
A candidate state passes the composition test only if one of the following is
proved:

1. the feature class is closed under the operation;
2. the state is enlarged by a controlled tensor/feature algebra whose metric
   entropy remains sub-landscape scale; or
3. a synchronization theorem makes every newly created feature a function of
   the retained state.

Merely writing down the infinite hierarchy is not compression.

## 5. Support resolution is part of the topology

For a normalized sequence, define upper-tail complexity

```math
\Sigma_H^\uparrow(e)=
\liminf_{n\to\infty}{1\over n}
\log\#\{x:H_n(x)\ge e\}.
```

One must distinguish:

- an empty level (`log 0=-infinity`);
- a subexponential but nonempty level (rate zero); and
- a positive-rate cloud.

For homogeneous quadratic forms, a theorem in [`theorems.md`](theorems.md)
shows that the maximum is the closure of the positive-rate upper tail.  For a
general landscape this need not hold.  Any compactness notion must state
which edge resolution it preserves.

## 6. Rooted and unrooted information are different experiments

An energy histogram and global pair-overlap law average over the roots and
coordinate locations.  A restriction, puncture, inserted vertex, block
bridge, or external field singles out an apparatus.  Its response depends on
a rooted joint profile such as

```math
(H(x),\phi_{\rm root}(x))
```

or its upper response roof.

The code and Curie examples prove that exact unrooted pair data can coexist
with a leading rooted response gap.  Therefore rooted data should be added
only when the future query forces it, not hidden inside an ambiguous phrase
such as “overlap geometry.”

## 7. Compactness and realization are separate obligations

At fixed feature dimension and bounded energy, downward convex response
bodies are compact in Hausdorff distance.  Every such body is approximable by
an unrestricted finite landscape, and bounded response values are continuous.

This does not imply realization inside a constrained model such as complete
sign matrices, linear codes, or a fixed CSP.  A useful limiting theory needs
both:

1. compactness of its abstract extremal state; and
2. a model-specific recovery theorem at the required finite sizes.

Failure of the second is not repaired by strengthening the state until it
contains a finite optimizer.

## 8. Current survival table

| Candidate | Verdict | Reason |
|---|---|---|
| Scalar energy entropy | survives for recovering a homogeneous quadratic maximum; fails composition | Hamming noise thickens the edge, but labels are absent |
| Full global pair-overlap support | exact for every global two-replica query built from energy and total overlap; fails labeled queries | fixed-half example |
| Finite `k` global overlap hierarchy | rejected as a universal invariant | for every fixed `k`, parity-half-cube code pairs agree through `k` replicas but have a scalable covering-radius gap |
| Upper response roof for a fixed interface | survives | exact duality and one-step composition theorems |
| Full-spin response roof | rejected as compression | cube vertices recover every `H(x)` |
| Query metric entropy | survives as information definition | exact packing/covering theorem and quadratic-rate example |
| Posterior response width `Gamma(R)` | survives as a fixed-embedding rate certificate | gives a sharp mutual-information curve, but same-space cancellation destroys composition |
| Bounded-dimensional universal state | rejected as a default assumption | scalable planted resonances evade fixed-arity tests |
| Boundary response kernel | survives for bounded separators | coarsest exact quotient for all endpoint fields; its universal worst-case cost is quadratic in the boundary-state count |
| Syndrome coset-leader profile | survives for fixed labeled binary syndrome interfaces | exact min-plus/union algebra, strict quotient of the code, sharp exact cost `Theta(2^w)` bits, and a positive macroscopic linear rate on a block family |
| Root-averaged outer code spectrum | rejected for appended-fragment composition | equal outer polynomials can have different response to the same labeled future fragment |
| Average conditional overlap variance | rejected for zero-temperature synchronization | rare exposed fibres can retain a fixed response gap |
| Linked deterministic overlap profile | survives conditionally | mixture ultrametricity plus a checkable cross-root path condition gives uniform scalar synchronization |
| Convex reachable body | survives approximately in fixed effective dimension | Shapley--Folkman makes nonconvexity cost at most the largest `r` component diameters; growing `r` can carry a leading gap |
| Robust tropical crossing rank | survives at uniform lattice scale | a four-cell gap protects channel count, but can disappear in normalized mean-square loss |
| Query-weighted tropical exposure | survives as a finite lower-bound certificate | it detects witness mass exactly; canonical code transversals make it exponentially small under diffuse queries |
| Mixed-relation holonomy | survives as an exact algebraic gluing state | composition creates `Hom(Z/Z_loc,W)`, but its gauge dimension alone need not be response-visible |
| Presented carrier | survives for distance-transform/min-plus landscapes | response entropy is carrier Hausdorff entropy up to presentation radius; the query law controls witness exposure |
| Metric-quotient synchronization | survives conditionally | small fibres and lift defect give a uniform quotient decoder, but strict compression also needs lower projected entropy and descended composition |
| Raw carrier/relation dimension | rejected as a response proxy | surjective maps, redundant alphabets, and two-scale metrics collapse exponentially many gauges |
| Sparse-flat spectrum of a quotient leader ball | survives exactly for one directed Grassmannian ball; rejected as a two-sided response state | it gives an exact ball identity, but isometric quotient norms can have a linear gap in rooted Hausdorff response |
| Puncturing/anticode quotient as entropy-optimal | rejected even for binary line carriers | line packing is unrestricted binary coding, and an explicit probabilistic line cover has strictly smaller rate than puncturing |
| Rooted metric extension data | forced, not yet compressed | quotient geometry and kernel geometry must be coupled to answer arbitrary future carrier queries; no strict sufficient quotient is yet proved |
| Directed response alphabet `r(a,b)` | survives exactly for product-response distances | its two oriented sums give the exact uniform metric after arbitrary product depth; it forgets the full local functions but not this declared query |
| Low-word count of `C+C'` plus injection distance | rejected as an improved asymptotic lower route | it yields a valid all-rank packing, but its exponent is always bounded by the common-host Gilbert exponent |
| Deterministic Parisi-like overlap state | narrowed | a finite synchronization theorem exists, but no natural deterministic hypothesis is yet known to force its cross-root linkage |

## 9. Ten questions every candidate must answer

1. Which examples force the definition?
2. What is the declared normalization and query experiment?
3. Which previous failures become immediate in its language?
4. What does the state forget?
5. What does it retain exactly or approximately?
6. What is its packing/covering or description complexity?
7. What is its algebra under the intended composition?
8. In what topology is it compact, and is the queried extreme continuous?
9. Which finite constrained objects realize its limits?
10. What theorem about it is not the tautology “if the state converges, the
    maximum converges”?

The present roof answers these questions for fixed linear interfaces and
unrestricted finite landscapes.  Its bi-affine closure theorem also gives an
exact polynomial-state dynamic program for fixed-rank Curie--Weiss/Potts-type
mean-field ground states.  It does not yet answer them for growing bridge
interfaces in dense sign matrices.

## 10. Three distinct complexity coordinates

The second investigation shows that “state size” is too coarse.  Three
coordinates must be reported separately.

1. **Exact algebraic size:** the number of labels needed for exact closure.
   A boundary kernel has `Q^2` independent entries; a syndrome fragment over
   `F_2^w` has `2^w-1` possible support bits.
2. **Uniform response resolution:** the factor or covering complexity at a
   declared sup-norm error.  Tropical crossings can protect the exact state
   throughout a fixed lattice-scale neighborhood.
3. **Query-weighted resolution:** the information or factor complexity under
   a declared distribution of interventions.  Posterior width controls this
   for Hilbert response embeddings; weighted tropical exposure controls one
   class of factorized tables.  Uniform hardness need not survive here.

These coordinates are inequivalent.  The graph-code distance table has
exponential exact tropical rank and uniform robustness below one half unit,
yet after normalization it has a rank-one approximation with vanishing
uniform mean-square error.  A claim of extremal complexity must therefore
state its query law and distortion scale.

The syndrome block theorem shows how a lattice-scale coordinate effect can
nevertheless become macroscopic without separately paying channels.  One
future fragment selects many direct-sum blocks, and their covering-radius
contributions add before the response error is charged.  This **joint
exposure** is an operation in the feature algebra, not merely a different norm
on the old one-bit queries.

## 11. Current selected abstraction

A **query-generated feature algebra** begins with the observable exposed to
the environment, closes under the declared compositions and contractions,
and identifies landscapes with the same resulting response.  Its value is
not the definition but two possible theorems:

- an algebraic quotient theorem proving that the closure is strictly smaller
  than the landscape, as for syndrome supports; or
- an information theorem proving that a stated distortion requires a stated
  rate, as for posterior response width.

Synchronization and fixed-effective-dimension convexification are two ways
the closure can collapse.  Tropical crossings and response packing are ways
to prove that it cannot.  No single one of these objects is presently a
universal deterministic analogue of the Parisi order parameter.

## 12. Presented carriers: the first two-sided growth law

The third investigation isolates a class in which both information growth
and compression follow from one state.  A **presented carrier** in a query
metric `(X,d)` is a nonempty set `C` with an access cost
`pi:C->[0,p]`; its response is

```math
F_{C,\pi}(x)=\min_{c\in C}\{d(x,c)+\pi(c)\}.
```

This object is forced by mixed relations.  Individually valid gauges become
incompatible on a new relation space; their holonomy image is `C`, while the
number of relation letters needed to reach a point is `pi`.

Four quantities must be kept distinct.

1. **Latent relation rank** counts possible compatibility coordinates.  It is
   only a supply of potential information.
2. **Carrier entropy** is the Hausdorff packing/covering entropy of the image
   sets in the declared endpoint metric.
3. **Presentation radius** `p` is the price of accessing the carrier.  It is
   subtracted once from uniform separations.
4. **Exposure mass** is the query-law mass near separating witnesses.  It is
   necessary for diffuse, rather than uniform, information bounds.

Theorem 12.3 proves that items 2--4 determine response complexity inside this
class.  Raw relation dimension does not: it can map surjectively to one
carrier, or live inside metric fibres invisible at the chosen scale.

The complementary compression primitive is a metric quotient
`varpi:X->Y`.  It is operationally useful only if:

1. fibres and lift defect are subscale;
2. projected carriers have a proved smaller covering rate;
3. the declared carrier composition descends to `Y`; and
4. the presentation-radius certificate stays controlled.

Under those checks, Theorem 12.5 gives a strict approximate feature algebra.
This is the first current framework in which a single object predicts both
`Theta(D kappa)` growth and exponential collapse in different models.  Its
scope remains distance-transform/min-plus response families; treating every
extremal landscape as a “carrier” without proving such a representation would
be vocabulary rather than theory.

For linear carriers, Theorem 13.1 replaces the qualitative four-part test by
two scale-dependent ranks.  The separated rank `s_W(Delta)` is a lower
capacity certificate; a synchronizing quotient rank is an upper compression
certificate.  Fibre geometry forces the generalized Singleton constraint

```math
s_W(a)\le\dim Y.
```

This relation is structural rather than definitional: it becomes the
classical Singleton bound under Hamming puncturing and the rank-metric
Singleton bound under row projection.  Gabidulin hosts then turn it back into
an extremal-response information theorem.  The optimal quotient rank is
exactly `dim W-A_W(a)`, where `A_W(a)` is the largest linear anticode
dimension.  Thus the unresolved quantity is the code--anticode gap

```math
\gamma_W(a)=\dim W-A_W(a)-s_W(a).
```

It vanishes in the two-scale and rank-metric examples but is linearly positive
in binary Hamming space by sphere packing.  The leading question is whether
actual Grassmannian carrier entropy fills that gap or requires a third
invariant between code packing and anticode quotienting.

The Hamming Grassmannian audit resolves part of that question.  The
one-channel entropy is, up to one codeword, the unrestricted binary coding
number; the puncturing quotient is therefore exponentially nonminimal.
Directed balls admit an exact sparse-flat formula: count linear flats inside
the quotient coset-leader ball.  But this is not a new sufficient state.
Two carriers can have isometric quotient normed spaces, hence identical
sparse-flat spectra at every scale, while a rooted Hausdorff query differs
linearly after direct sums.

The surviving abstraction is consequently a **rooted metric extension**:
the quotient leader geometry, the metric on the carrier/kernel, and the way
quotient representatives couple to it.  This is a requirement extracted by
a scalable counterexample, not yet a claimed compression theorem.  At
`k=Theta(D)`, systematic coordinates turn the coarse problem into a
`2^k`-ary column code, but Reed--Solomon codes already saturate its Singleton
bound.  A useful next theorem must control the coherent same-input
recoupling that this column shadow forgets; adding more ordinary coding
statistics cannot do so.

One strict quotient nevertheless survives.  For a finite alphabet of local
responses,

```math
r(a,b)=\sup_x(f_a(x)-f_b(x))
```

is an exact feature algebra for pairwise uniform distance under direct
products: product distances are the maximum of the two oriented sums of
`r`.  A two-sided carrier gap `d` and presentation radius `p` imply local
weights at least `d-p`; an outer code amplifies them without charging
matching channels.  The Hamming simplex and rank-metric multiplication
families validate the same theorem.  This is a generative composition law,
not a universal response state: it controls the declared product-distance
query and says nothing about arbitrary nonproduct couplings.

## 13. Benchmark correction: three tests, not one

The separator benchmarks force a sharper test for any proposed extremal
state.  Three logically independent facts are required.

1. **Context factorization.**  Every allowed future must meet the stored
   object through a declared interface.  Conditional maxima on that interface
   are then sufficient.
2. **Query exposure.**  Allowed futures must actually distinguish the claimed
   coordinates.  In pure Max-Cut, positive-edge gadgets expose projective
   boundary assignments, not oriented assignments; global flip is a genuine
   contextual gauge.
3. **Language realizability.**  The model must generate a large separated
   family of interface profiles.  Exposure alone does not give an information
   lower bound.  The gauge-anchored lookup theorem supplies this missing step
   for pure Max-Cut.

Approximation adds a fourth datum: regularity of the realizable profiles in
an interface metric.  The same exact projective table has full cube entropy
without a promise, but a one-Lipschitz promise has Hamming-covering-code
entropy.  Therefore exact contextual minimality does not determine lossy
state growth.

The benchmark states fit one profile calculus but have different growth
laws:

| model | contextual interface | exact state | growth mechanism |
|---|---|---|---|
| width-`w` Max-Cut | projective boundary assignment | conditional cut profile | full profile cube without normalization; full Lipschitz ball at unit boundary load |
| width-`w` pairwise Ising | oriented boundary assignment | conditional energy profile | `2^w` independently realizable coordinates |
| fixed-rank mean field | aggregate feature in `Z^d` | roof/profile over reachable sums | only `O(n^d)` interface types |
| fixed weighted automaton | reachable suffix experiment | residual response function | quotient of the presented forward vector |

This is explanatory rather than a new name for dynamic programming.  Its
generative content is in model-specific realizability and distortion
theorems, such as Theorems 16.1--16.2.

The directed table now also survives one nonproduct operation.  Under a
common min-plus continuation, an oriented response can be lost only after all
of its exact directed witnesses cease to be exposed.  More strongly, distance
cones on any finite metric close exactly: isometry holonomies multiply and
response strengths combine by a bottleneck minimum.  Uniform entrywise
perturbations accumulate additively.  This is the concrete gauge/cocycle
mechanism anticipated by the earlier product theorem.

The Max-Cut normalization benchmark adds a resource warning.  A private
universal compiler followed by a single resource-bounded distance shell
realizes every Lipschitz response at the shell's boundary cost.  Hence a
local sensitivity bound need not imply low response complexity: under unit
load the response class still has exponentially many bits at macroscopic
distortion.  The next compression theorem must charge a global resource
(component size, description length, internal precision, or grammar), not
only exposed-interface load.

These facts suggest one reusable mechanism rather than two unrelated
examples: a metric distance kernel is an idempotent tropical projector on
Lipschitz profiles, while isometries act as a holonomy group and unequal
strengths compose by a bottleneck.  This still does not imply closure for
arbitrary transfer kernels; the next extension must recognize or approximate
this metric algebra without storing all kernel rows.

## 14. Benchmark revision: three distinct complexity resources

The validation campaign separates three quantities that should no longer be
called simply “state size.”

1. **Interface sensitivity** controls regularity.  In Max-Cut, boundary load
   is exactly the weighted-Hamming Lipschitz seminorm of the response.
2. **Presentation complexity** controls how many optimizer cells a shared
   parameter grammar can generate.  An `m`-parameter binary max-affine
   presentation has only `exp(O(m^2))` normal-fan cells even with
   arbitrary-precision weights.
3. **Derivative-compatible quotient complexity** controls indefinite reuse.
   A small cover at one time is not an automaton; the quotient must also
   descend under every future derivative.  Tropical lumpability and metric
   semilattices are positive mechanisms, while the transition-toll example
   shows how a vanishing local defect can accumulate.

These resources are independent.  Unit interface load admits the full
Lipschitz ball if an exponential private compiler is allowed.  A
polynomial-size compiler has polynomial response-bit entropy, but that alone
does not supply an exact repeated transition.  Conversely, a tiny holonomy or
lumped transition algebra can transport a profile whose own response
description remains large.

The resulting candidate law is operational:

> Future-response compression requires both a small realizable response
> image and a congruence for the future semigroup.  Regularity bounds the
> image metrically; shared-parameter presentation bounds its entropy;
> synchronization/lumpability makes the same quotient reusable.

This is stronger than relabeling dynamic programming because each clause has
a separate theorem and a separate counterexample.  It is not yet a universal
rate formula: the gap between query-net upper bounds and robust-exposure
lower bounds remains model-dependent.

## 15. Depth-uniform reuse must declare its perturbation quantifier

The benchmark campaign first suggested a gauge--reset dichotomy and then
falsified its narrow coherent form.  Three finite certificates are now
proved, but they apply under different quantifiers.

1. **Cohomological cancellation.**  If every transition residual is an
   incoming interface potential minus an outgoing one, the internal
   potentials telescope before min-plus optimization.  Rectangle defects
   plus adjacent-interface circulations give an exact finite recognition
   test on complete interfaces.  On a repeatable context graph, every
   projective cycle holonomy must vanish exactly; any nonzero holonomy grows
   linearly under repetition.
2. **Projective reset.**  A continuation whose entire quotient image has
   diameter `rho` erases all error from before that continuation.  If such a
   reset of length at most `L` recurs with tails of length at most `L`, a
   one-step semiconjugacy defect `epsilon` costs at most
   `rho+2L epsilon`, independently of total depth.
3. **Algebraic absorption.**  If two coherent nonexpansive continuation
   families factor through the same finite semigroup, every long word has a
   bounded normal form.  Generator error `epsilon` then costs at most
   `L epsilon`, where `L` is the largest shortest representative length.

The third item is a finite realization of stable zero-increment recurrence,
not a universal category beyond every possible notion of zero holonomy.
Nearby idempotent max-plus clamps prove only that entrywise kernel gauges plus
small **full-image** resets are incomplete.  Their paired orbits become
stationary after one transient.

Against fresh arbitrary residuals, there is a sharper converse.  On a
coordinate-selector cell, depth-uniform `O(epsilon)` stability holds exactly,
up to explicit `r(r-1)` constants, when rank-one **tangent** selector products
occur with bounded gaps after endpoint gauges are removed.  A tangent reset
may kill every transported error direction even when the full nonlinear
state image remains large.  On a recurrent affine-selector cell, the
translation must be a twisted coboundary; equivalently, all functional-graph
cycles have the same translation mean.

These are not two descriptions of ordinary contraction.  On the full
finite-dimensional projective domain, every all-finite max-plus linear map
has global Hilbert coefficient either zero (additive rank one, hence a full
reset) or one.  A weak Ising bond has coefficient one but small image
diameter, so it is useful as a reset rather than as a multiplicative
contraction.

The transition-toll and translated-clamp counterexamples fail these coherent
relations and drift, despite small one-step error.  Thus a proposed stability
theorem must say whether defects are fixed by a common algebra or may vary
adversarially at each use.  A paired-selector theorem across switching cells,
not a bare static cover or global contraction coefficient, is the next
canonical target.

## 16. Static response width has an upstream two-sided certificate

For a metric query interface, response entropy can be bounded without first
covering the response class itself.

- Values on a query net, with scalar quantization, give a landmark upper
  bound.
- A balanced collection of queries that realizes every half-positive,
  half-negative margin pattern gives an external-cover lower bound.
- For the full Lipschitz response ball, this balanced exposure dimension is
  exactly the size of a query packing at twice the margin, rounded down to an
  even number.

The certificate separates two costs that a one-number fat dimension would
conflate: the number of independently exposed queries and the scalar
precision stored at each one.  A one-dimensional response interval has
constant fat dimension but logarithmic precision entropy.  In weighted
automata, robust suffix pins expose an isometric coordinate cube, while the
subgraph class is controlled by upper-orthant VC dimension.  In normalized
Max-Cut, the full Lipschitz compiler turns query packing directly into
response exposure.

Presentation complexity supplies a second limitation.  Binary max-affine
grammars with `m` shared parameters have an `O(m^2)` radius-bounded entropy
ceiling, but a single high-facet `0/1` polytope already forces
`Omega(m log m)` robust response bits.  The remaining gap concerns common
normal fans of many support polytopes.  This is a static law only: it must be
paired with the congruence mechanisms of Section 15 before the same summary
can be reused at arbitrary depth.

## 17. Feature growth is controlled by exposure and congruence, not interface cardinality alone

The mean-field benchmark supplies the missing contrast to separator tables.
A finite interface with `q` independently realizable lookup coordinates can
force `q` response degrees of freedom.  A much larger microscopic system can
instead have a one-parameter query that exposes a structured profile whose
slopes form a commutative multiset.  Its exact merge is histogram addition,
so a fixed-grid approximation has polynomially many states and charges each
microscopic field once rather than once per merge depth.

The quadratic benchmark sharpens this further.  The visible state is not the
raw conditional table but its concave envelope, because linear futures only
probe supporting faces.  Separate affinity of the cross interaction makes
that observational quotient a congruence.  Sufficient curvature then forces
all profiles onto one common exposed chord and synchronizes the state to a
single aggregate.

This suggests a benchmark-tested two-stage law.

1. **Exposure geometry** determines which quotient of a conditional profile
   is semantically observable (all coordinates, a Lipschitz function, a
   concave roof, or a residual class).
2. **Derivative congruence** determines whether that quotient remains valid
   after composition (lookup gluing, histogram addition, bilinear roof
   convolution, tropical lumpability, gauge cancellation, or reset).

Neither stage determines state growth by itself.  Boundary cardinality can
be exponential while a shared presentation is small; aggregate cardinality
can be polynomial while arbitrary real slopes make exact state cardinality
uncountable.  Exact table width, number of equivalence classes, lossy metric
entropy, and update complexity must therefore remain separate resources.

This is now more than vocabulary: the separator compiler, automaton
lumpability theorem, and mean-field roof theorem give different realizable
response images and different congruences, while the transition-toll and
changing-curvature examples independently falsify either clause in
isolation.  What remains missing is a general necessary structural theorem
for depth-uniform approximate congruences.

## 18. Finite additive feature algebras have an arithmetic growth exponent

When every component response is a nonnegative integer combination of a
fixed finite set of atom responses and composition adds those counts, the
exact contextual quotient is itself a finitely generated abelian monoid.
Its mass-`n` state count grows as

```math
Theta(n^(r_Z)),
```

where `r_Z` is the integer rank of the atom-response differences after any
projective baseline is removed.  This is the first benchmark law that
predicts a state-growth exponent directly from the feature algebra.

Arithmetic rank is not a robust dimension.  Irrational responses on one
real query can have high integer rank while integer combinations approach
each other arbitrarily closely.  Approximate response complexity instead
depends on the smallest real singular value of the query map and on its
lattice margin.  These three quantities answer different questions:

1. `r_Z` controls exact polynomial growth under unlimited precision;
2. real conditioning controls macroscopic packing and covering; and
3. lattice margin controls the resolution below which every exact state is
   still distinguishable.

The separator, Ising, mean-field, and automaton benchmarks now fit one
operational scheme without sharing one numerical dimension.  Arbitrary
separator lookup futures expose exponentially many independent semantic
coordinates.  Fixed-width transfer and lumped-automaton states have bounded
carrier rank.  Finite-atom mean field has a growing histogram simplex whose
degree is `r_Z`.  In every case the state is useful only because the same
quotient is also a congruence for future composition.

The same law has a continuous lossy form.  Quantize each microscopic atom in
its **response metric** before any composition and store the empirical type
histogram.  An atom-net of size `D(eta)` gives response error `n eta` and a
binomial number of histogram states.  Because the quantizer is applied once
per atom, the update remains exact and the error is extensive in mass rather
than in recursion depth.  This separates genuine atomic compression from
repeatedly rounding an already aggregated state, which can drift.

The sharp stars-and-bars entropy matters: an atom net with `D_n=o(n)` types
already has `o(n)` histogram bits, even when `D_n log n` is not `o(n)`.

## 19. Switching transports are priced by suffix-row memory

A max-plus secant across optimizer switches is row-stochastic. Its exact
gain from fresh projective residuals is the cumulative total-variation
separation of terminal suffix rows. This gives two distinct forgetting
mechanisms:

1. a finite consensus product kills old response directions exactly; and
2. repeated fractional mixing makes them geometrically summable without an
   exact reset.

A face label or tie flag is therefore not a sufficient stability state.
What matters is transported row distinguishability. Quantitative scrambling
gives deterministic synchronization, while an unweighted active-face graph
can miss arbitrarily weak mixing. The exact lower formula is adversarial;
fixed coherent kernels may have smaller gain because residuals and switches
are dynamically coupled.

## 20. Cycle algebra is complete only after symbolic realization

For an exact finite affine-selector presentation, one ordered coordinate
witness turns coherent response error into an additive graph cocycle. Its
positive cycles are precisely directed drift, and zero cycle holonomy is
precisely two-sided boundedness. Paired channels first require the joint
cross-difference carrier; diagonal error is not closed when selectors differ.

The hard information resource is therefore not the cycle test but the
control quotient on which paths are genuine. Nonempty adjacency between
local faces can create spurious repeatable words. A useful finite switching
state must be both a response quotient and path-realizing. For robust reset
queries, this principle itself compresses the old suffix-set state to the
kernel partition of the total selector product, which is the coarsest
partition component at a fixed control vertex for the full selector alphabet.

## 21. Dynamic compression has three noninterchangeable certificates

The switching benchmark forces a refinement of the static two-stage law.
There are three ways a small future state can remain valid.

1. **Exact contextual congruence.** Pull back the observation partition by
   all future words. Finite stabilization is equivalent to a finite exact
   path-realizing quotient. A finite oriented affine arrangement closed under
   branch pullback is a checkable sufficient certificate.
2. **Metric forgetting.** An approximate quotient or fresh quantization may
   be reused when every old residual is attenuated with summable suffix gain.
   Block contraction changes the resolution at which the quotient's metric
   entropy must be paid.
3. **Cocycle cancellation.** A projective control can return exactly while
   its scalar baseline accumulates error. Uniform absolute response requires
   every repeatable reward-discrepancy cycle to vanish, or a quantitative
   contraction which absorbs it.

These mechanisms solve different obligations. A static response cover gives
none of them. The compact Cantor encoder--decoder has only five affine pieces
and a bounded one-step response cover, yet its horizon-`T` predictive memory
is exponential. Conversely, an exact Ising transfer profile has a large
exposed response image but no future amplification of a one-time sup error.

The finite pullback certificate must use full **oriented affine forms**.
Unoriented normal fans forget side exchange, affine offsets, and scaling.
Tie-value paths and tangent-selector paths must also remain separate: a fixed
tie rule is exact for values, while independent tangents may have no common
perturbation cone.

The quantitative static--dynamic interaction is now a scale transform, not a
formal product. If `G_T` is the sum of suffix Lipschitz gains, a finite
sequential response simulator pays the internal covering number at scale
`epsilon/G_T`. For the entire branching context tree, local response entropy
at diameter `D rho^k` is paid `q^k` times. Both laws are sharp, but neither
alone supplies an exact semiconjugacy.

## 22. Compact arithmetic closes symbolic switching, while mixing prices lossy reuse

Unit-selector dynamics have an additional exact closure mechanism. On a
common compact projective carrier, selector pullbacks move every labelled
normal through a finite set. If affine offsets lie in a discrete additive
group, only finitely many translated walls can meet the carrier. Rational
polyhedral selector systems therefore admit an effectively constructible
finite path-realizing refinement for every finite polyhedral observation
coloring. Compactness alone is insufficient: an irrational circle rotation
has dense wall pullbacks and no finite autonomous predictor below half the
response diameter.

For approximate reuse, the intrinsic metric on a finite predictor is not its
one-step output distance but the supremum of output distance over all common
suffixes. In this behavioral metric every transition is nonexpansive, and
the encoder misses physical one-step recoupling by at most twice the original
response error. A strict behavioral margin turns this into exact
semiconjugacy.

Quantization then obeys a precise static--dynamic law:

```math
response error
<= original predictor error
 + (behavioral net radius)(suffix gain).
```

A weighted legal-control graph computes the suffix gain. It is bounded
exactly when every reachable coefficient-one cycle is absent, with the sharp
constant obtained from max-times path products or a maximum-cycle-mean
potential. This makes contraction a quantitative information resource: a
finite response cover can have two states while every reusable predictor
needs linearly many states if its forgetting time is linear.

Symbolic closure, metric forgetting, and scalar reward compatibility remain
separate. Rational compact refinement controls which branches can occur; it
does not by itself prevent an accumulated reward cocycle from drifting.

## 23. Recurrent germs, not raw periodic points, control reward memory

Compact rational selectors have a finite semigroup of realized projective
affine germs. Lifting the exact symbolic quotient by these germs turns each
affine transition reward into an affine function on one fixed seed atom.
Modulo scalar tolls, cumulative reward has a complete finite dichotomy:

- zero label on every lifted directed cycle leaves only a simple-path
  residual and hence a uniform depth-independent error;
- a nonconstant lifted cycle label is a genuine repeatable response and
  pumps linearly.

Exact state-potential cohomology is stronger than bounded response. It also
requires agreement across transient coterminal paths. Noninvertible selector
diamonds show that all ordinary periodic sums may vanish while this exact
potential fails, even though cumulative error remains bounded.

The recurrent label is itself an exposed query. If one return word has a
continuum of cycle means, bounded absolute simulation needs infinitely many
states; at positive per-step distortion, the packing number of the cycle-mean
image lower-bounds dynamic memory. This is the first direct law converting
holonomy geometry into response rate--distortion.

For a finite path-realizing germ graph this extends to a two-sided law. Give
two hidden seeds the maximum, over reachable directed cycles, of their
average cumulative-reward separation. Packing this cycle-response
pseudometric lower-bounds every deterministic simulator, while an internal
cover gives a simulator within the same radius, up to the separate finite
control factor and one bounded transient toll. Thus the reusable asymptotic
state is exactly the recurrent response image; acyclic history affects only
the nonextensive constant.

## 24. Exact invariant nets are a third alternative to congruence or forgetting

Fresh projection onto a static net generally drifts, but a net preserved
exactly by every transition introduces no fresh residual. For deterministic
continuous rational unit-selector maps on compact convex rational projective
polytopes, with a fixed whole-fibre target for each control/input pair,
sufficiently fine compatible ambient rational grids are finite invariant
nets. Selector nonexpansiveness then keeps the initial shadow error unchanged
under every switching word, even when no strict contraction occurs.
State-dependent target controls must first be carried by an exact symbolic
refinement.

This arithmetic phase-locking differs from exact contextual congruence: two
nearby raw points need not re-encode to the same grid state after a step. It
also differs from metric forgetting: old error need not shrink. The simulator
evolves one genuine nearby raw trajectory forever. Irrational rotation marks
the boundary—without a discrete invariant grid, equicontinuity alone supplies
no finite autonomous predictor.

## 25. Finiteness and compression are quantitatively different

Compact rational selector systems admit finite exact refinements and
invariant grids, but their smallest exposed recurrent carrier can still have
`2^(r-o(r))` states. Two permutation inputs and one repeatable scalar probe
already expose the constant-weight orbit of `S_r`. The correct quotient
removes a factorial stabilizer, yet remains exponentially large. Structural
closure must therefore be accompanied by response metric entropy; a finite
presentation alone is not an information bound.

A useful sharp boundary is generator count. One globally active selector on
`r` coordinates has only `exp(O(sqrt(r log r)))` recurrent iterate germs on a
compact carrier. Two switched selectors can have exponentially many exposed
classes. Interaction creates the information growth.

## 26. Approximate congruences are selected by cycle incidence, not pair distance

For a proposed finite input congruence, optimizing scalar tolls is exactly a
minimax cycle-mean linear program on the raw transition graph. Raw-cycle
defects vanish iff total error stays bounded. This is the finite quantitative
compatibility law missing from a static response cover.

Feasible bounded-error partitions are not closed under join and may have no
unique coarsest member. Even when every pair of raw states has zero
same-input asymptotic response distance, merging them can create a quotient
cycle whose required letter tolls are inconsistent. Thus a metric quotient
and a transition congruence cannot be optimized independently: the
composition law changes which response combinations are recurrent.

## 27. Recurrent exposure is word-synchronized

A critical graph compresses one fixed max-plus product. It does not in
general preserve responses to each aligned switching word. When the
projective product semigroup is finite, its weighted Cayley graph is the
path-realizing replacement: synchronizing two such graphs turns the exact
supremum over wordwise spectral-response gaps into a finite cycle-mean
problem. The carrier is intrinsic and checkable, but output-sensitive; its
size may itself be exponential.

## 28. A small circuit can hide an exponentially large response orbit

Polynomial presentation size does not bound dynamic response information.
One `O(r)`-gate continuous rational lattice-PWA selector map can implement a
binary counter of period `2^(r/2)`. A single bounded identity probe exposes
all phases at a fixed gap, forcing the full exponential state count below
that gap.

This also sharpens the query principle. With only the counter evolution
letter, every phase sees a rotation of one reward cycle and one mean toll has
bounded error. Adding the ability to freeze and repeat one phase observation
turns latent orbit complexity into extensive response separation. State
complexity is a property of dynamics *and its future interface*, not of the
dynamical map alone.

## 29. Optimal response congruence is a clustering problem, not refinement

Once a candidate partition is fixed, its best asymptotic toll error is a
polynomial cycle-mean LP. Selecting the smallest feasible partition is
different: it is NP-complete even with identity transitions and fixed error
`1/2`, by an exact reduction from graph coloring. Rank-one reset dynamics
retain the hardness.

Thus contraction controls how approximation errors propagate but does not
make the response image easy to cluster. Positive structural theorems need
low-dimensional reward geometry, an algebraic path lift, or another special
property; generic monotone partition refinement cannot find the optimum.

## 30. Path realization is stronger than equality of optimized responses

A coherent path-lifting relation is a powerful sufficient certificate, but
it is not intrinsic to scalar word responses. Exact equality of every
aligned-word spectral radius can coexist with a linear-size lower bound on
every low-defect path lift. The optimum may choose a different microscopic
critical witness after seeing each full word.

This is a quantifier obstruction:

```math
for every future word, there exists an optimizer
```

does not imply

```math
there exists one finite relation continuing all represented states.
```

A converse needs deterministic synchronization, exposed uniqueness stable
under prefixes, or another hypothesis that permits this interchange.

For `{0,-C}` max-plus kernels the distinction is exact. Equality with the
scalar zero response says that every word-composed zero relation contains a
directed cycle (**periodic completeness**). A scalar coherent path lift says
that every one-letter zero relation is left-total (**local continuation**).
Failure of periodic completeness supplies one finite word whose repetition
pumps a linear response gap, but periodic completeness does not imply local
continuation: the fixed-binary de Bruijn family separates them exponentially.
This is an intrinsic response-versus-simulation boundary, not another choice
of quotient terminology.

Even maximal metric forgetting is not enough. In the same de Bruijn family,
every sufficiently long word product has projective rank one and one critical
node, but the identity of that node is the word's initial window. Contraction
forgets an incoming message *after the word is fixed*; reusable compression
must choose a state *before the future word is known*. The missing property
is therefore cross-word coherence of exposed images, not merely scrambling,
primitivity, or unique optimization within each product.

## 31. Query-matched path realization can be anticipatory

The correct finite path certificate depends on the query quantifiers. For a
scalar unrooted tropical spectrum, one need not continue from every raw
representative. It is enough to carry a nonempty subset of possible raw
witnesses whose next subset lies inside the current good-edge image. Once a
finite future is known, a terminal witness can be lifted backwards to a
genuine path. Repeating a word then closes a raw cycle by finiteness.

For exact scalar tolls this **survival carrier** is complete: it exists iff
the good-edge subset automaton never dies. Its canonical state has at most
`2^r-1` endpoint subsets, not `2^(r^2)` full relations. Failure produces a
finite mortal word and hence a pumpable response gap. This is a genuine
finite tropical lumpability law, but it is not ordered with rowwise
bisimulation in general. In the free-tail de Bruijn example its one-state
subset carrier answers the unrooted query while rooted continuation costs
`2^m` states.

## 32. Weighted survival is a finite support cocycle

For a nonflat coarse response, mere survival can follow a subcritical loop.
The correct finite extension attaches selected endpoint supports to coarse
states and requires backward-surjective near-optimal lifts of every declared
coarse edge. Local shortfalls need not vanish: a potential on support states
absorbs their transient part, and only the residual cycle mean accumulates.

This yields a quantitative all-word theorem and separates symbolic coarse
paths, endpoint realization, response error, and support-state count. It can
be strictly smaller than rowwise simulation because the predecessor within a
support can depend on the future coarse edge: width-two Ising has `1` scalar
output, `2` anticipatory weighted states, and `4` forward transfer states.

The reverse separation is exponential. A deterministic de Bruijn shift has
one scalar response state and one rowwise state, yet every exact
backward-surjective support carrier has `q^m` states; its optimal certificate
toll at size `N` is `Theta(C/(1+log_qN))`. Anticipatory and forward carriers
are therefore incomparable proof architectures. Neither state count is
intrinsic response information unless the declared future queries expose its
hidden witness phase.

## 33. A bridge is measured by its exposed response transform

For an interaction `B(x,y)`, the canonical past state under arbitrary future
weights is

```math
P_Bh(y)=max_x\{h(x)+B(x,y)\}.
```

The future-context metric is exactly the sup distance between these tables.
This one transform unifies separator tables, low-rank upper roofs, and
symmetry-orbit responses. Rank, sparsity, and graph treewidth are only useful
when they actually shrink the realizable image of this transform.

The distinction is quantitative. Fixed rank gives an exact
finite-dimensional roof algebra and fixed-error feature nets independent of
the microscopic state count, although rank one can still have `2^n` exact
exposed atoms. Growing rank forces `2^(Omega(r))` response bits in an
unrestricted bounded class. Conversely, a degree-one matching bridge can
expose an exponential code and require exponentially many macroscopic bits.
The relevant sparse parameter is live interface geometry, not local degree.

## 34. Deterministic synchronization is a common optimizer section

A quotient closes under interaction when every tuple of quotient labels has
one microscopic representative family simultaneously realizing all discarded
pair features. For permutation-invariant spin blocks with bridges
`alpha I+beta J`, nonnegative identity coefficients have such a section:
nested plus sets maximize every pairwise overlap at once.

For signed identity channels the exact criterion is cycle balance. A vertex
gauge makes all signs positive iff every cycle sign product is positive.
An unbalanced unit cycle loses exactly `2n` against the sum of separately
optimized pair responses on balanced blocks. Thus gauge holonomy is an
observable obstruction to the separable pair-potential algebra.

This is not a universal lower bound on joint tables; one can always retain
more compatibility information. Its value is generative: it gives a
checkable full-rank dense class with an exact magnetization factor algebra and
checkable full-rank dense class with an exact joint magnetization reduction
and a finite-dimensional thermodynamic limit, and it identifies the precise
cycle witness when that reduction fails.  Counts are not automatically a
serial state when an old microscopic alignment is frozen and exposed to a
new future; the common-section theorem assumes joint reoptimization in the
declared synchronized family.

## 35. Extremal replacement is scale-sensitive

Cut-norm proximity controls every finite-label pair-energy maximum uniformly
over arbitrary conditional futures, including a future that pins one rare
labeling. This makes weak regularity a genuine extremal replacement theorem
at the dense `n^2` scale, not merely a statement about bulk statistics.

The same statement explains its own failure near the motivating scale.
Preservation at scale `L_n` requires cut error `o(L_n)`. Generic weak
regularity at `L_n=n^(3/2)` needs accuracy `o(n^(-1/2))` and may introduce
exponentially many blocks. A useful renormalization theory must therefore
carry the normalization in its definition; a compact state at one leading
scale need not retain any information at a finer extremal scale.

## 36. Proof memory becomes semantic only through observability

A backward-surjective support can remember exponentially many phases while
the optimized response and a forward path lift each have one state.  Carrier
size is therefore not an information lower bound by itself.

The missing condition is two-sided and query-relative.  A backward support
must sit inside a forward endpoint envelope, and the declared probes must
take the same maxima on the two sets.  Under a uniform energy gap, bounded
probes then filter exactly to the endpoint max-profile, turning separated
phase profiles into a genuine response packing.  A cyclic navigation group
and one local symbol readout already expose every deterministic de Bruijn
phase; arbitrary lookup tables are unnecessary.

Finite leakage imposes a hard ceiling: max-plus continuation is
nonexpansive, so endpoint vectors within `C` remain within `C` under every
common future.  Persistent rate separation requires hard reachability, a
phase-preserving reset/filter, or a bounded probe horizon.  Observability,
not proof architecture, decides whether hidden witness phase is information.

## 37. Bridge compression has several incomparable structural resources

The canonical object is the realizable image of

```math
P_Bh(y)=max_x\{h(x)+B(x,y)\}.
```

Low algebraic rank compresses it through an upper roof.  A small live vertex
cover gives a `2^k` conditioned table and this dependence is worst-case
sharp.  A common optimizer section can compress a dense full-rank family by
orbit counts.  None of these follows from bounded degree, density, or
spectral norm alone.

At the negative endpoint, a dense sign bridge with operator norm
`O(sqrt n)` can expose exponentially many independent response coordinates
with `n^(3/2)` margins for unrestricted children.  Even weighted linear
children retain an `Omega(n)`-bit projective packing. Complete sign
quadratics themselves now retain `Omega(n)` bits: planted gauge-ferromagnetic
poles turn the random bridge correlation matrix into an exponential
projective packing. Thus quadratic syntax alone is not a compression
resource.

The qualification is structural, not cosmetic. Those pole landscapes have
`Theta(n^2)` cap and spread. A theory near the signing problem must use
rigidity of the bounded-cap or near-minimizing subclass. Coefficient-Hamming
covering gives a strict universal saving. More strongly, universal random
thinning discards a constant `Theta(epsilon^2)` fraction into a sparse
weighted surrogate while retaining all Boolean energies at target error.
Both remain `Theta(n^2)` bits; the present information bracket is only
`Omega(n)` to `O(n^2)`.

## 38. Rare-event states may be created only after renormalization

Finite-port response quotients are not the only compression mechanism. In a
boundary-case branching random walk, the exact finite-depth derivative
state is `(W,Z)`, but critical centering makes `W` vanish and the limiting
unmarked extremal process becomes a decorated Cox process controlled by one
scalar derivative mass `Z`. The scalar composes by a smoothing transform.

This suggests a genuinely orthogonal design principle: a macroscopic state
can emerge only after a model-specific rare-event renormalization, even when
no finite-depth exact quotient exists. It does not weaken query relativity.
Genealogy-marked futures distinguish branchwise derivative masses with the
same total, so the correct marked state is a measure rather than a scalar.

## 39. Finite query mass permits universal feature thinning

For a signed linear combination of a fixed bounded feature dictionary, the
relevant concentration parameter is not the raw number of coefficient
vectors. It is

```math
V_Phi=max_x sum_e phi_e(x)^2
```

together with the logarithm of the number of exposed state rows. Importance-
weighted Bernoulli thinning can erase a positive fraction of the features,
and only a linear-size public list of masks is needed to cover every
coefficient signing. Sup-norm contraction then protects all shared max-type
futures.

This is a semantic compression law and a useful general upper bound. It is
not automatically a dynamic congruence: the weighted sparse surrogate may
leave the model class, and repeated copies of the approximated landscape pay
repeated error. Thus finite-dictionary sparsification quantifies response
image complexity while sharply exposing the separate realization and reuse
obligations.

## 40. Bounded extremal scale does not imply bounded response memory

The cap-`1/2` Walsh family has the smallest spectral upper scale supplied by
the classical Hadamard construction, yet its responses through one fixed
dense sign bridge contain `Omega(sqrt(n))` bits at `n^(3/2)` resolution. The
hidden state is a bent switching table; flat dual fields expose it, while a
low-bias condition keeps cross responses in a strictly smaller spectral
sector.

Thus a bound on the unperturbed maximum controls neither contextual entropy
nor dynamic memory. A positive compression theorem for bounded-cap
quadratics must impose an additional synchronization or realization law. The
Walsh example also shows what a useful strict state can look like: the
`sqrt(n)`-bit switching table is exponentially smaller than the full energy
landscape and sufficient for this structured family, even though no bounded
state exists.

## 41. Near-top entropy deficit is an information amplifier

For a Boolean landscape, the number of states within one target-scale unit
of the maximum is not merely a descriptive statistic. If that set has an
entropy gap `kappa n`, then almost-unbiased weighted sign tests supplied by a
random bridge can simultaneously avoid every one of its neighborhoods. The
coordinate-switching orbit consequently has a positive contextual
information rate.

This principle is deliberately one-way:

```math
\text{near-top entropy deficit}
\Longrightarrow
\text{rooted response packing},
```

not “extremal entropy is a sufficient state.” It explains in one theorem why
both a spectral-cap Walsh landscape and a high-cap pole landscape expose
linear switching information. The missing information is the root of the
rare extremal set, not its scalar cardinality.

At the exact Walsh cap the lower rate is matched by storing the switch, so
the response state is `Theta(n)` bits under coordinate-pinning futures. This
is still strictly smaller than quadratic coefficient data and exponentially
smaller than an explicit energy table. A positive composition theorem must
now show how that gauge label is transported or quotiented jointly with the
interface; a cap bound alone cannot erase it.

## 42. Relative operator algebra can be composition-created information

Static summaries of each child and pairwise label overlaps may omit how an
on-site operator sits relative to the bridge operator. In the Walsh family,
the missing invariant is the sign in

```math
FC=(-1)^\omega CF.
```

It is invisible to bias and within-word pair correlations but changes a
long bipartite composition by a fixed leading amount per block: commuting
involutions share transported maximizing eigenspaces, while anticommuting
ones obey a smaller Clifford-type spectral ceiling.

This is a stronger form of holonomy than an unbalanced scalar edge-sign
cycle. The edge connection itself can be flat while its compatibility with
the on-site extremal geometry carries a cocycle. A reusable state must retain
the relative algebra needed by future products, not merely each operator's
isolated spectrum or low-order correlations.

## 43. A scalar asymptotic state can exist without a periodic optimizer

Finite exact response carriers are not the only compositional mechanism. For
finite-width adversarial chains, submultiplicativity of minimum transfer-
product norms creates a scalar lower-spectral growth rate even when no finite
word attains it. Under strict positivity, a projective cavity potential and
contraction refine that scalar into a polynomial approximate mean-pressure
carrier.

This separates three claims which should not be conflated:

1. an asymptotic scalar value exists;
2. a finite or compact state computes it approximately;
3. one bounded-period microscopic construction realizes it.

The first two can hold while the third fails. Conversely, the mechanism is
controlled by interface width: at a dense quadratic cut the bilinear response
is already of leading order and the local transfer dimension is exponential.
A useful dense analogue would need a nonlocal multiplicative quotient, not a
larger finite-width transfer matrix.

## 44. Relative algebra gives certificates before it gives a quotient

For involutory child and bridge operators, commutators and anticommutators
control two different extremal mechanisms. A small commutator preserves a
transported common section with quadratic loss, while a small
anticommutator imposes a global spectral ceiling. This information can force
a leading compositional separation even though it does not compute the
intermediate Boolean maximum.

The distinction is methodological. A low-dimensional invariant can be a
powerful **certificate state** without being a semantically sufficient
response state. Promoting the former to the latter requires a realizability
law for Boolean poles and closure under varying children.

## 45. Orbit sufficiency and compositional congruence are different

For linear Walsh labels in every ambient dimension, a finite presented form--binary
Gram data, all label relations, and one characteristic-root fibre--classifies
every ordered tuple up to ambient orthogonal symmetry. It therefore answers
all equivariant Boolean extremal queries on the tuple. This is a strict
quotient of the raw labels when the word is short relative to the ambient
dimension.

Yet two such isolated orbit states do not determine their joint state:
cross-pairings and mixed relations appear only after gluing. Exact orbit
classification is therefore a static response theorem, not automatically a
dynamic congruence. The missing cross-form measures composition-created
information in a concrete algebraic family.

## 46. Accumulated-span memory is the closure of relational gluing

For rooted bilinear presentations, the relative gluing datum has three
logically independent parts: a cross form, an intersection correspondence,
and a root fibre. Their pullback formulas are associative once the composite
retains its accumulated presented span. Pairwise edge data is not closed:
one ternary dependency can be invisible on every pair.

This gives a quantitative version of composition-created information.
Cross-form memory can cost `rs` bits and intersection memory `r^2+O(1)` bits,
even when isolated states are fixed. The resulting carrier is still a strict
symmetry quotient in the low-word-length regime. Thus "new information" need
not mean reconstructing the microscopic landscape; it can mean updating a
small but growing relational presentation.

## 47. Interface rank is relative to the declared error scale

Algebraic rank counts even singular directions whose entire Boolean
contribution is below the target error. The relevant linear interface is the
number of singular values above `L_n/n` at objective scale `L_n`. Operator-
norm truncation then gives an all-state pointwise replacement, and the
retained directions feed the exact upper-roof algebra.

This does not erase the response/congruence distinction. Quantizing the roof
answers only bounded fields through retained features, and local low-rank
factors can create a large separator state after elimination. Numerical rank
therefore predicts **local response image**, while dynamic memory still
depends on the composition graph and allowed future queries.

## 48. Rare-event compactness must charge descriptor complexity and speed

Hypograph convergence of exact-fibre log counts preserves the best recovered
fibre, but total mass can hide in exponentially many distinct descriptor
values. A valid probability or pressure theorem therefore needs either a
subexponential descriptor image, a coarse-bin Laplace principle, or local-
ball entropy. This is the rare-event analogue of charging query complexity.

At a fixed exponential speed, supremal convolution is closed and finite-rate
branches survive. Subexponential multiplicities and spacings do not. The
choice of speed is therefore part of the state, just as the allowed future
query family is part of a contextual response quotient.

## 49. Complete spectral data can be a strict unrooted quotient

In the linear Walsh family, all power traces conditional on a marked tuple
factor through its Gram form and relation kernel. The characteristic-root
fibre, although necessary for the full orthogonal orbit and visible to a
rooted Boolean query, disappears from every such unrooted spectrum.

Thus “all moments” is not synonymous with “all extremal information.” The
missing datum is not a higher spectral moment but a root supplied by the
future interface. Any proposed moment hierarchy must declare whether its
queries carry such roots or other coordinate anchors.

## 50. Declared roots are symmetry stabilizers, not intrinsic coordinates

The same structured family can have two exact orbit states because its query
languages have different symmetry groups. For unrooted Walsh graphs, the
full ambient orthogonal group makes the label-space characteristic root
invisible and `(Gram,relations)` conjugates the whole landscape. A fixed pole
or coordinate field restricts the group and makes the root fibre observable.

This gives a practical reduction rule: compute the symmetry group of the
**complete declared continuation**, not merely the coefficient
parameterization. Information surviving only under a smaller, artificially
rooted group must be omitted from an unrooted carrier. Conversely, gluing
requires compatibility data only on label sets that can occur together in
one connected future; local orbit charts can use independent gauges across
disconnected supports.

## 51. Approximate congruence is weighted by interaction mass

When exact semantic charts are available on connected pieces, one can build
an approximate reusable state by omitting weak interactions.  The natural
error is not the number of omitted compatibility variables but the uniform
oscillation they can contribute.  For Walsh bridges this is exactly the
deleted absolute edge mass times `n^(3/2)`.

This principle has two limits.  It gives linear memory on bounded-incidence
graphs at extensive error, but a dense graph with uniformly nonnegligible
edges keeps quadratic component mass under the same deletion architecture.
It is therefore an approximation law, not a claim that all weakly encoded
holonomy can be ignored by a more global decoder.

## 52. Cycle flux may be scalar semantic information

A compatibility cocycle need not require a rooted or vector-valued query to
be observable.  In the linear Walsh family, an ordinary triangle maximum
distinguishes two tuples that agree in every self-pairing, relation, and root
flag but have opposite off-diagonal Gram flux.  Thus the correct unrooted
quotient may discard coordinate anchors while retaining cycle holonomy.

The open quantitative question is how many independent fluxes a family of
scalar future queries can expose at fixed distortion.  Orbit-state dimension
alone is only an upper bound; scalar visibility must be established by
separating queries.

## 53. A scale-preserving witness lift creates response compactness

An exact microscopic optimizer need not be consistent across scales. It is
enough that every Boolean witness at one scale has a feasible lift at the
next scale preserving all declared normalized response coordinates. Then
the convex response images are nested, and compactness turns their increasing
union into a limiting state whose support functions converge.

Regular-Hadamard amplification realizes this principle exactly. The useful
object is the whole fixed-port cross-correlation set, not one maximizing
sequence. Approximate applications require a quantitative lift or retraction
with summable response defect; assuming Hausdorff convergence itself would
merely restate the desired conclusion.

## 54. The finest retained speed must pay branching entropy

A finite hierarchy of exponential scales composes lexicographically only
while the logarithm of the competing decomposition count is negligible at
the smallest retained scale. This is a worst-case exact boundary: tied terms
convert branching multiplicity directly into the last response coordinate.

Smooth leading entropy does not remove the issue. A saddle can have
polynomially many near-maximizers whose logarithmic mass survives after all
pointwise coefficients have been recorded. Any genuinely multiscale compact
state must therefore carry tangent density, Hessian data under suitable
regularity, or another object that accounts for exposed-fibre mass.

## 55. Directed recovery, not reverse reconstruction, is enough for a limit

If every old response point has a common Boolean lift to the next scale with
summable all-pairs distortion, the response sets converge even though new
points may appear indefinitely. Compactness converts one-sided approximate
nesting into a limit; a backward map is needed only for a quantitative
innovation-side rate.

The all-pairs qualifier is essential. Preserving only self-quadratics can
lose cross-block queries, and separate optimizer-dependent lifts do not
define one composable response map. Finite total drift is also a real
threshold absent additional cancellation: vanishing square-summable steps
can still trace a bounded nonconvergent scalar path.

## 56. Connectedness does not erase independent compatibility flux

Compatibility memory is not merely an artifact of disconnected local
queries.  Independent Walsh relation triangles retain independently
observable off-diagonal Gram fluxes even when every query uses one common
connected support.  Arbitrary nonnegative bridge padding is harmless when a
favorable witness saturates every added bridge and the same operator bound
caps it in all competing states.

The scaling qualification is part of the principle.  This proves additive
memory per flux at a fixed port scale; it does not by itself prove an
extensive information rate after the number of ports is included in the
normalization.  A general dense-interface lower bound must control both the
number of exposed fluxes and their response gap at the total-system scale.

## 57. Recovery can leave a continuous scale phase

Exact scale-preserving witness lifts make every fixed outer tensor template
converge, but they need not identify different outer prefixes.  For any
regular-Hadamard automatic prefix hierarchy, the residual all-order state is
a continuous base-`h` mantissa profile.  Convergence of the scalar response
is equivalent to constancy of that profile, not a consequence of phasewise
compactness.

The Walsh hierarchy shows this residual state can be nontrivial by a fixed
gap.  Hence an all-order recovery theory needs a phase-synchronization or
phase-selection principle in addition to realizability at every order.  The
phase is still a strict compact quotient of the Boolean landscape, so the
counterexample diagnoses the missing information without reinstating the
full optimization.

## 58. Residual synchronization does not imply witness synchronization

A bounded reset delay can collapse all normalized max-plus rows to finitely
many residual profiles.  This is enough for exact rooted response updates and
an exact cycle-potential law.  It is not enough to select one support per
residual phase whose individual edges attain the residual scalar toll.

The greatest incoming-image core exactly decides that stricter fixed-context
presentation.  Its failure is certificate failure unless the query exposes
the missing endpoint geometry.  This separates three resources cleanly:
semantic residual state, locally path-realizing witness support, and the
terminal query language that can distinguish them.

## 59. Phase mixing is useful only relative to its transfer defect

Cross-scale recovery can be promoted to phase synchronization when the
recovery maps themselves mix phase labels.  The controlling quantity is not
full support alone but

```text
uniform recovery error + accumulated operator defect
----------------------------------------------------- .
                 refreshed phase mass
```

A vanishing ratio triggers a maximum principle and collapses every continuous
limiting phase response.  A fixed positive ratio can preserve a macroscopic
phase forever.  This gives a checkable positive counterpart to the automatic
Walsh phase obstruction without assuming equality of phase carriers or their
maximizers.

## 60. A tangent law becomes a theory only on a closed regularity class

The polynomial mass of an exposed convolution fibre is part of extremal
response.  At a nondegenerate `d`-dimensional saddle it contributes `d/2` to
the logarithmic exponent and a Hessian/lattice-density amplitude.  Pointwise
multiscale roofs omit this term.

Arbitrary tangent profiles still have unbounded functional information.  A
genuine compression theorem requires a class closed by the tangent law.  The
Gaussian class supplies one: mean and covariance add, mass amplitude
multiplies, and finitely many parameters realize every future convolution.
Quartic and flat saddles mark the boundary rather than inviting an
unrestricted hierarchy of formal coefficients.

## 61. A public interaction cannot amplify state it does not contain

Optimization may create joint cancellation, but it is still one-Lipschitz in
the underlying landscape.  If hidden information enters only through
bounded local terms, every state-independent interaction cancels pointwise
when two hidden states are compared.  The total response diameter is then
controlled by the sum of local oscillations, not by the size or density of
the public coupling.

At the `N^(3/2)` scale, disjoint `n`-variable ports lose a factor
`1/sqrt(k)`.  A compatibility bit can remain macroscopic only if it is
broadcast into a superlinear number of state-dependent cross terms, its
normalized coefficient grows, or the query reads a different local/vector
scale.  This is a general composition-created-information threshold, not a
Walsh-specific spectral fact.

## 62. Averaged scale statistics are query quotients, not synchronization

A continuous discrete-scale phase always has a unique logarithmic empirical
law, but that law forgets the cyclic ordering of response values.  Positive
power-weighted futures retain the phase and even reconstruct its entire
profile from their subsequential mean function.

Hence existence of a canonical averaged limit proves only sufficiency for
that averaging query.  It cannot replace an all-order recovery or phase-
mixing theorem when later composition is sensitive to the newest scale.

## 63. Forgetting the endpoint profile does not forget accumulated reward

A bounded-delay near-reset pays terminal residual error only once.  The
scalar gained when one residual hands off to the next is a directed
compatibility cocycle, and any non-potential cycle defect is paid at positive
rate under repetition.

This distinction is scale-sharp: an arbitrarily small one-profile shell can
contain a scaled copy of an arbitrary finite weighted response algebra.
Dynamic compression therefore needs two certificates--terminal-state
forgetting and reward-cocycle compatibility--not merely a small static image
or a strong projective contraction.

## 64. Mixing a semantic phase has an inverse-defect information price

Uniform stationary scrambling turns an observable excess over the phase
average into a quantitative alternative: either the transfer accumulates a
comparable toll over one mixing window, or the phase space grows
exponentially in the inverse one-step toll.  This converts local mixing time
into a lower bound on reusable response memory.

Gauge labels are exempt because their responses were never distinct.  A
large coordinatewise operator orbit can refresh with constant description
when signed conjugacy makes it semantically one state.  The lower law charges
response excess, not syntactic phase count or operator diameter.

## 65. Spectral rank and response rank are different resources

At bounded `sqrt(n)` operator scale, a dense sign bridge has linearly many
singular directions visible at fixed `n^(3/2)` accuracy.  Thus spectral
truncation cannot provide a subextensive bridge state there.

This does not preclude compression by symmetry or synchronization.  A valid
low-information dense interface must quotient the visible directions
nonlinearly; algebraic or numerical rank alone cannot explain it.

## 66. Closure of the leading roof does not imply closure of tangent mass

Infimal convolution may preserve a finite family of rate functions while
the normalized saddle shapes generate an ever richer convolution family.
The power-exponential example makes the separation exact, and Gaussian
stability identifies the exceptional closed tangent law.

Multiscale extremal compression must therefore test every response speed
used by future composition.  A closed leading exponent is not evidence that
the next entropy or amplitude coordinate has finite state.

## 67. Total-scale compatibility memory requires broadcast incidence

Public interactions cannot amplify hidden state absent from their
coefficients, but state-dependent dense coefficients can.  An alternating-
form evaluation code broadcasts only linearly many hidden bits across a
quadratic number of constrained edge phases and attains a positive
`N^(3/2)` response rate while every child remains spectrally flat.

Hidden description length alone is therefore not an incompressibility
parameter.  One must also charge how widely each hidden coordinate changes
the interaction, as well as shared public advice and the declared context
language.

## 68. Broadcast information is limited by total hidden incidence

For a fixed public future language, optimization cannot make one bounded
atom respond to more hidden coordinates than its coefficient actually
reads.  Summed contextual influence is bounded by the total hidden
coordinate--atom incidence.  With only `O(N^2)` bounded pair atoms, fan-in
`t` therefore supports at most `O(t sqrt N)` coordinates having an
`N^(3/2)` neighbouring response gap.

This law is sharp for unrestricted quadratic signings, including a matching
pairwise response packing.  Spectral flatness and concise public advice are
separate resources.  A flat positive-rate family can evade the local regime
only by making its hidden coordinates densely incident, by changing the
normalization, or by leaving the bounded fixed-child presentation.

## 69. Forgetting scale and cycle memory add; they do not substitute

On a certified common-law contracting carrier, centred terminal information
is forgotten geometrically and fresh centred reward error accumulates only
to its geometric resolvent.  The invariant scalar reward direction is a
different resource: its persistent response is exactly recurrent graph
cohomology, measured by directed cycle means or stationary flows.

Thus dynamic response complexity separates into a static dictionary sampled
at the forgetting scale and a cycle-space memory sampled at the target rate.
More contraction cannot compensate for scalar holonomy.  Conversely, once
the cohomology class is retained, arbitrary switching inside the contracting
fibres creates no additional extensive expected-reward state.

## 70. Public randomness and semantic broadcast are different resources

An exponentially large family may be simultaneously spectrally controlled
by a short pseudorandom seed when the proof tests only bounded moments.  Its
semantic response rate can remain extensive because hidden coordinates are
broadcast through the resulting dense coefficient code.

Therefore a quadratic public edge table is not intrinsic to dense response
memory.  One must separately charge hidden incidence, seed/description
length, query-language size, and closure under the declared composition.
Derandomizing one resource does not silently solve the others.

## 71. Persistent broadcast must enter the invariant dynamic channel

On a fixed strictly contracting carrier, hidden terminal and centred-reward
coordinates have only finite forgetting-scale value.  Their all-depth
response cannot grow linearly.  The exact positive-rate channel is the
visible scalar reward cohomology, and its total hidden capacity is bounded
by the incidence of atoms that alter that cohomology.

Mixing therefore does not turn dense compatibility information into free
memory.  It either forgets the information, pays it repeatedly as fresh
reward, or exposes it as a cycle coordinate.  Static response compression
and dynamic persistence are linked by this resource conversion, not by one
interchangeable complexity number.

## 72. A compiler must preserve both semantics and normalization

An exact optimization identity can be useless after composition if its
auxiliary order changes the macroscopic normalization.  For unit quadratic
systems, compiling a `k^(3/2)` response at positive total scale requires
`N=O(k)` total vertices, not merely a finite exact realization.

Independent selectors make this tradeoff quantitative through their pair
Fourier mass: bounded fan-in forces superlinear order, and complete fan-in
forces quadratic order when both orientations are required.  Correlated
selectors replace that atomwise obstruction by an exposed-set covering
problem.  Hence the next invariant is not “number of auxiliary variables”
alone but the joint triple of total order, slope geometry, and entropy of
the source states that future optimization must expose.

## 73. Exact contextual closure can hide a macroscopic calibration

A query-dependent exact-sign future can expose every coordinate of a child
landscape with only linear order and preserve its entire response metric.
The price may be a common quadratic energy baseline that cancels between
responses but remains visible to an absolute ground-state objective.

This distinction is structural.  A fixed future robustly pinning one state
against the full complete-sign child cube must have quadratic oscillation and
quadratic parent cap.  Consequently a response compiler should be charged
separately for semantic distortion, order blow-up, and absolute calibration.

## 74. Extremal entropy is orientation-sensitive

Mean-zero low-degree structure forces each **signed** side of a nonzero
quadratic landscape to occupy constant cube mass.  Once one extremum is
oriented as a minimum, this gives a full-entropy antipodal exposed set and a
selector covering lower bound.

Taking absolute values can destroy that conclusion even under spectral
flatness and target-scale maximum.  A compiler that pays both orientations
separately may therefore face linear information cost while a genuinely
joint absolute channel evades that particular proof.  Same-switch
cancellation is an algebraic resource, not just a constant-factor issue.

## 75. Exact pullback rigidity can expose a locking no-go

Requiring a Boolean encoding to preserve an entire algebra of pair queries
may leave no nonlinear compression: the only maps are a coordinate
permutation, signs, and a common gauge invisible to those queries.  Once this
rigidity is proved, a freely optimized exact-sign bridge has a leading
advantage over every intended copy.

This does not invalidate query-dependent contextual compilation.  It says
that a proof which first transports the complete query algebra to a second
shore and then pays a dominant public equality lock has already lost the
target scale.  Any escape must preserve a smaller query algebra, correlate
the auxiliary selectors, or make cancellation occur before the lock is
optimized separately.

## 76. Bounded cap charges response bins per witness, not fibre size

Exponentially many bounded-cap children may share one optimizer under one
fixed future.  What one witness cannot carry is exponentially many
macroscopically separated scalar responses: its evaluation lies in a bounded
interval with only constantly many target-scale bins.

Accordingly, a response packing forces an extensive approximate witness
dictionary, but that cost has two realizations.  It may be a public library
of common query pins, or it may be child-dependent optimizer switching.
Any theorem claiming the latter must first rule out the former.

## 77. Extremal tail mass controls orbit-query complexity

For a transitive translation language, preserving a threshold extremum is a
set-transversal problem.  The reciprocal mass of one extremal witness set is
the query complexity up to a logarithmic factor.  Under independent product
composition, its exponent is the Legendre dual of the additive log-moment
state.

This is a genuine bridge between rare-event and contextual theories.  It
closes exactly for product landscapes and covering codes, while
overlap-dependent interaction is the structural event that invalidates the
one-point rate state.

## 78. A bridge penalty must dominate the local modulus, not merely the cap

An orthogonal sign bridge supplies an exact quadratic penalty for leaving its
Boolean pullback relation.  Whether that relation survives a future is
controlled by the future's modulus near the relation, not by its global
Boolean maximum alone.  A cap-scale child may concentrate its entire gain on
a `sqrt(n)`-coordinate departure and defeat the relation at leading scale.

Thus a synchronization theorem needs a pullback cancellation, local-field or
operator bound, or another alignment hypothesis.  “Both pieces have the
correct cap” is not a compositional invariant.

## 79. Metric compilation is weaker than dynamic congruence, but genuinely so

A context family may preserve an exponential contextual metric with bounded
cap even when no query reconstructs a coordinate response and no fixed
witness is pinned.  Child--query-dependent optimization can create an
anti-diagonal gap whose two signs are recovered only after taking the
projective oscillation across queries.

This makes one-layer metric compilation a real intermediate resource, not a
failed form of pointwise simulation.  Its next obligation is reusable
closure: after another interaction, the response quotient must still be a
congruence.  A metric embedding without such closure cannot by itself support
cross-scale recurrence.

## 80. Port features compose by amalgamation, not by product

A finite port family may admit a small Gram or resolvent state, while two
such states cannot be glued from their separate coordinates alone.  The
cross pairings are new semantic variables.  If the model supplies no
synchronization law determining them, every amalgamation enlarges the state
even though each isolated port remains simple.

This is composition-created information in its most elementary form.  A
claim of reusable feature closure must specify the cross-data map; taking the
cartesian product of marginal feature states silently sets precisely the
unknown part to zero.

## 81. Visibility has an interface-scale threshold

Whether a quotient is valid can depend on how much future system is allowed
to interact with it.  For quadratic absolute caps, the orientation
`A` versus `-A` is indistinguishable up to twice the internal cap of the new
shore.  It is therefore invisible to every continuation of width
`o(n^(3/4))`, yet a dense shore of comparable order exposes it at the full
`n^(3/2)` scale.

Approximate contextual equivalence should consequently be indexed not only
by error but also by a continuation budget.  A state can be a true
congruence for one scale class and false for another without contradiction.

## 82. Interface size and interface cap are different resources

A sublinear continuation may still carry a leading old-scale response if its
internal landscape is supercritical relative to its own order.  For
quadratic sign systems, width `n^(3/4)` supports a clique cap `n^(3/2)` and
is already enough to expose a hidden orientation bit.  Requiring every
component to have its natural `m^(3/2)` cap moves the same threshold to
linear width.

Compositional complexity should therefore charge an interface by both its
number of variables and the extremal budget it can spend internally.  An
order count alone can mistake a small but high-cap calibration device for a
negligible future.

## 83. Reconstruction dimension can exceed macroscopic information dimension

An exact feature algebra may store every pairwise compatibility coordinate,
yet positivity can prevent those coordinates from varying independently at
fixed amplitude.  In a Gram/Rayleigh state the exact table is quadratic in
the port count, while every full independently toggleable macroscopic cube
has only linear dimension.

Counting entries is therefore neither an upper nor a lower law for reusable
information.  The operational question is the metric entropy of the
realizable image under declared collective responses; dense correlated
small coordinates must be tested separately from sparse independent ones.

## 84. Quotient gluing has a compatibility fibre

Compressing each component modulo its internal symmetries does not make
composition cartesian.  The joined object lies over the tuple of marginal
quotients with a fibre measuring the relative frames and new holonomies that
the marginals discarded.

In the regular-Hadamard orientation model this fibre is exact: `r` joining
edges create `s-1` relative antipodes and `r-s+1` cycle bits.  This suggests
a concrete dynamic-complexity law: reusable state equals marginal semantic
state plus the dimension of the amalgamation fibre, with the latter charged
only when interaction makes it observable.
