# Extremal landscape information

This directory starts a theory-building program about a question that occurs
across combinatorial optimization, coding theory, and disordered systems:

> What information about an exponentially large energy landscape is necessary
> and sufficient to preserve its extreme values under future perturbations or
> composition?

The motivating example is a quadratic Boolean landscape

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
```

but the definitions here are not tied to the open problem about
`M_n/n^(3/2)`.  This directory deliberately does not modify the main project
state.  It asks for general theorems first and will reconnect to that problem
only after the framework has a nontrivial information inequality, a
composition theorem, an application outside sign matrices, and a credible
finite-realization theory.

## Present status

The program has reached **Level 2 (explanatory theory)** and has several
Level-3 generative results.  It is not yet a finished theory and it is not a
proof of the motivating convergence statement.

The first investigation gave four precise conclusions.

1. For homogeneous quadratic Boolean forms, positive-rate upper-tail entropy
   reaches the maximum.  A Hamming-noise cloud around any maximizer contains
   exponentially many states at every fixed strict fraction of the maximum.
   Consequently two such sequences cannot have the same resolved upper-tail
   entropy but different limiting normalized maxima.
2. Entropy is not compositionally sufficient.  Two weighted quadratic
   landscapes can have exactly the same finite energy histogram and the same
   complete energy--energy--global-overlap histogram, yet have fixed-block
   coupling responses separated by `1/2` asymptotically.
3. Pair data can fail for the same reason in coding theory.  Two explicit
   code landscapes have identical exact energy--energy--distance counts at
   every tensor power but covering radii `2r` and `3r`; a rooted intervention
   separates them by order `r`.  More generally, for every fixed replica
   order `k`, explicit code pairs agree in their complete unrooted data through
   `k` points while Cartesian powers retain a positive normalized
   covering-radius gap.
4. Preserving all counterfactual optima can be information-heavy even when
   preserving the unperturbed maximum is free.  There is a
   `2^binom(n,2)`-element family of quadratic landscapes, all shifted to have
   maximum zero, for which mean-square accurate pinned-query answers require
   `binom(n,2)[1-h_2(D)]` bits.

These results force the current candidate concept: a summary must be declared
relative to the **query interface** through which a future environment can
interrogate the landscape.

The second investigation tested whether that principle generates mathematics
rather than only terminology.  It produced three core project-level results:

1. a sharp posterior-width information inequality for an arbitrary Boolean-
   indexed Hilbert response embedding, with exact applications to quadratic
   responses, nearest-code responses, Max-Cut, and boundary kernels;
2. a syndrome-support semilattice for binary parity-check fragments: it is a
   strict quotient of the code, composes exactly under concatenation, answers
   every unrestricted appended-fragment covering-radius query, and has sharp
   worst-case exact response complexity `Theta(2^w)` bits.  A block-direct-sum
   subfamily also has a positive macroscopic rate: for every fixed
   `epsilon<1/8`, error `epsilon*w` requires `Omega_epsilon(w)` bits, matched
   by an `O(w)` exact state on that family; and
3. a finite deterministic synchronization theorem: approximate mixture
   ultrametricity plus a cross-root monotone-linkage condition forces every
   species overlap to be a uniformly controlled scalar function of total
   overlap.  A rare-fibre construction proves that averaged conditional
   variance cannot replace the linkage condition at zero temperature.

Two additional boundaries are now quantitative.  Shapley--Folkman
convexification controls all Lipschitz aggregate queries when the effective
response rank is fixed, but a same-zonotope example has a leading discrepancy
gap when that rank grows.  Four-cell tropical crossings protect exact channel
count at lattice-scale uniform error, but canonical code transversals lose
exponentially all such exposure under diffuse normalized mean-square queries.

The third investigation moved from isolated examples to a composition law.
For arbitrary syndrome supports it proved both a subexponential approximate
response net and an exactly closed hard-core quotient; a geodesic/anticode
argument gives the latter only
`2^((1-2 epsilon)w+O(log w))+O(w^2)` description bits at error `epsilon w`.
More importantly, it identified the missing information created by
composition.  In a binary group extension, kernel offsets are removable
gauge on quotient-independent columns.  They become observable precisely on
new mixed quotient cycles, where their holonomy controls the extremal word
response.  The mixed relation space has dimension

```math
\kappa=\sum_j\dim U_j-\dim\sum_jU_j,
```

and its labeled gluing freedom is exactly `Hom(F_2^kappa,F_2^D)`.  Circuit
defects amplify by at most cycle-space nullity, with equality on disjoint
blocks; bounded-arity composition tests fail at every arity; and even one
new mixed channel carries `Omega(D)` rooted-response bits at fixed
macroscopic accuracy.  Conversely, cycle contraction synchronizes a dense
source to a linear graph at constant all-context error, while exact linear
graph sources of affine rank `r` have a closed feature state with `O(r)`
response error.

The fourth investigation forced the framework through classical benchmark
problems instead of extending its vocabulary.  It produced the following
conclusions.

1. Contextual response equivalence independently recovers the standard
   separator table, zero-temperature transfer state, fixed-rank aggregate,
   and weighted-language residual for the operational reason that every
   allowed future factors through that interface.
2. Pure weighted Max-Cut futures expose exactly the projective boundary
   profile.  More strongly, a private lookup construction realizes every
   projective table, so contextual exposure is matched by language
   realizability rather than inferred from dynamic programming alone.
3. Under unit boundary load, the exact Max-Cut response class is still the
   entire projective one-Lipschitz ball.  Its radius-`epsilon w` covering
   number has logarithm between
   `2^((1-H_2(2epsilon)+o(1))w)` and
   `2^((1-H_2(epsilon)+o(1))w)`.  Local sensitivity is therefore not a
   compression promise.  Polynomial component size is a genuinely stronger
   resource, quantified by the shared-parameter bound below.
4. Finite-metric distance kernels form a nonproduct interacting algebra:
   strengths combine by a bottleneck minimum and isometry labels by holonomy.
   Anisotropic projective-Hamming shells strengthen this to a transported
   coordinatewise-minimum lattice.
5. The same distance shell is the exact nearest Lipschitz response to an
   arbitrary profile.  Repeated shells pay precisely the weakest-layer
   Lipschitz defect, not one approximation loss per layer.
6. Approximate idempotence recognizes a nearby metric with a sharp
   dimension-free repair, but it is not by itself stable under long
   composition.  A bounded-diameter line family has vanishing one-step
   defect and fixed long-depth response drift because a small
   per-transition toll accumulates.
7. Global presentation size supplies the missing resource bound.  Unit-load
   Max-Cut components with `m` edges have only
   `O_epsilon(m^2+m log(w+m))` response bits at error `epsilon w`, even with
   arbitrary real weights.  Approximating the full unit-load response ball
   therefore needs exponentially many edges.
8. The weighted-automaton benchmark yields a second strict composable state:
   tropical lumpability replaces `p` raw coordinates by `r<p` nonlinear
   block maxima, closes at arbitrary depth, and has matching
   `Theta(r log(1+B/epsilon))` response complexity when quotient suffixes
   expose every retained coordinate.
9. Depth-uniform approximate reuse depends on the perturbation model.
   Compatible interface gauges telescope and recurrent small-diameter
   continuations reset old error.  For fixed coherent maps, exact finite
   semigroup relations also give bounded normal forms: nearby idempotent
   clamps falsify the narrower “kernel gauge or small full-image reset”
   dichotomy.  Against fresh arbitrary residuals, however, syndetic tangent
   resets are quantitatively necessary and sufficient on selector cells;
   recurrent translations obey a twisted cycle-mean criterion.  For every
   regular tie-free selector language, the minimal kernel-partition lift
   decides which side holds and returns either a uniform bound or a pumpable
   drift cycle.
10. Static response complexity now has a cross-benchmark two-sided law.
    Query landmarks bound response covers from above, while balanced robust
    exposure bounds them below.  For the full Lipschitz language the exposure
    dimension is exactly the packing number of the query interface; for
    weighted automata it is controlled by exposed coordinates and orthant VC
    dimension.
11. Shared-parameter presentation has a nontrivial robust lower bound.
    Facets of high-facet `0/1` polytopes produce tie-free binary max-affine
    response families with `Omega(m log m)` bits at fixed distortion.  The
    general radius-bounded upper remains `O(m^2)`; raw normal-fan cells do not
    settle that gap.
12. The heterogeneous mean-field benchmark now closes exactly.  A scalar
    chemical potential exposes the discrete-concave top-occupancy profile;
    its slopes merge by multiset union, and a common quantized histogram has
    simultaneous sublinear bits and depth-independent sublinear error, for
    example with grid scale `eta=B/sqrt(N)`.  A fixed bilinear
    pair interaction can strictly synchronize this profile to its concave
    envelope; at mass `n`, `J>=4B/n` sharply collapses the state to total
    field alone.
13. A cross-benchmark feature-algebra theorem now predicts the polynomial
    branch rather than merely cataloguing it.  For `d` additive response
    atoms, the exact mass-`n` quotient has `Theta(n^r_Z)` states, where
    `r_Z` is the arithmetic rank of the atom-response differences.  Robust
    conditioning gives two-sided response covers; on an equally spaced
    mean-field grid the exact cover below half-grid error is
    `{n+d-1 choose d-1}`.  This complements the separator's exponentially
    exposed lookup carrier and the automaton's fixed lumped carrier.
14. The polynomial branch extends beyond finite alphabets.  If individual
    atom responses have an `eta`-net of size `D(eta)`, their type histogram
    is an exact additive update with error `n eta` and at most
    `{n+D-1 choose D-1}` states.  Whenever
    `D(eta_n)=o(n)` for some `eta_n->0`, both the summary rate and
    normalized response distortion vanish without a merge-depth penalty.
    One common root-scale net must be used throughout the merge tree; an
    adaptively refined parent cannot recover types collapsed by a child.
15. Switching and ties now have an exact fresh-residual law. A max-plus
    secant is row-stochastic, and its worst accumulated projective error is
    the sum of total-variation separations of terminal suffix rows. Uniform
    block scrambling gives depth-independent error even without a finite
    tangent reset. This does not yet decide which secant paths are
    dynamically realizable by fixed coherent kernels.
16. For coherent switching, an exact regular affine-selector presentation
    has a finite necessary-and-sufficient witness-cycle test: positive cycles
    pump directed drift, while zero cycle holonomy characterizes two-sided
    boundedness. Different left/right selectors force the joint cross-
    difference carrier. The earlier suffix-set reset state also compresses
    sharply to the kernel partition and is worst-case minimal at a fixed
    control vertex for the full selector alphabet. A naive local
    face graph is not path-realizing, so finite tropical lumpability—not more
    cycle notation—is the remaining structural obligation.
17. The forward finite-lumpability obligation now has an exact criterion and
    a checkable nontrivial subclass. Iterated pullback of the observation
    partition stabilizes exactly when a finite deterministic path-realizing
    quotient exists. A finite oriented affine
    arrangement closed under branch pullback is a checkable certificate; for
    selector coordinate differences, its unrestricted closure is finite
    exactly when every reachable unit-transport cycle has zero affine
    holonomy. Normal directions without offsets or orientations do not
    suffice.
18. Approximate block lumpability plus uniform switching contraction gives a
    finite depth-independent response simulator without active-cell
    enumeration. More generally, static response entropy and dynamic memory
    interact through scale: suffix gain changes the resolution at which an
    internal cover is paid, while `q`-ary branching transforms local response
    entropy across the context tree. A finite affine tree shift attains the
    resulting bound at every level.
19. The converse cannot be only a drift-cycle theorem. A compact five-piece
    encoder--decoder has bounded one-step response entropy but exponential
    horizon memory; irrational rotation has unbounded exact refinement
    without amplitude drift. On the positive benchmark side, arbitrary
    strict-strip futures expose the complete Ising boundary table with an
    exact sup-metric isometry, while a restricted rational width-two alphabet
    genuinely minimizes from seven reachable profiles to two weighted states.
20. Compact rational unit-selector systems now supply a genuine finite
    tropical lumpability theorem. Selector pullbacks have a finite normal
    orbit; rational offsets form a lattice; compactness leaves only finitely
    many relevant translated walls. Rational LP saturation therefore builds
    an exact path-realizing quotient for every finite polyhedral observation
    coloring. Irrational rotation proves that compactness without a discrete
    offset group is not enough.
21. Every finite predictor has a canonical future-behavior pseudometric. Its
    transition maps are nonexpansive, its failure to recouple with physical
    re-encoding is at most twice its response error, and a strict behavioral
    gap promotes this to exact semiconjugacy. Quantizing this metric gives the
    sharp law `new error <= old error + net radius * suffix gain`.
22. A finite weighted control graph makes the dynamic factor checkable:
    depth-uniform reuse holds when every reachable unit-coefficient cycle is
    absent, with constants from max-times path products. A saturating-chain
    family has a constant-size static cover but linear optimal predictive
    memory, proving that the mixing factor cannot be removed.
23. The rational compact theorem extends to accumulated affine rewards after
    a finite lift by realized affine germs. Vanishing of the nonconstant part
    of every lifted cycle label is equivalent to a finite-state scalar-toll
    simulator with depth-independent error; failure returns a genuine
    pumpable cycle. Transient merging diamonds can obstruct an exact potential
    without causing drift, so raw periodic-orbit tests do not characterize
    exact cohomology in noninvertible systems.
24. Repetition turns cycle-mean geometry into an information lower bound:
    any simulator with per-step error `epsilon` needs at least the
    `2epsilon`-packing number of the exposed cycle-mean image. A continuum of
    means therefore forbids finite bounded-error memory. On an arbitrary
    finite path-realizing germ graph, the maximum average separation over all
    reachable cycles is a pseudometric whose packing and covering numbers
    give matching lower and upper dynamic-memory bounds up to the finite
    visible-control factor and a bounded transient error.
25. Deterministic continuous rational selector systems with fixed whole-fibre
    control transitions also admit finite invariant projective grids at
    arbitrary resolution. Exact evolution of one nearby grid point shadows
    every switching word forever by nonexpansiveness, so terminal-response
    error is depth-independent without contraction or repeated rounding.
    State-dependent target controls require prior exact refinement; this
    arithmetic phase-locking is also absent for irrational rotation.
26. Finiteness is not compression. With two permutation selectors and one
    repeatable coordinate probe, a fixed three-letter system requires exactly
    `binom(r,floor(r/2))` states below per-step error `1/2`, despite compact
    rational isometric dynamics and globally unique branches. One selector
    alone has only Landau-subexponential recurrent iterate complexity.
27. A candidate finite reward congruence has an exact cycle-mean LP for its
    best asymptotic distortion. Such feasible congruences need not have a
    unique coarsest member: pairwise asymptotic response equivalence can lose
    the cycle incidence created by merging controls.
28. Finite projective max-plus semigroups give an intrinsic all-word carrier.
    Their synchronized weighted Cayley graph realizes the exact extreme
    aligned-word spectral gaps as finite graph cycle means; a letterwise
    critical graph or switching envelope does not.
29. A single `O(r)`-gate continuous rational lattice-PWA selector map can
    implement a `2^(r/2)`-period binary counter. Constant-weight block digits
    raise this to `2^(r-o(r))` exposed phases with a polynomial-size circuit.
    One identity probe exposes exactly all phases below error `1/2`; without
    that probe, phase rewards differ only by bounded cycle remainders. This
    cleanly separates orbit complexity from future-response complexity.
30. Optimizing the finite reward quotient is NP-complete even for identity
    dynamics, three target states, error `1/2`, and rewards in `{-1,0,1}`:
    the optimum is exactly graph chromatic number. A fixed quotient remains
    easy to audit by its cycle LP. Construction and verification are therefore
    distinct resources.
31. A coherent path-lift relation gives a strong all-word theorem but is not
    necessary for scalar spectral response. Exact wordwise equality may use a
    different critical raw witness for each word while every subscale path
    lift retains all microscopic states. The remaining positive problem is a
    synchronization theorem that legitimately interchanges these quantifiers.
32. For `{0,-C}` switching kernels, scalar zero-response equivalence has an
    exact relational form: every word-composed zero relation must contain a
    directed cycle. Failure gives a finite word whose repetition has linear
    drift. This periodic-completeness condition is strictly weaker than the
    left-total relations required by a coherent path lift.
33. The fixed-binary de Bruijn separation persists after blocking until every
    generator is projectively rank one with one critical node. Thus even
    contraction coefficient zero and wordwise unique optimization do not
    synchronize a reusable witness. Positive converses must control how
    exposed images vary *across words*, not only mixing within each product.
34. Scalar tropical response admits a strictly weaker exact path certificate:
    a finite carrier of nonempty endpoint subsets stable under backward
    good-edge lifting. It yields an all-word quantitative spectral bound;
    at zero defect it is equivalent to nonmortality of the subset automaton.
    The canonical carrier has at most `2^r-1` states, and failure supplies a
    bounded-length mortal word whose repetition pumps linear drift.
35. The multi-state extension uses finite endpoint-support families over a
    coarse system. Backward-surjective near-optimal edges realize critical
    paths, while a support potential turns transient shortfalls into a
    cycle-mean error. This gives a generator-checkable quantitative all-word
    theorem strictly weaker than rowwise path lifting.
36. A genuine width-two Ising alphabet needs two anticipatory support states
    but four forward path states. After an interacting weight perturbation,
    the two-state carrier exactly computes the order-sensitive response
    `2N_ca^cyc-N_a`; equal-Parikh words rule out one state.

The classical ingredients are max-plus dynamic programming, Myhill--Nerode
residuals, McShane envelopes, and tropical distance projection.  The
project-level generative content is their resource-complete synthesis:
private compiler plus distance bridge plus sensitivity converse gives an
exact response class, exact distortion, and an interacting composition law;
the mean-field roof adds a distinct strict quotient whose complexity changes
with the future interaction rather than with microscopic size alone.  The
arithmetic feature-algebra law is the first general theorem here to predict a
state-growth exponent from response generators and an exact update algebra.
The experiments do not justify reconnecting to the motivating signing
problem yet.

## Candidate object: the upper response roof

For a finite landscape `H:Omega -> R` and a declared feature map
`phi:Omega -> R^d`, define

```math
\widehat H_\phi(u)=
\max\left\{\sum_x\lambda_xH(x):
\lambda\in\Delta(\Omega),\ \sum_x\lambda_x\phi(x)=u\right\}.
```

This is the upper boundary of

```math
\operatorname{conv}\{(\phi(x),H(x)):x\in\Omega\}.
```

It is the minimal exact quotient for all linear response queries

```math
V_H(\theta)=\max_x\{H(x)+\langle\theta,\phi(x)\rangle\}.
```

The roof has an exact algebra:

- additive composition becomes sup-convolution;
- lifted response bodies add by Minkowski sum;
- every bi-affine cross energy and bi-affine parent feature is recovered
  exactly from the two child roofs, with a non-amplifying query-distortion
  inequality;
- fixed-dimensional bounded response bodies form a compact Hausdorff space
  with unrestricted finite recovery sequences.

This is useful only when the declared feature algebra closes under the next
operation.  If the next coupling asks for an omitted block, root, or
correlation, a larger roof is required.  If `phi(x)=x` on the Boolean cube,
every vertex is exposed and

```math
\widehat H_\phi(x)=H(x).
```

Thus exact sufficiency for arbitrary spin-pinning queries necessarily retains
the complete landscape.  A viable compression theorem must either restrict
the query family, exploit special structure, or prove a synchronization law
that makes omitted features functions of retained ones.

The positive side is not empty: for fixed-rank Curie--Weiss/Potts-type
mean-field Hamiltonians, total feature vectors form a fixed-dimensional
bi-affine algebra.  Their roofs have polynomially many attainable feature
values and give an exact, bracket-independent ground-state dynamic program.
This is the first application in which the candidate state is both exact and
genuinely smaller than the landscape.

The name “upper response roof” is descriptive, not a claim that a new
mathematical field has already been established.  It is ordinary convex
duality used as a query-relative extremal sufficient statistic.

## What the obstruction atlas says

Across the motivating repository and the new examples, failed summaries fall
into three mechanisms.

| Mechanism | What survives | What is lost |
|---|---|---|
| Bulk compression | spectra, fixed moments, bounded-temperature pressure, or positive-density states | an exposed zero-entropy face or planted Boolean resonance |
| Unrooted compression | energy histograms and global overlap laws | which coordinate block, code root, or interface carries the extreme response |
| Nonclosed compression | child scalar channels or a fixed feature set | correlations created by the next coupling; enlarging until closure may recover the full landscape |

The detailed, evidence-labeled atlas is in
[`drafts/obstruction_atlas_report.md`](drafts/obstruction_atlas_report.md).
No single scalar statistic currently survives all three mechanisms.

## Research-director judgment

The selected architecture now has two independent costs.

1. **Algebraic reuse and compatibility:** an approximate state must either be
   a congruence/idempotent retraction or control its reuse defect.  Locally
   valid gauges must also be glued across newly created relations.
2. **Response width and exposure:** omitted information matters only to the
   extent that declared queries separate its response fibres.  Metric entropy,
   posterior width, and rooted packings measure that price.

The framework is **Level 3 for binary extension/Cayley word landscapes and
Level 2 globally**.  In that class the relation space, its holonomy, the exact
gluing dimension, a robust defect-amplification bound, and both positive and
negative response-rate theorems belong to one mechanism.  Across roofs,
matroids, convex bodies, spin-glass overlaps, and codes, the same two costs are
visible but no single theorem yet contains all of them.

The central law is therefore narrower and more precise than a universal
“extremal information” slogan:

> Local summaries compose without new information exactly when their gauges
> glue.  Mixed relations are the obstruction; their holonomy is the new
> feature, and query exposure decides how much of it must be retained.

The unrestricted syndrome dichotomy is no longer the selected target: its
positive side is resolved, although the optimal macroscopic rate remains
open.  The mixed-holonomy test is also resolved in a nontrivial linear regime.
For `kappa <= D/32`, two individually shear-trivial fragments over
`F_q^D` can create at least

```math
q^{3D\kappa/16}
```

kernel-endpoint response profiles separated by `D/16`.  Thus fixed
macroscopic accuracy can require `Theta(D kappa log q)` bits, matching the
exact gluing dimension up to constants, even though the queries carry no
channel labels.

The continuation now supplies the general law inside the class actually
generated by this composition.  A mixed holonomy produces a **presented
carrier** `C` in the endpoint metric, with response

```math
F(x)=\min_{c\in C}\{d(x,c)+\pi(c)\}.
```

If `0<=pi<=p`, its response metric differs from carrier Hausdorff distance by
at most `p`, sharply.  Under diffuse queries, witness-neighborhood mass gives
the corresponding exact exposure charge.  This transfers the lower theorem
to Lee, flag-ultrametric, and rank-metric Cayley realizations, while recording
that the rank proof currently uses an equilateral host.

The same object predicts collapse.  Surjective maps, redundant alphabets, and
a linear-diameter two-scale metric disprove any law based only on `D kappa`.
A metric quotient with fibre diameter `a`, lift defect `b`, and presentation
radius `p` decodes every carrier from its projected carrier to error
`a+b+p`; the error is nonamplifying under fixed min-plus continuation.  This
is strict compression only when projected carrier entropy is smaller and the
composition descends to the quotient.

The framework is therefore **Level 3 for presented-carrier/min-plus
landscapes and Level 2 globally**.  The selected next theorem is no longer a
finite-field packing.  It is to predict, from local fragment data, whether
composition-created carriers are Hausdorff-rich or synchronize through a
smaller metric quotient.  That is the missing step from a capacity diagnostic
to a predictive composition theory.

The first predictive step is now proved for linear metric carriers.  The
largest dimension `s_W(Delta)` of a linear host with minimum distance above
`Delta` gives the lower response exponent, while any synchronizing quotient
of dimension `r` gives the upper exponent.  They obey

```math
s_W(a)\le r
```

whenever quotient fibres have diameter at most `a`.  This single inequality
recovers Hamming and rank-metric Singleton bounds.  A Gabidulin host converts
the rank-metric version into `Omega(D^2 k log q)` response information, so the
non-Hamming validation now uses genuine rank geometry rather than only an
equilateral subspace.  The optimal synchronization rank is exactly the
codimension of a largest linear anticode, yielding the general
code--anticode inequality.  Binary Hamming space then gives a decisive
negative result: at relative scale `delta`, sphere packing forces a linear
duality gap of at least `H_2(delta/2)-delta`.

The Hamming Grassmannian audit now resolves the unrestricted version of that
target negatively.  At one channel, carrier packing is the classical
unrestricted binary coding number up to one codeword; an actual metric cover
also proves that puncturing is exponentially nonminimal.  Sparse-flat counts
inside quotient leader balls give an exact formula for directed
Grassmannian balls, but isometric quotient norms can still have a linear gap
in rooted symmetric response.  At general channel rank, a low-word sum-code
plus injection-distance construction is valid yet provably cannot beat the
common-host Gilbert exponent.  These facts force rooted extension geometry,
not another scalar rank.

The positive synthesis is a new exact product algebra.  For a finite local
response alphabet, retain only

```math
r(a,b)=\sup_x(f_a(x)-f_b(x)).
```

After arbitrary direct-product composition, the uniform distance between
two word responses is exactly the maximum of the two oriented sums of this
table.  Two-sided local carrier gap `d` and presentation radius `p` therefore
give a per-differing-channel margin `d-p`; matching channels pay nothing.
Outer codes amplify this to a macroscopic response rate whenever `d>p`.
Binary simplex lines and rank-metric `F_8` multiplication lines provide two
independent nontrivial validations.

This moves the framework to **Level 3 for product-composed presented
responses and structured tropical continuations**, while it remains Level 2
globally.  The benchmark campaign now separates a reusable state into a
small realizable response image and a future-semigroup congruence.  Unit-load
Max-Cut shows that regularity alone need not make the image small;
shared-parameter presentation bounds its coarse entropy.  Exact tropical
lumpability and metric shells supply congruences; the transition toll shows
that one-step approximation does not.

The recognition/stability target is now substantially resolved and corrected.
Approximate idempotence can drift, but a *shared exact* idempotent relation can
also absorb a coherent defect forever.  Thus no converse may treat fixed
kernel families as fresh adversarial noise.  On a fixed factorial tie-free
selector language with fresh arbitrary residuals, endpoint gauges plus
syndetic tangent resets are complete, and twisted cycle means give the exact
recurrent criterion. The remaining unifying target is a robust converse for
finite selector-affine dynamics: failure of a finite oriented pullback
certificate should yield either a realizable response cocycle or a
quantitative packing of finite-horizon response trees. A static response
cover alone cannot decide this, and syntactically nonzero holonomy may be
dynamically infeasible.

The fifth investigation began the bridge hierarchy and forced a further
revision. A rank-`r` bilinear bridge has an exact associative upper-roof
algebra, and fixed-error response compression costs exponentially in `r` in
the worst case. Exact and approximate complexity sharply separate: rank one
can expose every microscopic atom at exponentially small margins, whereas a
fixed response scale admits a finite-dimensional net. Conversely, neither
sparsity nor bounded degree is a compression promise: a degree-one matching
bridge can expose exponentially many macroscopic response bits for arbitrary
internal landscapes.

Algebraic rank is not the whole positive story. For dense full-rank bridges
`alpha I+beta J`, nonnegative identity channels admit one common nested
optimizer section, while signed cycle balance is the exact gauge condition
for retaining that section. This yields a polynomial magnetization quotient
and a genuine thermodynamic limit for the resulting restricted model; an
unbalanced cycle gives an extensive holonomy loss against separately
optimized pair responses. At the dense `n^2` scale, labeled cut-norm
replacement supplies a different all-future compression theorem, including
rare pinned optimizers. Its generic regularity complexity becomes
exponential at the motivating `n^(3/2)` scale, so it is a calibrated adjacent
theory rather than a solution to the signing problem.

Finally, deterministic de Bruijn dynamics prove that anticipatory-support
size is not itself semantic information. One example has one semantic state
and a one-state forward lift but needs `q^m` exact anticipatory-support states,
with optimal certificate toll `Theta(C/(1+log_q N))`. Thus forward path
simulation and backward anticipatory support are incomparable proof
architectures. Hidden witness memory becomes part of the reusable state only
when the declared future queries can observe it.

The sixth investigation validated and sharpened these claims on harder
benchmarks. Future equivalence independently recovered the Viterbi survivor
vector, the minimal partial-syndrome trellis quotient, the clipped Potts
cavity message, and the discounted-control value vector. Their shared law is
now explicit: a fixed compact polyhedral response image of dimension `d` has
small-error codebook size `Theta(epsilon^(-d))`; indefinite reuse separately
requires an invariant congruence or contraction. Discounting changes the
cover scale by exactly `lambda^h` and gives the sharp repeated-error factor
`1/(1-lambda)`.

The bridge ladder is also calibrated on both sides. A live vertex cover of
size `k` gives a universal `2^k`-entry bridge table, and a matching proves
that exponent sharp. A common multitype partition yields a compact joint
optimization grid under a common-section condition, but not a serial state
when a past microscopic alignment is frozen. At the opposite endpoint, one
dense sign bridge with operator norm `O(sqrt n)` carries exponentially many
independent response bits at `n^(3/2)` accuracy for unrestricted children;
weighted linear children already retain a linear information rate. Thus an
eventual positive theory for dense quadratic signings must exploit rigidity
specific to that child class.

Literal rank is now replaced by a scale-sensitive parameter. Uniformly
truncating every singular value below `epsilon sqrt(n)` changes a balanced
Boolean bridge response by at most `epsilon n^(3/2)`. The surviving numerical
rank feeds the exact upper-roof algebra, with
`exp(O(r_epsilon log(1/epsilon)))` local feature cells under bounded
feature-visible ports. This explains in one law why fixed-rank and
`alpha I+beta J` bridges compress at the target scale, bounded-degree bridges
are subleading there, and typical dense sign bridges do not. It is a local
factor theorem rather than a global quotient across large separators.

That endpoint now reaches the full syntactic class of complete sign
quadratics. Gauge-ferromagnetic pole forms yield an `exp(Omega(n))` response
packing and force `Omega(n)` bits at the same scale, while Hamming covering
and simultaneous discrepancy rounding give `O(n^2)`-bit universal upper
states. A stronger Bernoulli-thinning theorem discards a fixed
`Theta(epsilon^2)` fraction of the coefficient bits into a sparse weighted
surrogate while preserving every Boolean energy to target accuracy. The
lower witness has `Theta(n^2)` cap, so it does not determine the bounded-cap
rate or near-minimizer rigidity.

The bounded-cap endpoint is now quantitatively sharp for a full switching
orbit. A general amplification theorem says that a positive entropy deficit
in the target-scale near-top set forces a positive contextual response rate
under one low-operator-norm sign bridge. On `n=2^(2m)`, Hanson--Wright applies
this to the exact-cap regular-Walsh child. Its full orbit therefore has
matching `Theta(n)` response bits: the lower bound is an exponential
projective packing, while the switch itself is an `n`-bit exact state. Thus
bounded extremal scale alone does not imply sublinear memory.

An intermediate fully algebraic benchmark is also available. Permutation-
valued Maiorana--McFarland switches and the deterministic Walsh bridge give
a matching `Theta(sqrt(n)log n)` response state. Its proof identifies
high-bias approximate self-isometries of the Boolean inner product and shows,
by a Fourier-rank tail, that they are sparse at speed `sqrt(n)log n`.

Repeated composition supplies the next, different obstruction. The
structured Walsh family has an exact `k sqrt(n)`-bit Kronecker presentation
on a `k`-block graph, but isolated bias and all within-word pairwise
truth-table correlations are not reusable. A linear-label commutation bit
between the child and bridge involutions changes the path optimum by
`((3-sqrt(5))/2+o(1))n^(3/2)` per block. This is composition-created
information even though the edge connection itself is flat.

The linear-label subfamily now has a positive exact endpoint. In every label
dimension, binary Gram data, the full relation kernel, and one
characteristic-root relation coset classify every ordered tuple up to Walsh
coordinate symmetry. This gives an `O(k^2)`-bit exact Boolean-extremal
quotient for every `k`-block graph, independent of the `sqrt(n)`-bit ambient
label dimension. It is static rather than independently composable: gluing
creates cross-Gram entries and mixed relations.

That exact gluing fibre is now characterized. It consists of a descended
cross bilinear form, an intersection correspondence between the two
presented spans, and the combined characteristic-root fibre. Pullback of
these three objects reconstructs the joined state and is associative on the
accumulated span. The cross form can carry `rs` bits and the intersection
correspondence `r^2+O(1)` bits with isolated states fixed; pairwise data can
still miss a ternary relation. The next issue is whether an approximate
quotient reduces this relational memory, or whether a nonlinear Walsh-label
family has an analogous closed presentation without performing the complete
Boolean maximization.

There is also an exact query-relative collapse, now at the level of the whole
landscape rather than only its spectrum. Embed a label as `(0,a)` in the full
`2m`-dimensional Walsh coordinate space. The ambient characteristic vector
never lies in the label span, so Gram plus relations extend to one ambient
orthogonal coordinate permutation conjugating every child and common bridge.
Thus every unrooted real weighted graph landscape factors through `(G,R)`;
the label-space root fibre appears only after a pole or field pins the old
coordinate splitting. The earlier Weyl-word theorem is the corresponding
spectral shadow.

Declared query locality reduces the remaining gluing state further. On each
maximal connected support one retains only the restricted cross form and
coincidence relation (plus a root fibre only for root-sensitive components).
The cost is `O(sum_C|C|^2)=O(wL)`, hence linear for bounded-size,
bounded-incidence support families. This rate is sharp: independent
three-block path queries expose one coincidence bit apiece at a fixed
`n^(3/2)` gap. A separate ordinary triangle exposes an off-diagonal
Gram/flux bit, even after self-parities, relations, and root fibres are fixed.
What remains open is the metric rate of these fluxes on a large connected
family, not whether they can be scalar-visible at all.

Approximate locality has a matching scale law. Deleting cross-component
interaction of total weight `delta` costs at most `delta n^(3/2)` and leaves
the exact local orbit states. A unit path therefore has an `O(t/eta)`-bit
carrier at error `eta t n^(3/2)`, whereas this total-variation deletion
architecture stays quadratic on a dense graph with edge weights bounded
away from zero.

The sparse upper argument also survives abstraction. For any finite public
bounded feature dictionary, a linear-size list of importance-weighted masks
gives an explicit response code controlled by row variance and `log|X|`.
This one theorem covers Littlewood/CSP dictionaries and code correlations as
well as quadratics. It is a one-shot semantic replacement theorem, not yet a
compositionally invariant quotient.

An orthogonal benchmark adds a different kind of compression. In a
boundary-case branching random walk, critical renormalization turns the
finite-depth pair `(W,Z)` into one limiting derivative mass `Z`; conditional
on it, the unmarked extremal process is a decorated Cox process and composes
by a smoothing transform. This is a limiting distributional state, not an
exact finite-port quotient, and genealogy-marked futures force a larger mass
measure. It is retained as a scoped rare-event branch.

The same branch now has a deterministic leading-rate compactness theorem.
Upper-semicontinuous microcanonical log-count hypographs compose by supremal
convolution and have recovery sequences; finite count convolution inherits
the limit when decomposition and descriptor complexity are subexponential at
the declared speed. A bounded-temperature pressure family can miss an
isolated maximum of positive normalized height, while the hypograph retains
it. This is a speed-sensitive response roof, not a full extremal process:
subexponential multiplicity and spacing remain invisible.

A second orthogonal benchmark keeps disorder discrete and adversarial. For a
finite-state nearest-neighbor chain, minimum transfer-product norms are
submultiplicative, so fixed-temperature pressure is a lower spectral radius;
a uniform soft-max sandwich then proves convergence of the adversarial
ground-state density. Strictly positive transfers have a contractive Hilbert-
metric cavity and polynomial finite nets for asymptotic mean pressure. This
does not transfer locally to a dense sign split: every dense sign bridge has
Boolean bilinear response at least `n^(3/2)/sqrt(3)`, already the leading
scale, and the ordinary transfer interface has `2^n` states.

## Extremal rate--distortion

For a landscape class `H`, query set `Theta`, and response metric

```math
d_\Theta(H,G)=
\sup_{\theta\in\Theta}|V_H(\theta)-V_G(\theta)|,
```

the deterministic information price of uniform error `epsilon` is, up to the
usual factor-two gap, the metric entropy of the response class:

```math
\operatorname{Pack}(\mathcal H,d_\Theta,2\epsilon)
\le K_\epsilon
\le\operatorname{Cov}(\mathcal H,d_\Theta,\epsilon).
```

Under an ensemble, the corresponding quantity is Shannon's rate--distortion
function for response loss.  The query family is part of the definition.
Encoding one maximum, one optimizer, all bounded fields, and all possible
pair couplings are different information tasks and can have scalar, linear,
or quadratic rates.

This formulation connects convex response geometry, approximate sufficient
statistics, sketching, and information theory without identifying them.  The
finite quadratic lower bound in [`theorems.md`](theorems.md) is the first
nontrivial example.

## Canonical model tests

| Model | Landscape | Natural extreme scale | Interface that must be declared |
|---|---|---|---|
| Dense quadratic signs | `sum_(i<j) a_ij x_i x_j` | `n^(3/2)` for optimized signs | fields, restrictions, or bridge features |
| SK model | variance-normalized random quadratic form | order `n` | overlap/external-field experiment |
| Dense Max-Cut/CSP | satisfied-constraint count | order `n^2` before centering | vertex prizes and added constraints |
| Code/coset distance | negative distance or membership reward | block length | root, puncture, and amalgamation queries |
| Random CSP ground state | centered clause Hamiltonian | model-dependent extensive scale | clause additions and overlap constraints |
| Littlewood polynomial | Boolean evaluation of a unimodular polynomial | degree/support dependent | restrictions and product substitutions |

An invariant is suspect if it only restates one of these models.  The current
results already apply independently to quadratic forms, vertex-prize Max-Cut,
and code covering-radius queries.

## Evidence discipline

The surface files use four labels:

- **Proved**: a complete finite or asymptotic proof is in `theorems.md` or a
  cited repository artifact and has been independently audited.
- **Exact computation**: exhaustive integer enumeration with a reproducible
  script, but no extrapolation beyond the checked range.
- **Imported**: a published or repository theorem whose hypotheses are stated
  and mapped explicitly.
- **Open**: a falsifiable target, not a conjectural fact.

Machine-learning prediction is not evidence.  The experiments are intended
to find collisions and invariants that can be translated back into exact
mathematics.

## Files

- [`examples.md`](examples.md): compact obstruction atlas and canonical
  counterexamples.
- [`axioms.md`](axioms.md): definitions and principles that survived the
  first falsification pass.
- [`theorems.md`](theorems.md): rigorous results only.
- [`open_questions.md`](open_questions.md): minimal next mathematical targets
  and their stopping tests.
- [`experiments/`](experiments/): exact finite programs, outputs, and a small
  order-eight landscape laboratory.
- [`drafts/`](drafts/): detailed specialist reports retained for audit, not
  promoted automatically as theory.

## Literature coordinates

The framework is grounded in, but not reducible to, several mature theories:

- Shannon's [rate--distortion theory](https://ieeexplore.ieee.org/document/5311476);
- Blackwell's [comparison of statistical experiments](https://doi.org/10.1214/aoms/1177729032)
  and Le Cam's [approximate sufficiency](https://doi.org/10.1214/aoms/1177700372);
- Kolmogorov--Tikhomirov
  [metric entropy](https://www.mathnet.ru/php/archive.phtml?jrnid=rm&option_lang=eng&paperid=7289&wshow=paper);
- Derrida's [REM](https://doi.org/10.1103/PhysRevLett.45.79) and
  [GREM](https://doi.org/10.1051/jphyslet:01985004609040100);
- the zero-temperature Parisi formula of
  [Auffinger--Chen](https://arxiv.org/abs/1606.05335); and
- Panchenko's [ultrametricity](https://doi.org/10.4007/annals.2013.177.1.8)
  and multi-species
  [synchronization](https://arxiv.org/abs/1310.6679); and
- Madaule's [decorated extremal process for boundary-case branching random
  walks](https://arxiv.org/abs/1107.2543);
- Birkhoff's [projective contraction theorem](https://doi.org/10.1090/S0002-9947-1957-0087058-6);
  and
- Guglielmi--Zennaro's [lower-spectral-radius and antinorm
  review](https://arts.units.it/retrieve/handle/11368/2972807/338571/GZ2019_finale.pdf).

The research question is which compression mechanisms behind these theories
remain valid for deterministic, adversarial landscapes without importing
random-model identities that fail there.

A finite tower of separated exponential scales now also has an exact
lexicographic count algebra on fixed descriptor sets. Its sharp failure on
Vandermonde convolution shows why subleading extremal compactness must retain
saddle multiplicity, not only pointwise rate coefficients.

## Reconnection rule

Do not force this framework back onto the original `M_n` problem.  Reconnect
only after it has all of the following:

1. a nontrivial general information or composition theorem;
2. an application outside dense sign matrices;
3. a state demonstrably smaller than the full Boolean response landscape for
   the required bridge queries; and
4. a constrained compactness/realization theorem for that state.

Items 1--3 now exist for restricted models: fixed-rank mean field, bounded
separators, fixed-label syndrome fragments, and fixed-effective-rank additive
response sets. A constrained version of item 4 now exists for exact
regular-Hadamard tensor amplification: its fixed-port response sets are
nested and yield normalized limits for genuine dense hollow sign hierarchies.
This is only a geometric-order structured benchmark; it neither realizes
near-minimizers at arbitrary orders nor closes the general dense-sign bridge
interface, so the campaign still does not return to the motivating problem.
