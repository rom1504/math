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
problems instead of extending its vocabulary.  It produced five conclusions.

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
   compression promise; polynomial component size or bounded description is
   a genuinely stronger unresolved resource.
4. Finite-metric distance kernels form a nonproduct interacting algebra:
   strengths combine by a bottleneck minimum and isometry labels by holonomy.
   Anisotropic projective-Hamming shells strengthen this to a transported
   coordinatewise-minimum lattice.
5. The same distance shell is the exact nearest Lipschitz response to an
   arbitrary profile.  Repeated shells pay precisely the weakest-layer
   Lipschitz defect, not one approximation loss per layer.

The classical ingredients are max-plus dynamic programming, Myhill--Nerode
residuals, McShane envelopes, and tropical distance projection.  The
project-level generative content is their resource-complete synthesis:
private compiler plus distance bridge plus sensitivity converse gives an
exact response class, exact distortion, and an interacting composition law.
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
responses and metric-shell continuations**.  The directed table now survives
a genuine cross-block minimization on the finite-metric family, with robust
sublinear cumulative loss whenever the total entrywise perturbation is
sublinear.  The anisotropic shell state is one load vector plus a monomial
holonomy, demonstrably smaller than a general transfer table.

The result is deliberately narrow.  It transports arbitrary profiles but
does not compress them, and an exponential private compiler can hide behind
unit interface load.  The strongest next theorem is therefore a recognition
and stability result: characterize when a structured continuation is close
to the metric-shell semigroup using less information than its full kernel,
or prove a scalable obstruction showing why approximate tropical
idempotence alone cannot provide such a quotient.  In parallel, determine
which global resource on separator components yields nontrivial response
compression beyond raw description counting.

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
  [synchronization](https://arxiv.org/abs/1310.6679).

The research question is which compression mechanisms behind these theories
remain valid for deterministic, adversarial landscapes without importing
random-model identities that fail there.

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
response sets.  Item 4 exists only for unrestricted finite landscapes at
fixed interface dimension.  None of these states yet closes the dense-sign
bridge interface, so this campaign does not reconnect to the motivating
problem.
