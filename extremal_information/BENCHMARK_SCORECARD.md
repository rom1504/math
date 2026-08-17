# Benchmark scorecard

Compact scorecard for the campaign begun at `266d101`.  “Independent” means
the operational state was frozen before comparison with the named classical
representation.  Detailed proofs and exact checks are in `drafts/` and
`experiments/`; this file does not promote an unaudited draft by itself.

| Benchmark | State predicted from future responses | Classical or comparison state | Independent? | Exact / approximate complexity | Verdict |
|---|---|---|---:|---|---|
| Viterbi / finite-state best path | endpoint survivor vector modulo one additive scalar | Viterbi survivor metric; max-plus transfer vector | yes | exact `q-1` real coordinates; `Theta((D/epsilon)^(q-1))` projective codebook; global zero-temperature contraction is `0` only at max-plus rank one and otherwise `1` | pass; predicts a sharp mixing obstruction and near-rank-one reset |
| Binary linear-code trellis | reachable past modulo the past-supported subcode; equivalently the partial syndrome in `im H_P intersect im H_F` | minimal trellis quotient `C/(C_P direct-sum C_F)` | yes | `2^r` labels, `r=dim C-dim C_P-dim C_F`; rich weighted tables need `Theta(2^r log(B/epsilon))` absolute bits | pass; strict quotient and raw-interface no-go both predicted |
| Potts tree / finite-state graphical model | baseline plus clipped projective message `C_K={r in [-K,0]^q:max r=0}` | max-sum / junction-tree message | yes | exact projective dimension `q-1`; sharp codebook `Theta_q((K/epsilon)^(q-1))`; an exact finite lattice closes after one-time factor rounding | pass; a new scoped rate and reusable approximation accompany the classical state |
| Discounted deterministic control | full terminal value vector, projectively modulo one scalar | Bellman value function / discounted semigroup | yes | exact depth-`h` metric is `lambda^h` times sup or half-oscillation; bounded-reward rate is `Theta(n log(1+B_H lambda^h/epsilon))`; fresh errors cost `1/(1-lambda)` | pass; contraction produces a sharp forgetting-time information law |
| Finite-width Ising | projective boundary response or smaller rational anticipatory carrier under restricted futures | max-plus transfer matrix | earlier solution-hidden run | exact size depends on declared futures; unrestricted futures expose the whole table | pass; calibrated the query dependence |
| Bounded-separator Max-Cut/CSP | projective boundary response table | treewidth dynamic-programming table | earlier solution-hidden run | pure futures expose the full one-Lipschitz response ball; logarithmic cover number is exponential in width | pass; local sensitivity is not approximate compression |
| Fixed-rank Boolean bridge | concave upper roof over the aggregate feature | finite-rank mean-field / support-function state | derived directly | exact associative roof algebra; fixed-error cells `exp(O(r log(1/epsilon)))`; worst cases need `2^(Omega(r))` bits | pass; exact and approximate complexity sharply separate |
| Scale-sensitive spectral bridge | roof over singular features above the target-scale threshold | truncated SVD plus low-rank factor model | derived and independently audited | tail error `n sigma_(r+1)` is uniform over all futures; bounded feature-visible ports need `exp(O(r_epsilon log(1/epsilon)))` local cells; separator size still controls global memory | pass as a multi-benchmark local law, not a global dense-graph quotient |
| Degree-one matching bridge | bridge response transform on the live interface | separator DP with extensive interface | no classical state supplied | arbitrary children require `exp(Omega(n))` response bits at error `Theta(n)` | negative; bounded degree is not a compression parameter |
| Vertex-cover bridge width `k` | `2^k`-entry conditioned bridge-response table, quotiented further only by numerical envelope collisions | separator / bucket-elimination table | yes | worst-case projective rate `Theta(2^k log(1+D/epsilon))`; matching attains it | pass pending final audit; explains why degree-one matching is hard when its live cover is extensive |
| Structured full-rank dense bridge | one magnetization per common permutation orbit, plus signed cycle gauge | rearrangement / mean-field count state | derived directly | polynomial joint search for fixed orbit count; a restricted thermodynamic limit follows; not an arbitrary frozen-fragment quotient | pass; rank is not the controlling resource |
| Random dense sign bridge | no small all-landscape state: exponentially many exposed response coordinates | no proposed classical compression | theorem independently audited | arbitrary children require `exp(Omega(n))` bits at error `Theta(n^(3/2))`; linear children already require `Omega(n)` bits | strong negative endpoint; special quadratic child classes remain open |
| Complete sign-quadratic children across a dense sign bridge | coefficient signing, with response collisions allowed | quadratic coefficient description | lower and upper developed by opposing agents; lower independently audited | projective packing forces `Omega(n)` bits; a sparse weighted surrogate stores at most `(1-epsilon^2/4)binom(n,2)+O(log n)` bits, still `Theta(n^2)` | pass as an information bracket; this planted lower has `Theta(n^2)` cap and does not explain the bounded-cap rate |
| Exact cap-`1/2` Walsh sign quadratics | coordinate switch for the full orbit; smaller truth table/permutation for explicit subfamilies | regular Hadamard/Walsh switching representation | yes; all three packing theorems independently audited | near-top entropy amplification gives matching `Theta(n)` bits for the full orbit; explicit Boolean-table and permutation families give `Theta(sqrt(n))` and `Theta(sqrt(n)log n)` states with deterministic Walsh bridge | strong near-original pass; bounded cap alone does not collapse response memory |
| Repeated Walsh-bridge graph composition | unrooted `(Gram,relations)` state; rooted relation form; local `(kappa,J[,Z])` amalgamation; Kronecker carrier for arbitrary truth tables | no classical extremal quotient supplied | yes; relative-involution, orbit, gluing, Weyl/ambient-symmetry, and query-local arguments independently audited | exact `O(k^2)` global orbit quotient; the entire unrooted landscape forgets `Z`; bounded-incidence connected supports cost `O(sum_C|C|^2)` and independent paths force one `J` bit each; unrestricted gluing retains `rs` cross-form and `r^2+O(1)` intersection orbit memory | strong pass: first strict Boolean-extremal quotient with an exact nonlocal composition law, a genuine rooted/unrooted symmetry transition, and a sharp local semantic rate; minimal scalar state on one large connected support remains open |
| Finite signed feature dictionary | one public sparse mask plus retained signs | importance sampling / uniform sparsification | generalized only after quadratic theorem | saves `Omega(min{m,E^2m/(V_Phi log|X|)})` bits and protects all shared max-type futures | new multi-model upper law; not an invariant internal quotient |
| Generic random weighted dense Max-Cut at scale `n^2` | density/block representative in cut norm | weak regularity graphon state | literature-grounded | exact table exponential, but one block gives all-future `O(n^(3/2))=o(n^2)` error | pass at dense scale; fails to resolve `n^(3/2)` scale |
| Mean-field BEG microcanonical landscape | usc log-count hypograph over occupation/energy descriptors | large-deviation microcanonical entropy | literature-grounded and independently audited | exact sup-convolution recovery under subexponential decomposition count; bounded-temperature pressure loses a rare maximum; exponential descriptor images and subexponential decorations are sharp falsifiers | scoped pass for leading-rate rare-event compactness |
| Boundary-case branching-random-walk extremes | derivative mass `Z` for unmarked limiting Laplace queries | derivative martingale / decorated Cox process | primary-literature scout, theorem hypotheses checked | one scalar after critical renormalization; finite depth requires `(W,Z)` and marked futures require a mass measure | orthogonal pass; a macroscopic state can emerge only in the limit |
| Finite-width chain with discrete adversarial disorder | lower spectral growth rate; under positivity, normalized projective cavity plus additive potential | transfer matrix / lower spectral radius | literature-grounded and independently audited | pressure and ground-state density converge; a `delta`-net has `O(delta^{-(q-1)})` states and mean-pressure error `delta/(1-kappa)` | orthogonal pass; nonconvex disorder retained, but dense cuts have leading `n^(3/2)` response and exponential transfer width |
| Deterministic de Bruijn carrier | semantic response is scalar; anticipatory certificate retains suffix phase | finite automaton / subset carrier | derived directly | one semantic and one forward state versus `q^m` exact anticipatory states; bounded rooted probes recover `q^m` states only under a gap/filter condition | sharp falsifier and repaired observability theorem |

## Current prediction

Across the independently frozen finite-port benchmarks, the framework
recovers the classical state by the same operation: take the realizable image
of the conditional response vector and quotient only directions invisible to
all declared futures.  Its small-error exponent is the dimension of that
realizable projective image, not raw separator width.  Reuse under arbitrary
depth additionally requires a congruence or invariant arithmetic carrier.

The bridge ladder introduces a different transition.  Low rank and common
optimizer sections produce strict quotients; sparsity without a small live
separator does not; a generic dense sign bridge exposes extensive information
at the same scale as the motivating problem. Complete sign quadratics still
carry an extensive response rate, although they have a strict polynomial-bit
coefficient representation. The full exact cap-`1/2` Walsh switching orbit
also carries a linear response rate. Under repeated Walsh composition, the
smaller bias/pair-overlap state fails by an extensive commutation-holonomy
gap. The remaining question is whether a joint gauge-covariant state can
transport the structured truth table without reconstructing the Boolean
maximum, and whether any such structure is forced for near-minimizers.

The branching benchmark records a second kind of successful state. It is not
an exact finite contextual quotient: critical centering and a vanishing
martingale coordinate create a scalar limiting Cox intensity. This keeps the
rare-event branch mathematically distinct from finite-port response geometry.
