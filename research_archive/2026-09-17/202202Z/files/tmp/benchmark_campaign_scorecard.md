# Final research-director audit: benchmark campaign

**Scope.** I audited commits `660ad18` through `6dec554` and the current
working tree, including the compact theory files, all benchmark proof drafts,
the independent mean-field and carrier-growth audits, and the current
uncommitted Theorems 16.15--16.16. I reran the ten benchmark verifiers. All
advertised finite checks passed. The judgments below rely on the written
proofs; enumeration is treated only as a formula falsifier.

## Executive judgment

The campaign meets the user's severe success criterion in several ways:

- it independently recovered the exact composable state in the separator,
  Ising, and weighted-automaton benchmarks;
- it proved genuinely quantitative approximate-state theorems beyond the
  standard exact algorithms, most strongly for unit-load Max-Cut and
  finite-atom mean field;
- it produced a general arithmetic feature-algebra growth law which predicts
  the polynomial exponent of additive response states; and
- it extended directed response beyond product composition through metric
  shells, then found sharp obstructions and finite certificates for their
  depth-uniform approximation.

The severe qualification is that much of the exact-state layer is classical
dynamic programming, convex duality, or weighted Nerode theory. The framework
has become **coherent and locally generative**, but not a universal theory of
extremal compression. Its strongest content is the joint insistence on two
separate obligations:

```math
\boxed{
\text{usable future state}
=\text{small realizable response image}
+\text{a congruence for the continuation semigroup}.}
```

Neither clause follows from the other. The benchmark counterexamples now
make that a theorem-guided research principle rather than vocabulary.

## Benchmark scorecard

| Problem | State predicted by theory | Known classical state | Was discovery independent? | Exact / approximate minimality result | Verdict |
|---|---|---|---|---|---|
| **Bounded-separator Max-Cut / binary CSP** | The conditional optimum profile on boundary assignments; for pure Max-Cut this is a projective table on `X_w={+-1}^w/{+-}`, with the absolute offset stored separately. Under a boundary-load declaration it becomes the corresponding projective Lipschitz profile. | Treewidth-DP separator table, modulo global spin flip for Max-Cut. | **Qualified yes.** The solution-hidden derivation fixed the response state before archive comparison. The basic conditional table is nevertheless the standard and nearly forced answer. Positive-edge projective lookup, exact unit-load realization, and the presentation-size separation were derived afterward and were not fed as the standard solution. | Arbitrary CSP boundary fields expose every one of `2^w` oriented coordinates. Pure Max-Cut futures expose every one of `2^(w-1)` projective coordinates, and an unrestricted private compiler realizes every table. With unit boundary load, the response class is exactly the full projective one-Lipschitz ball. At error `epsilon w`, its response-bit count lies between `2^((1-H_2(2epsilon)+o(1))w)` and `2^((1-H_2(epsilon)+o(1))w)`. An `m`-edge shared-parameter presentation has only `O_epsilon(m^2+m log(w+m))` response bits, so a universal macroscopic net needs exponentially many edges. | **Strong success (A/B).** This is the best benchmark. Caveats: the exact entropy exponent is not closed; the universal compiler is exponential; nothing comparable is proved for polynomial-size Max-Cut components beyond the presentation bounds. |
| **1D / finite-width Ising** | Width one: baseline plus the cavity gap `d=f(+)-f(-)`; a two-ended fragment: a projective `2 by 2` max-plus kernel. Width `w`: the conditional table on `2^w` exposed spin assignments. | Zero-temperature transfer matrix / cavity field and the usual width-`w` boundary table. | **Yes, but low-surprise.** A solution-hidden agent derived the gap recurrence, two-ended kernel, and general boundary profile before repository comparison. These are classical objects that the operational definition almost forces. | Exact pinning proves coarseness, max-plus multiplication proves repeated composition, and the width-one clipped recurrence is explicit. A bounded projective `2 by 2` kernel body has cover size `Theta((B/epsilon)^3)`. For gadget-complete pairwise width-`w` fragments the table requires `Theta(2^w log(B/epsilon))` bits. Signed chains further close on a baseline, sign holonomy, and bottleneck magnitude. | **Clean sanity-test success (A).** The framework earns credit for deriving the state for the right contextual reason. It does **not** yet give a sharp approximate lower bound for a strict nearest-neighbour finite-width strip or for bounded-degree translation-invariant strips; the full-cube lower construction uses exponentially many private gadgets in a general tree decomposition. |
| **Heterogeneous / fixed-rank mean field** | For anonymous union plus one chemical potential: the discrete-concave top-occupancy profile, equivalently its sorted slope multiset or hinge-response function. On a finite grid: the field histogram. With fixed bilinear pair score: the least concave roof; at sufficient positive curvature this synchronizes to total field alone. | Empirical field distribution / sorted local fields; magnetization or color counts in homogeneous Curie--Weiss/Potts; concave-envelope/Legendre response under linear fields. | **Not cleanly solution-hidden.** Fixed-rank aggregate closure was already present in the repository and the benchmark was explicitly prescribed. The heterogeneous slope state, exact roof congruence, sharp collapse threshold, and rate theorem were agent-derived during this campaign, but should not be advertised as blind rediscovery of a classical abstraction. | Exact biconjugacy proves coarseness and an exact response-metric identity. Sorted multiset union is associative. Common-grid histograms add exactly and incur at most `eta N/2` root error on any merge tree, using `O((1+B/eta)log N)` bits; choosing `eta=B/sqrt(N)` gives sublinear bits and error. On `d` fixed atom types the exact number of states is `Theta(n^r_Z)`, and for equally spaced hinge atoms the exact cover below half-grid error is `{n+d-1 choose d-1}`. The bilinear roof is a strict congruence; `J>=4B/n` sharply collapses every mass-`n` state to total field. A matching macroscopic rate-distortion law for arbitrary continuous fields is not proved. | **Strong success (A/C), with a rate gap.** This is the clearest evidence that the framework predicts different growth regimes rather than merely cataloguing boundary tables. The claims require anonymous contexts, one uniform field, and fixed same-`J` normalization; separately labelled blocks, changing `J/n`, or non-biaffine interactions can resurrect discarded information. |
| **Weighted automata / weighted languages** | The restricted future residual `rho_x:y mapsto L(xy)`, not the raw forward vector. Under tropical block lumpability, the strict state is the vector of nonlinear block maxima `P_a=max_(i in I_a)(p_i+c_i)`; for affine suffixes it is the visible upper hull. | Weighted Myhill--Nerode/Hankel residual, max-plus transfer vector, equitable/tropical lumping, support-function hull. | **Qualified yes.** The residual and explicit four-to-two lumping were derived solution-hidden before comparison. The collision with classical weighted Nerode and tropical lumpability is exact and is correctly acknowledged. | Residual equality is the coarsest exact contextual equivalence, derivatives compose at every depth, and the lumped aggregate updates exactly. Query nets give an external-cover upper bound; robust coordinate pins give `Theta(r log(1+B/epsilon))` bits on an exposed `r`-aggregate box. A one-dimensional suffix family exposing all `p` raw coordinates falsifies affine dimension as a compression certificate. For one fixed automaton, full minimality still requires reachability of the aggregate box; without it only the reachable residual set is controlled. | **Successful harder validation (A/B), mostly classical exact algebra.** The response-entropy/exposure sandwich and low-dimensional-query falsifier are the useful additions. A general sharp residual rate-distortion theorem remains open. |
| **Directed response under interacting continuation** | For metric-isometry kernels: bottleneck strength plus isometry holonomy; anisotropic projective-Hamming shells retain a transported load vector. More generally, one must retain exposed witnesses and a depth-stable quotient congruence, not merely the static directed table. | Tropical distance transforms, transfer kernels, gauge/cocycle cancellation, finite-state semigroup actions. | **Not a solution-hiding benchmark.** This was the agent-authored live theorem target. Classical metric projection and graph cohomology were used openly; the exact combination, falsifiers, and perturbation-quantifier split were developed in the campaign. | Metric shells compose exactly under genuine min-plus interaction: strengths take a minimum and isometries compose. Directed row gaps and weakest-layer distortion are exact. Compatible gauges telescope; small recurrent images reset old error. Nearby exact actions of one finite semigroup have bounded-normal-form stability. Against fresh residuals on a fixed selector language, syndetic tangent resets are necessary and sufficient up to `r(r-1)`, twisted cycle means are the exact recurrent obstruction, and the finite suffix-product lift returns either a uniform bound or a pumpable drift cycle. | **Strong specialized success (D/E), not the generic target.** It extends beyond products, but only for structured shell/selector semigroups. Entrywise perturbations still cost `sum eta_i` unless an additional gauge, reset, or algebraic relation is proved. Switching perturbations that change active cells remain unsolved. |

## New general theorems that survive severe audit

### 1. Arithmetic feature-algebra growth

For a fixed finite set of additive response atoms, the exact mass-`n`
contextual quotient has

```math
Theta(n^{r_Z})
```

states, where `r_Z` is the integer rank of the atom-response differences.
Real conditioning gives coarse two-sided response covers, while lattice
margin gives exact microscopic covers. This is the first result in the
campaign that actually **predicts** a polynomial state-growth exponent from
the response generators and the update algebra.

The result is elementary finitely generated abelian-monoid geometry, so it
should be called generative within this project, not a field-level novelty.
Its conceptual contribution is the necessary separation of exact arithmetic
rank, robust real conditioning, and lattice-scale distinguishability.

### 2. Exposure and presentation control different information prices

Landmark nets upper-bound response entropy, while balanced robust exposure
lower-bounds it. For the full Lipschitz language the exposure number is
exactly an interface packing number. Independently, shared-parameter
max-affine presentations have an `O(m^2)` coarse upper bound and a robust
`Omega(m log m)` facet lower bound. These theorems explain why semantic
profile dimension and syntactic grammar size are not interchangeable.

### 3. Genuine interacting composition

Metric-shell kernels form a closed nonproduct algebra with exact bottleneck,
holonomy, and weakest-layer distortion. The anisotropic shell is a
coordinatewise-minimum lattice. This supplies the required example in which
the directed response is preserved because cancellation occurs inside the
joint continuation, before separate absolute-value bounds are paid.

### 4. Depth stability depends on the perturbation quantifier

The campaign now has four rigorously separated mechanisms/results:

1. endpoint gauges and zero recurrent holonomy telescope;
2. recurrent small images reset prior error;
3. coherent exact semigroup relations absorb a fixed defect through bounded
   normal forms; and
4. against fresh selector residuals, absence of syndetic tangent resets
   forces adversarial linear growth.

The finite suffix-product lift makes the fourth alternative decidable for a
regular tie-free selector language. Twisted functional-cycle means replace
ordinary untransported holonomy on one recurrent selector cell.

## Failures and counterexamples that forced principled revision

1. **Local sensitivity is not compression.** Unit-load Max-Cut realizes the
   entire projective Lipschitz ball, with exponentially many response bits at
   linear distortion.
2. **Low query-parameter dimension is not compression.** One affine line of
   weighted-automaton suffixes can robustly expose every raw coordinate.
3. **One-step approximate congruence is not reusable.** A transition-toll
   family has vanishing local idempotence defect but fixed long-depth drift.
4. **Global tropical contraction has no useful middle regime.** On the full
   projective domain an all-finite max-plus linear map has Hilbert coefficient
   zero or one. Weak links help through bounded image diameter, not a generic
   coefficient below one.
5. **Kernel gauge plus small full-image reset is too narrow for coherent
   maps.** Nearby idempotent clamps stay close forever by a shared exact
   semigroup relation despite nonzero rectangle defect and no small image.
   This does not refute broadly defined zero-increment recurrence; it forces
   the coherent/adversarial distinction.
6. **Arithmetic exactness is not robust exposure.** Rationally independent
   atoms may have high exact integer rank while their real conditioning and
   lattice margin vanish.
7. **A visible roof is context-dependent.** Changing the Curie--Weiss
   normalization or allowing non-biaffine cross interactions can re-expose
   points discarded by the fixed-`J` concave roof.
8. **Abstract table realizability can hide exponential syntax.** The
   Max-Cut/CSP lookup compilers validate semantic lower bounds but do not
   establish the same response class for polynomial-size or strict-lattice
   instances.

These failures improve the theory because each kills a proposed implication,
not merely one construction.

## Overclaim and missing-obligation audit

1. **Max-Cut:** the macroscopic cover exponent retains the
   `H_2(epsilon)` versus `H_2(2epsilon)` gap. The compiler proving the lower
   semantic class has exponential private size. Do not infer a polynomial-
   size lower theorem from unit boundary load alone.
2. **Finite-width Ising:** the full profile-cube lower bound is for general
   treewidth-`w` pairwise fragments with exponential private gadgets. A
   strict nearest-neighbour strip may have a smaller realizable profile
   class; its sharp approximate response complexity was not determined.
3. **Mean field:** the continuous-field macroscopic lower bound does not
   match the grid-quantization upper bound. The roof theorem is only for the
   declared anonymous, same-`J`, linear-field continuation semigroup.
4. **Weighted automata:** robust pins prove minimality on a formal aggregate
   box. A fixed automaton needs a reachability/product-packing hypothesis.
5. **Directed response:** the exact shell theorem is highly structured. No
   theorem yet preserves a generic directed table through arbitrary weak
   cross-block couplings with sublinear cumulative loss. The statement
   `sum eta_i=o(leading scale)` is an extra hypothesis, not a consequence of
   entrywise closeness.
6. **Finite-semigroup absorption:** it is a useful finite certificate but is
   close to quotient compatibility and zero-increment recurrence. It should
   not be advertised as disproving every broad gauge--reset decomposition.
7. **Feature-algebra theorem:** the equally spaced specialization should
   explicitly assume `d>=2`; its projective constants differ by a factor two
   from the anchored sup-norm constants. The current working tree otherwise
   incorporates the independent audit's fixed-mass and external-cover
   qualifications.
8. **Selector decision theorem:** the current finite suffix-product proof is
   correct for a declared start set and a fixed tie-free selector language.
   Its one-vertex exhaustive verifier does not by itself test arbitrary
   graph start/co-reachability conventions, which should be stated when the
   theorem is finalized.
9. All finite experiments passed, including the newly added roof-congruence,
   exact grid-state, and suffix-product checks. They do not replace the
   hypotheses in the analytic statements.

## Did the framework predict something, or was it fed the answer?

Both occurred.

- The exact conditional profile, transfer matrix, and weighted residual are
  classical and nearly immediate once the future query experiment is
  declared. Their independent rediscovery is a valid sanity check, not
  strong generative evidence.
- The framework did predict several facts not supplied as benchmark answers:
  the pure-Max-Cut projective lookup compiler, the exact unit-load Lipschitz
  class, the presentation-versus-semantic entropy separation, the
  heterogeneous mean-field slope/roof quotient with a sharp synchronization
  threshold, the arithmetic feature-growth exponent, and the precise split
  between coherent and adversarial depth stability.

The latter results justify calling the program generative in these model
classes. They do not yet amount to a universal representation theorem.

## Is the theory coherent?

Yes, at **Level 2 globally and Level 3 in several bounded-interface or
finite-feature classes**. Three objects now have stable roles:

1. **Contextual response equivalence** identifies the semantic quotient.
2. **Response metric entropy / exposure / presentation complexity** prices
   the information retained at a declared accuracy.
3. **Feature-algebra congruence** decides whether the same quotient survives
   future composition.

The additive atom theorem explains polynomial growth; separator lookup
explains exponential semantic width; weighted lumpability explains a strict
fixed carrier; metric shells and selector theorems explain depth reuse. This
is one operational theory rather than a taxonomy, although no theorem yet
derives all three layers from general hypotheses.

The correct recommendation is to **deepen this theory**, not branch to a new
one and not reconnect to the signing problem.

## Single strongest next theorem

### Paired-selector skew-product theorem with switching cells

Let two finite families of piecewise-affine projective max-plus
continuations be locally `epsilon`-close. Build the finite reachable automaton
of paired active selector cells, including ties through a declared finite
refinement. After factoring exact common semigroup relations, prove that on
every pumpable strongly connected component exactly one of the following
finite certificates occurs:

1. the affine residual is a twisted endpoint coboundary on every surviving
   selector cycle, and every remaining residual direction meets a rank-one
   selector product within a bounded gap, yielding a uniform `C epsilon`
   projective error; or
2. a pumpable paired selector cycle has unequal transported cycle means,
   yielding an explicit word and admissible residuals with
   `Omega(T epsilon)` drift.

The fixed-itinerary suffix-product theorem proves this when perturbations do
not change the active cell. The missing theorem must handle selector
switching and ties without replacing coherent fixed kernels by fresh
adversarial noise.

This is the strongest next theorem because it would link the static response
quotient to depth-stable reuse, unify the Ising, weighted-automaton, shell,
gauge, reset, and clamp examples, and produce either a finite verifier or a
finite falsifier. A mere `O(T epsilon)` estimate does not count.

## Why there should be no reconnection to `M_n` yet

None of the successful benchmarks supplies a state for the signing landscape
which is both demonstrably smaller than its full Boolean response and closed
under cross-order composition.

- The separator theorem is cautionary: when interface width grows, even unit
  local sensitivity can retain exponentially many future responses.
- Mean-field compression relies on a fixed finite atom algebra and anonymous
  aggregate coupling, neither of which has been derived for arbitrary near-
  optimal sign matrices.
- Weighted lumpability assumes an exact block congruence not known for the
  signing problem.
- Metric shells form a special tropical semigroup, and no verified mapping
  sends signing composition into that semigroup without a leading loss.
- The paired-selector theorem itself is not yet proved for switching cells.

Therefore the framework has earned further development on its own terms, but
it has **not** earned a cautious reconnection to `M_n`, much less a new attack
on convergence.
