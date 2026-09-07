# Three precise principles from the signing campaigns

2026-09-07. Synthesis-track working proof map, not a convergence claim.
Original problem:

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad Q(A)=\max_x|H_A(x)|,\quad
M_n=\min_A Q(A),
```

where A is hollow, symmetric, and has a sign on every off-diagonal entry.
The verified interval is

```math
0.4333221116640807\le\liminf M_n/n^{3/2}
\le\limsup M_n/n^{3/2}<0.493608094.
```

The new director selector/stability construction improves the old exact
upper `U=0.49360809358874865...` by a fixed strictly positive amount.
Its improvement is symbolic, not a newly certified useful decimal digit;
see `principle_director_selector_entropy_stability_2026_09_07.md`.

The user allows selectable optimal children and arbitrary changes to old
internal edges. No claim below silently strengthens the needed quantifier to
every representative or to every additive near-minimizer.

## Readback and provenance

The main synthesis inputs were README, ACTIVE_STATE, latest ledger 10.152--
10.157, `flatify_director_final_synthesis_2026_09_07.md`, the earlier
`retrieval_grounded_panel_synthesis.md`, `ar_director_synthesis.md`, and
`extremal_information/CAMPAIGN_SYNTHESIS_266d101.md`. The last is the preserved
ten-hour theory/bridge synthesis and is the plausible earlier Sol synthesis
mentioned by the user; the file itself does not identify the model author, so
that attribution is not asserted as verified.

Original proofs were read for the local-profile collision, exact insertion
collision, tensor classification, Cayley cap identity, fourth-moment bound,
residual-width lower theorem, projector localization, quenched Gaussian-sign
transfer, regular-core/refill entropy theorem, and direct-E upper construction.
The 2026-09-07 111937Z preserved diagnostic rubric was also read. Its demand
for every required exact minimizer is retained only when required by the
specific construction, not promoted to a universal user requirement.

## Principle 1: preserve the joint query until its optimization is paid

This is an operation-specific rule, not a demand that every successful proof
compress the entire future or retain the whole spin landscape.

For a fixed-child block extension with bridge R, the exact identity is

```math
Q\!\begin{pmatrix}A&R\\R^T&D\end{pmatrix}
=\max_{x,y}\bigl(|H_A(x)+H_D(y)|+|x^TRy|\bigr).            \tag{P1}
```

Reversing one block changes the bridge sign and leaves internal quadratic
energies unchanged. Thus favorable bridge cancellation cannot be assumed at
a fixed pair; both bridge polarities occur among legal parent spins. The
quantifier in (P1) is EVERY pair x,y for ONE chosen bridge R. Taking separate
maxima is an upper bound, not an equivalent joint objective.

Two consequences explain genuinely different archived failures:

1. **Insufficient state under a specified query.** Complete magnetization
   extrema of two order-seven EXACT minimizers agree, but their minimum
   one-vertex extension caps are 12 and 10. Thus that statistic is not closed
   under arbitrary signed-row insertion. The counterexample does not rule out
   choosing the favorable exact representative or changing old edges.
2. **Overpayment after forgetting joint geometry.** For EVERY bounded-cap
   actual child and EVERY rank-o(n) subspace, its near-orthogonal Boolean slice
   has half-width at least `(c_R-o(1))n^(3/2)`,
   `c_R=(2/3)sqrt(2/(3pi))`. Therefore a separately paid residual-norm parent
   certificate costs at least `(1+2c_R-o(1))n^(3/2)>1.614 n^(3/2)`, above the
   balanced target `2sqrt(2)M_n<1.397 n^(3/2)`. This lower-bounds the
   certificate, not any actual optimized bridge.

The ten-hour contextual-response synthesis identifies the same phenomenon
in exact finite examples: a one-step defect table need not remain a
congruence after tensor composition; individual product deficits need not
control a joint selector. Its positive algebraic closure theorems succeed
because the declared channels close jointly, with a composable relative
error. These examples are not near-minimizer theorems and do not force a
small response state on original optimizers.

The new Gaussian-sign comparison is a genuine success under this principle:
it preserves all deterministic child offsets and compares the quenched
maximum uniformly over the ENTIRE configuration set. Its scalar one-pair
predecessor could not pay the entropy of the central shells. However, the
favorable covariance-matched Gaussian parent bound remains unpaid.

The strict-upper weave is a second success. It does not separately pay each
child cap. The exact full-spin sign/permutation counting, graph Cauchy--
Schwarz, and information supersolution control one whole nonlinear
construction. It may erase a prescribed seed and still prove an upper bound:
seed preservation is not a premise of that upper theorem.

There is now an actual constant-multiplier test of this principle. In the
two-spin rank-one weave, preserve `S` exactly and set
`C_ij=S_ij eta_ji`, `D_ij=S_ij eta_ij eta_ji`. Independent port signs create
an iid opposite child `D`. After selecting a high-energy `D` state, its
optimized bridge field contributes `sqrt(2/pi)n^(3/2)`, while the old seed
energy has a favorable exact conditional expectation. A finite-polynomial
random-child witness gives a parent coefficient at least
`(2/pi+sqrt(2/pi))/(2sqrt(2))=0.5071738708...` with high probability.
Uniform balancing of each incidence row changes the parent by
`O_K(n^(5/4)log^(5/4)n)` in cap, so this also rejects independent EXACT
two-row orthogonal physical frames for every fixed bounded-cap seed,
including any selected exact minimizer. It is not a theorem about all
correlated completions. The joint-response proof and exact checks are in
`principle_synthesis_2026_09_07_constant_two_frame_obstruction.md`.

The construction agent's independently audited extension covers EVERY
fixed physical multiplier `k>=2` for the explicit independent Walsh-type
frame law. Independent macro-block Lindeberg comparison costs only
`O_k(n^(4/3))` in expected cap; a Gaussian opposite-child witness and
greedy bridge give
`c_k=[(2/pi)(k-1)^(3/2)+sqrt(2/pi)sqrt(k-1)]/k^(3/2)`.
This is increasing in `k`, with minimum `c_2>.50717`. An exact
category-balancing coupling has `o(sqrt(n))` operator cost and transfers
the obstruction to the stated independent orthogonal frames. This is
not a theorem about arbitrary orthogonal frame laws or growing `k`;
see `principle_construct_2026_09_07_fixed_frame_gaussian_obstruction.md`.

A NEW actual-sign operation now supplies a third success. Anchor an actual
negative maximizing spin y and flip only edges agreeing with y_i y_j. After
t flips its negative extreme is EXACTLY R+2t, regardless of the chosen edge
set. A uniform random prefix has mean `(1-p)A-p offdiag(yy^T)`, whose positive
extreme never exceeds the original Q. Uniform scalar and matrix concentration
are paid before selecting the first crossing of the two extrema. This produces
an actual same-order signing with

```math
|P(B)-R(B)|\le2,\quad
Q(B)\le Q(A)+O(\sqrt{n|P(A)-R(A)|}+n),\quad
\|B\|_{\rm op}\le\|A\|_{\rm op}+O_C(\sqrt n)
```

for every input with `Q(A)<=Cn^(3/2)`. The retained state is only the one
protected witness and the exact monotone negative response; the unprotected
positive maximum is uniformly bounded, not reconstructed or compressed.
See `principle_synthesis_2026_09_07_global_balancing.md` for the finite theorem.

A further positive joint-query theorem now treats the actual balanced
bulk's variance ratio as a single ferromagnetically completed objective.
Its exact spin-flip condition includes the `2b/q` hollow-diagonal
correction. Conditioning on arbitrary heavy incoming signs, a convex
minority-violation statistic has positive mean and controlled Lipschitz
constant. A common threshold `A rbar` then gives row penalties
`exp[-c q r_i/rbar]`, whose graph product is at the full physical scale.
Consequently sparse local maxima must carry almost all minority mass
on coherent incoming-tail fibres; a mean-weighted version permits the
other mass only near the balanced face. This holds for every balanced
transform array with bounded operator norm under actual independent
outer signs, not for Gaussian replacement fields. An explicit Walsh
example attains variance coefficient `1/2` and shows why its coherent
exception cannot be dropped from a balance-only theorem. See
`principle_synthesis_2026_09_07_balanced_bulk_sparse_stability.md`.

The coherent exception has a further positive partial resolution.
In the actual stratified and exactly balanced ensemble, uniform column
coherence is `mu=O_L(sqrt(m log m))`. Signed reciprocal-port factorization
and a restricted Gram bound give EVERY word on `ell>=2` nonconstant
fibres the relative bound
`|H|<=[q+(ell-2)mu]sum_i||P_i x_i||^2/2`, simultaneously for EVERY
outer seed. Thus `ell=o(sqrt(m/log m))` permits arbitrary bias scales
and coherent spectra at coefficient `sqrt(p)/2+o(1)<1/2`. This does not
pay the many-active mixed phase; see
`principle_synthesis_2026_09_07_balanced_bulk_sparse_active_bound.md`.

**Operational conclusion.** A new summary is useful if it closes the actual
queries made by one constructor, or if its loss is separately paid at the
original scale. Failure of a particular summary is not evidence that all
convergence mechanisms require a large state. No unpaid replacement of
`max(f+g)` by separately controlled f and g is allowed.

## Principle 2: law and objective scale matter more than an isolated statistic

Low cap gives a genuine uniform theorem,

```math
Q(A)\le Cn^{3/2}\quad\Longrightarrow\quad
\operatorname{tr}A^4\le16K_G^2 C^2n^3.                   \tag{P2}
```

It does not give bounded `||A||op/sqrt(n)`. Changing a clique of order
`n^(3/4)/log n` inside an exact minimizer costs `o(n^(3/2))` in cap while
making the normalized operator norm diverge. Consequently dividing a whole
covariance by the product of the two GLOBAL operator norms can erase its
leading effect on additive near-minimizers. The proved iid-like failure is
about this particular opposite-spectral law's expectation and typical outputs.
It is not a theorem about all covariances, rare selection, selectable exact
minimizers, or global old-edge changes.

Conversely, a repaired regular core and random refill give ANY low-cap input
a same-order full-sign replacement with cap cost `2sqrt(epsilon)n^(3/2)` and
operator norm `O_C(epsilon^-1)sqrt(n)`. The adaptive mask costs only O(n)
disorder-entropy bits, so it preserves the leading n^2-speed low-cap counting
entropy. This revives fixed-accuracy bounded-operator analysis despite
arbitrarily bad global spectral normalization. It does not yield one fixed
operator cutoff with a vanishing cap error and does not by itself give a
cross-order recurrence.

There is a NEW positive law theorem making the statistic/law distinction
quantitative. For ANY prescribed actual signing A, random relabelling and
vertex switching give a law supported entirely on its cap-preserving orbit.
The full edge covariance is exactly the identity, while for B=A+I its
k-vertex law L_(A,k) satisfies

```math
d_{\rm TV}(L_{A,k},U_k)
\le\sqrt{\frac{(2^{k-1}-k)\operatorname{tr}B^4}{2n^4}}
      +\frac{k(k-1)}{2n}.                               \tag{P3}
```

Thus (P2) constructs exact-cap-preserving laws that are uniformly locally iid
through `(1-delta)log_2n` vertices, for EVERY bounded-cap input, including EVERY
exact minimizer. Yet an orbit law of one fixed seed has support at most
`n^k2^k`, so its k-vertex law has total variation tending to one from iid once
`k>(2+delta)log_2n`. The logarithmic threshold is sharp in order, not in its
constant. Proof and exact finite checks:
`principle_synthesis_2026_09_07_logarithmic_local_law.md`.

This is an actual structured sign-law construction, not a Gaussian
certificate. It explains why covariance matching or even growing local
agreement cannot by themselves universalize optimizing disorder. Loop-space
dependence survives all the erased edge moments. It does NOT contradict
Gaussian-sign universality, whose specified input-law assumptions exclude
arbitrary gauge laws. Nor does it assert that these loop variables must be
compressed to prove convergence.

The next distinction is sharper: retaining genuine cycle information
still need not retain it in a PARTICULAR leading pressure. The new
mixed signed-input/ordinary-output weave works for every fixed macro
seed without any final column sign gauge. Its `K_5` moment distinguishes
seed switching classes even after conditioning on the good-cap event.
Yet its complete independent-port Gaussian-kernel pressure is leading
seed-independent for ALL row profiles of squared energy at most `Cm`,
including coherent spikes and arbitrary asymmetry. Descending signed
matching, category entropy, tagged-zero swaps, and a profile-chosen
low-energy annulus prove this without a tail-integrability hypothesis.

The director's aggregate residual-energy argument sharpens the seed
scope to `Q=o(m^2/log m)`; a matching-scale explicit Hadamard-block
counterexample shows that merely `Q=o(m^2)` is insufficient. The same
compiler proves minimum-port-defect universality and a uniform Laplace
principle in ENERGY units for every diverging temperature. Near-optimal
port placements have probability at least `exp(-O_epsilon(m^2))`, rather
than paying for an entire permutation description. None of these results
interchanges query-dependent placements with one parent's simultaneous
Boolean maximum. Full proofs are in
`principle_synthesis_2026_09_07_complete_seed_pressure_universality.md`
and `principle_director_sharp_port_universality_2026_09_07.md`.

An immediate original-input consequence improves the archived random-
restriction falsifier: uniform `k=(1-delta)log_2n` restrictions of EVERY
bounded-cap signing have normalized cap at least
`(2/3)sqrt(2/pi)-o(1)=0.5319230405...-o(1)` with high probability. Cap is
switching-invariant, so the switched-local theorem applies to ordinary
restriction caps. Rare favorable restrictions and different/global operations
are not excluded.

**Operational conclusion.** Any proposed universality theorem must specify
its law, conditioning, and uniformity. A statistic can be useful without
determining the objective, and an unfavorable law can be replaced when global
changes are permitted. Exact-minimizer statements cannot be inferred from
planted near-minimizer witnesses.

## Principle 3: distinguish orbit persistence from accuracy-first realization

Literal tensoring preserves a very specific algebra. For EVERY FIXED full
symmetric sign seed B of order d>=2, even self-tensor powers have normalized
cap tending to one half if B is Hadamard and to infinity otherwise. In the
Hadamard branch `vec(B)` supplies a Boolean endpoint eigenvector at the
second tensor power. In the non-Hadamard branch an explicit strict vector
value grows exponentially and even powers make all arcsine terms nonnegative.
This is a complete theorem for fixed-seed even self-tensoring, not for a
growing fresh seed, one-time tensoring, non-Kronecker fusion, or arbitrary
full-sign constructions.

The elementary-abelian Cayley obstruction is different and stronger in its
native class. For EVERY independently chosen Boolean function g at EVERY
native order, all characters are Boolean eigenvectors, so

```math
Q(A_g)=\frac n2\max_a|\widehat g(a)-g(0)|
\ge\frac n2(\sqrt n-1).                                \tag{P4}
```

It defeats scale-dependent functions within that algebra, not merely a fixed
tensor seed. It says nothing about arbitrary signings at those orders. A
selected arithmetic family cannot prove a universal high subsequence for M_n.

The actual strict-subhalf construction escapes these persistence hypotheses
and has a paid limiting schedule:

1. Fix a target cap margin and choose a finite recursive depth r through the
   conditional-variance supersolution and Gaussian-boundary stopping theorem.
2. At this FIXED r, the all-spin word/pair counting and Fock orbital estimate
   are uniform in EVERY terminal orthogonal matrix and every terminal vector,
   including large coordinates. Hence terminal choices may vary with order.
3. H2/H12 tensor orders supply ratio-dense terminal sizes. The fixed-depth
   construction and principal deletion then cover all sufficiently large
   original orders with vanishing loss. Only last is the target margin removed.

At p=24/25 and t=97/20 the exact ternary interval bound supplies

```math
U=\frac{97/20+(24/25)\log2-5151/6250}
        {(97/10)\sqrt{24/25}}
=0.49360809358874865\ldots .                             \tag{P5}
```

Tensor terminal orders do not mean that the nonlinear recursion depth grows:
the terminal-uniform theorem is what permits these distinct roles. Neither
the fixed-seed tensor no-go nor the Cayley Boolean-eigencharacter floor
applies to the resulting restricted nonlinear weave as a general operation.

**Operational conclusion.** Accuracy-first finite complexity plus uniform
all-order realization can prove an upper theorem without proving optimum
preservation. Convergence still needs an actual liminf-to-all-order mechanism;
for example a summable selected-child defect. Algebraic/tensor phase effects
prove nonconvergence only if one separately obtains an all-signings lower
obstruction on infinitely many other orders.

## Quantifier audit in one table

| Result | Exact covered class | Not established |
|---|---|---|
| Six-vertex profile collision | Two scalable bounded-cap actual families | Near-minimizer collision; selected optimal representative impossibility |
| Magnetization insertion collision | Two exact minimizers at order seven | Asymptotic obstruction; failure of selecting the better child |
| Even self-tensor classification | Every fixed full symmetric seed, literal even powers | Growing fresh seeds; arbitrary global fusion |
| Additive Cayley half-floor | Every scale-dependent g in the stated elementary-abelian algebra | Universal lower floor for arbitrary signings |
| Fourth moment and logarithmic local law | Every bounded-cap actual signing, hence all exact/near-minimizers | Objective determined by local law; bounded normalized operator norm |
| Residual envelope floor | Every bounded-cap child, every sublinear-rank projection; that envelope | Actual parent cap floor for a selected bridge |
| Opposite-spectral bad law | A specified normalized law on planted additive near-minimizers | All exact minima; every law; rare output; global changes |
| Quenched Gaussian-sign transfer | Every input of its specified bounded-spectrum elliptic Gaussian-sign law, arbitrary offsets/configurations | Favorable Gaussian parent bound; arbitrary covariance-matched laws |
| Strict subhalf upper | Existence of actual full signs at every large order | Optimality of ensemble; convergence of original minima |
| New paired-row anti-invariant weave | Existence at every half-order of the exact skew-cross family, normalized cap below .499 | Seed-sensitive completion; a better original upper endpoint |
| Fixed-seed mixed-orbit weave | Every fixed macro seed; true cycle-sensitive output laws; old strict upper | Leading dependence of the cap on the seed value |
| Complete port pressure and ground theorem | All per-row energy-bounded profiles; independent uniform ports; seed envelope `Q=o(m^2/log m)` | Actual parent maxima; seed-dependent profile laws; unscaled growing-temperature pressure |
| Fixed-multiplier independent frames | Every fixed bounded-cap seed, every fixed k>=2, then the explicit independent exact Walsh-category frames; typical actual cap above .50717 | Arbitrary orthogonal laws; correlated or rare selected completions; growing k; one simultaneous event for adaptive seeds |
| Sparse balanced-bulk stability | Every conditionally fixed balanced transform array with bounded operator norm, actual independent outer signs; all local maxima of both ferromagnetic completions | Whole mixed-profile cap; arbitrary seed-conditioned signs; controlling coherent tails by balance alone |
| Sparse-active balanced bulk | Actual stratified/repaired columns, every seed/mask simultaneously, every word on o(sqrt(m/log m)) active fibres with arbitrary bias scales | Many-active mixed phase; complete variance ratio; selected-child transfer |

## Positive result and remaining obligation

The first new theorem here is (P3), a cap-preserving actual-law construction with a
quantitative logarithmic local-iid window and an order-sharp support boundary.
It strengthens a universal restriction theorem on ACTUAL exact minimizers and
removes any inference that favorable signs require nontrivial unconditional
edge covariance. It is not itself a convergence recurrence.

There is now also a genuinely new global repair, proved in
`principle_synthesis_2026_09_07_global_balancing.md`: every bounded-cap input
can be balanced to polarity gap <=2 at same-order cap cost O(n^(5/4)) while
increasing operator norm by only O_C(sqrt n). The earlier clique repair gave
gap O(n) and could create n^(3/4)-scale operator spikes. The new distributed
path therefore yields simultaneous fixed-accuracy spectrally regular,
essentially exactly balanced near-minimizers at every order. Its edge count
also preserves n^2-speed low-cap disorder entropy. It does NOT identify the
width minimum with the absolute-cap minimum or supply a cross-order defect.

A third positive theorem now makes the certificate-versus-object distinction
constructive. In `principle_synthesis_2026_09_07_anti_invariant_weave.md`, a
paired-row, two-class weave gives actual matrices of the exact form
`[[A,C],[-C,-A]]`, `C` skew, with normalized cap below .498795641 at every
large half-order. Thus the family's spectral-certificate half-floor is NOT
an actual Boolean half-floor. The same fixed tournament is imposed before
independent row sampling; product-orbit Finner keeps the two child channels,
and density convexity reduces the root continuum to twenty certified scalar
points. A direct off-diagonal defect identity avoids an invalid leading-cost
triangle deletion. The completed directed certificate checked 100,466 boxes.
The construction is seed-blind in the same outer-sign sense as the original
weave, so it removes a family obstruction without supplying seed landing.

The later fixed-seed mixed-orbit theorem makes seed retention genuinely
nontrivial at the law level, but the complete port theorem identifies its
exact leading-pressure boundary rather than declaring all cycle information
useful for optimization. Its finite-accuracy entropy consequence is positive:
one can approximately realize the best response to each prescribed profile
without specifying all row permutations at `m^2 log m` cost. The remaining
issue is simultaneous-query realization, not lack of an approximate
single-query compiler. Conversely, the independent constant-two-frame law
keeps the original energy scale and still fails through an actual joint
Boolean witness; simply reducing the multiplier does not repair independence.

Finally the sparse-stability theorem and the independently proved
stratified packet theorem are actual positive components of the mixed
bulk problem. The former handles arbitrarily heterogeneous minority
scales but leaves quantified coherent/near-balanced mass. The latter
handles the coherent packet geometry uniformly within bounded-ratio
near-constant bands, reaching coefficient `sqrt(p)/2` in the one-hole
case. Uniform column coherence additionally handles every bias scale
and coherent spectrum when few fibres are nonconstant, at the same
matching-scale relative coefficient. Neither their union nor separate endpoint bounds pay the full
mixed cross interaction. The remaining multiscale phase obligation is
therefore explicit, not silently replaced by the original scalar
orbital certificate that was already shown to fail near constant bias.

The surviving useful revival is to design conditional, non-Gaussian, or global
operations while retaining their joint query, rather than insist on the
failed global spectral normalization. The older regular-core/refill theorem
already pays fixed-accuracy spectral repair; it is not claimed as new here.

For convergence the present selected-child sufficient target remains:
for comparable N=m+l, select actual minimizers A_m,A_l and construct a full
hollow signing C_N, allowing all old edges to change, such that

```math
Q(C_N)\le Q\!\left(\operatorname{diag}\left(
\sqrt{\frac{N-1}{m-1}}A_m,
\sqrt{\frac{N-1}{l-1}}A_l\right)\right)
+O(N^{3/2-\delta})
```

for a fixed delta>0. The common principles explain which proven failures do
and do not touch this target. They do not replace it by a covariance-only,
local-state, fixed-seed, or separately maximized certificate.
