# Minimal open questions

These are theorem targets selected after the carrier-capacity checkpoint, not
a catalogue of neighboring ideas.  The unrestricted syndrome dichotomy is no
longer the lead question: a composition-stable subexponential response net and
a quadratic-bit macroscopic packing are both proved.  Optimizing that gap is
secondary to understanding the new carrier law.

## 1. Find a small cycle-faithful quotient without partition enumeration

The exposed-carrier dichotomy is now sharp in general. One polynomial-size
continuous rational lattice-PWA selector circuit plus one identity probe has
`2^(r-o(r))` exposed recurrent states below per-step error `1/2`. A globally
affine selector and a strict irreducible max-plus-linear map remain
Landau-subexponential. Thus compactness, rationality, nonexpansiveness, a fixed
alphabet, and polynomial presentation do not guarantee compression; the
min/max switching class can already encode essentially every Boolean phase.

For a *proposed* finite input congruence, Theorem 17.1l computes its optimal
asymptotic reward error by a raw-cycle linear program. But feasible
bounded-error congruences are not closed under join, so ordinary monotone
Nerode/bisimulation refinement cannot discover a canonical coarsest one.
The finite projective-semigroup carrier is exact but output-sensitive.

Theorem 17.1n now gives a nontrivial structural certificate: a dominating
coarse system plus a coherent **path-lifting relation** preserves all aligned
word responses, even when the microscopic representative must switch along
the path. Pairwise commuting max-plus generators collapse further to a common
eigenprofile. In the opposite direction, finding the smallest general
cycle-compatible congruence is NP-complete even for identity or rank-one
reset dynamics.

The de Bruijn family now falsifies the most obvious response-to-structure
converse: after blocking `2m` letters, every product is projectively rank one,
has contraction coefficient zero and a uniquely exposed critical node, yet
every sub-`m`-defect block path lift retains all `2^m` states. Uniform
scrambling and wordwise exposed uniqueness do not synchronize the critical
seed across different future words.

The scalar `{0,-C}` case now has a complete answer, and Theorem 17.1s extends
it to arbitrary weighted coarse systems: finite endpoint-support families,
backward-surjective lifts, and a support potential give quantitative all-word
spectral control. The width-two Ising example proves this is strictly smaller
than forward transfer state and genuinely needs state-dependent tolls.

The support-certificate converse is now closed at its natural level. Theorem
18.5 identifies its optimal toll with a finite mean-payoff support game. But
the deterministic de Bruijn shift proves that the resulting carrier size is
not semantic response information: scalar response and rowwise state both
have size one while exact target-surjective support needs `q^m` states. Its
certificate tradeoff is `Theta(C/(1+log_qN))`.

Witness observability now has a first positive theorem. A backward support,
forward envelope, and saturated bounded orbit readout turn phase profiles
into an exact rooted-response packing under a uniform gap. Cyclic rotations
and one local symbol readout expose all deterministic de Bruijn phases, so
the theorem does not smuggle in an arbitrary lookup table. Finite leakage
also gives a sharp negative law: every common future remains nonexpansive in
the endpoint sup norm and cannot amplify a gap beyond `C`.

Bounded-delay weighted synchronization now supplies a second positive
structure theorem.  If every `D`-block has max-plus row rank one, normalized
suffix residuals form a finite exact response quotient, residual-context
cycles give exact accumulated toll, and a descending core decides whether
one locally thresholded support can be attached to each residual.  The
empty-core counterexample proves that this last condition is proof memory,
not automatically semantic response information.

The remaining target is **approximate residual/core synchronization**:

```math
\text{how do profile-covering error, cycle toll, and support leakage combine}
\quad\text{without paying the same approximation at every depth?} \tag{OQ.1}
```

Seek a theorem for products that become `epsilon`-row-rank-one after a bounded
delay.  It should return a finite residual cover plus a cycle potential and
either a depth-uniform rooted response bound or a pumpable pair-cycle
falsifier.  If supports are also claimed, charge the additional greatest-core
leakage separately rather than identifying it with residual error.

**Success:** a depth-uniform quantitative theorem in two model classes, with
profile covering number, forgetting delay, and terminal query radius all
explicit.

**Stop condition:** do not infer a scalar response gap from an empty support
core.  The exact two-matrix counterexample in Theorem 17.1u is mandatory.

The one-generator PWA question is closed: the binary/block counter proves
`2^(r-o(r))` exposed classes. Do not spend another campaign optimizing that
base unless it yields a positive quotient theorem elsewhere.

## 2. Prove a same-input recoupling theorem at linear channel rank

For mixed-relation models, composition creates a homomorphism

```math
h:K=Z/Z_{\rm loc}\longrightarrow W
```

and a presented carrier `C_h=h(K)`.  Theorem 12.3 reduces response complexity
to carrier metric entropy.  Theorem 13.3 identifies the optimal synchronizing
quotient rank with anticode codimension, while Theorem 13.4 proves that binary
Hamming separated rank lies a linear amount below it.  Universal scale-rank
duality is therefore false.

The actual Grassmannian packing is

```math
P_{D,k}(\Delta)=
\log_2\operatorname{Pack}
\bigl(\operatorname{Gr}_k(F_2^D),d_H,\Delta\bigr).
```

The unrestricted form of this question has now been reduced to a recognized
hard object.  At `k=1`, Theorem 14.3 proves

```math
A_2(D,\Delta+1)-1\le P_{D,1}(\Delta)\le A_2(D,\Delta+1),
```

so its packing exponent is the classical unrestricted binary coding rate.
An explicit probabilistic line cover, using the exact line-ball volume,
proves separately that puncturing retains exponentially too many states.
Thus neither a universal closed exponent nor anticode-quotient optimality is
a credible next target.

At `k=kappa D`, a systematic carrier has generator `[I_k|X]`, and

```math
d_H(C_X,C_Y)\le
\max_u wt(u(X-Y))
\le d_{\rm column}(X,Y).
```

The outer column alphabet has size `2^k`, large enough for Reed--Solomon
codes to meet Singleton.  Ordinary column coding therefore reproduces the
puncturing exponent and cannot narrow the gap.  The next theorem must charge
the **same-input recoupling** in the first inequality.

For fixed admissible `kappa,delta`, prove either

```math
\log_2P_{D,\lfloor\kappa D\rfloor}(\lfloor\delta D\rfloor)
\le
\bigl(\kappa(1-\delta-\kappa)-c_{\kappa,\delta}\bigr)D^2
+o(D^2)                                                        \tag{OQ.2}
```

for some explicit `c_(kappa,delta)>0`, or construct a matching
quadratic-bit family showing that no such recoupling gain exists.

**Success:** a positive recoupling exponent in one nontrivial linear regime,
or a construction saturating puncturing after the exact Hausdorff metric is
verified, transferred through the presentation toll to a response
rate--distortion theorem.

**Stop condition:** a bound only on the number of differing columns cannot
work: the growing-alphabet shadow already meets Singleton.  Likewise, the
sparse-flat spectrum exactly counts directed balls but is not two-sided
sufficient; isometric quotient norms can have a linear rooted-response gap.

This remains the sharpest model-specific target.  It asks whether
composition turns the microscopic coupling between a quotient and its rooted
lifts into quadratic-scale response information.

## 3. A natural strict metric synchronization

Theorem 12.5 shows that one uniform approximate submetry
`varpi:X->Y` compresses all presented carriers when four conditions hold:
subscale fibres/lift defect, smaller projected-carrier entropy, descended
composition, and controlled presentation radius.  The designed two-scale
metric satisfies all four; rank-row projection and puncturing presently prove
only the small-error factor.

Find a natural model with growing interface in which all four conditions are
forced by the model rather than installed by hand.

**Success:** a strict quotient whose state has asymptotically smaller metric
entropy than the full carrier class and whose error remains subscale under
unbounded composition depth.

**Stop condition:** an identity map, a projection with no entropy saving, or
a fixed-context approximation without a descended state update is not
compression.

## 4. Query-mass-sensitive carrier capacity

The uniform carrier law is sharp.  Under a query distribution, Theorem 12.3
only gives the local exposure bound

```math
(\Delta-2r-p)_+\,\mu(B(x_0,r))^{1/s}.
```

In Hamming and rank-metric ambient spaces this mass may be exponentially
small.  Determine when a carrier family has a *massive* collection of
separating witnesses rather than isolated Hausdorff witnesses.

**Success:** a geometric or algebraic condition implying an `L^2` or mutual-
information lower bound of the same exponential order as uniform carrier
packing, with an application outside the finite-field Hamming example.

**Falsifier:** a full uniform carrier packing whose normalized distance
transforms all collapse under the proposed diffuse query law.

Do not infer average hardness from one endpoint witness; its neighborhood
mass must be charged explicitly.

## 5. Bridge complexity between rank and regularity

The fixed-rank port is now exact: upper roofs compose associatively, explicit
feature quantization uses `exp(O(r log(1/epsilon)))` cells, and an unrestricted
bounded class needs `2^(Omega(r))` response bits. A matching bridge shows that
bounded edge degree alone can still carry exponentially many macroscopic
bits. Signed-balanced `alpha I+beta J` bridges give a different positive
mechanism: a common optimizer section defeats full algebraic rank.

There is no hidden subextensive SVD regime for bounded-operator dense sign
bridges.  Frobenius mass forces at least
`n(1-epsilon^2)/(C^2-epsilon^2)` singular values above
`epsilon sqrt(n)` when `||R||<=C sqrt(n)`.  Any intermediate positive result
must therefore use nonlinear synchronization or symmetry rather than a
finer spectral cutoff.

The intermediate sparse parameter is now calibrated. A live vertex cover of
size `k` gives a universal `2^k`-entry bridge table, and a matching proves the
worst-case projective rate `Theta(2^k log(1+D/epsilon))`. A common partition
into `K` coordinate types gives an `O(mK log n)` joint optimization grid when
all blocks remain jointly reoptimizable, but counts are not a serial state if
a past alignment is frozen.

At the negative dense endpoint, a random sign bridge with
`||B||_(2->2)=O(sqrt n)` carries `exp(Omega(n))` response bits at
`Theta(n^(3/2))` error for unrestricted children. Weighted linear children
already carry `Omega(n)` projective bits. Thus density, spectral norm, and
the target scaling alone cannot yield an all-landscape quotient.

The negative theorem now reaches genuine complete sign quadratics. Planted
pole forms `((x^Tz)^2-n)/2` give `exp(Omega(n))` response states and hence
`Omega(n)` bits at the same scale. This closes the question “does quadratic
syntax alone compress?”, but not the question that matters near the
motivating model: the planted forms have `Theta(n^2)` cap. A universal
Bernoulli-thinning code discards a fixed `Theta(epsilon^2)` fraction of the
coefficient bits into a sparse weighted surrogate, while an internal Hamming
cover and simultaneous discrepancy rounding give other `O(n^2)` upper
bounds. None is subquadratic; the internal coefficient-ball architecture
provably cannot become so.

The bounded-cap rate question is now answered for one full switching orbit
on an infinite subsequence. For `n=2^(2m)`, a regularized Walsh signing and
its switching orbit lie in

```math
\mathcal Q_C(n)=\{A in \{-1,1\}^{\binom n2}:
                  \max_x|H_A(x)|\le Cn^{3/2}\}              \tag{OQ.5a}
```

with `C=1/2`. Hanson--Wright gives a positive entropy gap in its fixed
`n^(3/2)` near-top set. Theorem 21.8 amplifies that gap into an
`exp(Omega(n))` projective packing through one sign bridge of operator norm
`O(sqrt n)`. Storing the switch gives the matching `O(n)` upper bound. Thus
this full orbit has `Theta(n)` response bits under coordinate-pinning futures,
even though every child has identical isolated cap, spectrum, and all other
switching invariants.

The first **joint gauge-covariant reuse** theorem is now proved for a strict
linear-label family. In every label dimension, the triple

```math
(\text{binary Gram},\text{relation kernel},
  \text{characteristic-root fibre})                      \tag{OQ.5b}
```

classifies every ordered label tuple up to the orthogonal group. It is an
exact `O(k^2)`-bit Boolean-extremal quotient on every `k`-block Walsh graph,
independent of the ambient `sqrt(n)` label bits. Gram alone fails at leading
scale, and even Gram plus relations misses the characteristic root under a
rooted future.

The exact **gluing fibre** is now known. It is the triple consisting of the
cross bilinear form, the coincidence correspondence between the presented
spans, and their combined characteristic-root fibre. It reconstructs the
joined state associatively. Fixed isolated states can carry `rs` arbitrary
cross-form bits or `r^2+O(1)` intersection bits, and pairwise amalgamation
data can miss a ternary relation. At least one such coincidence bit changes a
three-block scalar maximum at full `n^(3/2)` scale.

The first query-local reduction is now exact. On maximal connected supports,
only restricted cross forms and coincidence correspondences are needed for
unrooted queries; root fibres survive only where a continuation pins an
external pole or field. The carrier costs `O(sum_C|C|^2)`, and independent
three-block paths prove a matching linear lower bound for bounded-size,
bounded-incidence supports.

The ambient symmetry question is also closed: `(Gram,relations)` conjugates
the **entire unrooted weighted Walsh graph landscape**, not merely its
spectrum. The sharp next target is semantic minimality at the total-system
scale inside one large connected support.  There are now `2^h` states whose
`h` independent relation-cycle fluxes are separated by scalar queries on the
same connected bounded-degree or complete support at a fixed one-Walsh-block
`n^(3/2)` gap.  This proves `h` bits at fixed port accuracy, but the gap is
only `Theta(h^(-3/2))` in units of all `3hn` variables.  Determine whether a
different dense weighting exposes a positive-rate packing at total
`N^(3/2)` accuracy, or prove a further quotient at that scale.  The entire
bounded **state-local** architecture is now closed negatively: arbitrary
state-independent public bridges leave diameter at most `D/sqrt(k)` at total
scale, and unequal disjoint cells have bounded fixed-distortion packing
entropy.  To escape, a local bit must alter `Omega(k^(3/2))` unit interaction
atoms or enter a genuinely nonlocal state-dependent cross block.  The next
construction/lower bound must therefore target that broadcast regime; adding
more public connectors is no longer a legitimate experiment.

That broadcast regime is now attained in a sharply scoped model.  An
alternating-form evaluation code broadcasts `Theta(k)` hidden bits across a
constant density of constrained edges; one shared random dressing keeps
every child spectrally flat, and a predeclared negative-clone overlay yields
fixed `k^(3/2)` response separation.  The construction uses quadratic shared
public advice and its same-support overlay leaves the exact-sign class.  The
next discriminating target is to retain the positive rate under a disjoint
bridge composition or an exact-sign parent with subquadratic public advice,
or prove that one of those costs is unavoidable.
Interaction-mass deletion already gives an
`O(t/eta)` path carrier and an architecture-specific quadratic dense ceiling;
the next approximate theorem should beat or match that law without invoking
the original `2^(kn)` Boolean maximum.

A parallel replacement target remains:

```math
\text{Find }B_n\text{ with subexponential state and }
||A_n-B_n||_square=o(n^{3/2}),                                 \tag{OQ.5}
```

for one nontrivial structured dense family, or prove that every such
replacement needs exponential information. Generic Frieze--Kannan regularity
does not suffice: its state can already be exponential at this accuracy.

**Success:** an approximate reduction of the exact cross-information, a reusable joint quotient
for a nonlinear cap-bounded family, or a lower bound proving extensive
growth of every such quotient. Any positive state must retain the
commutation/compatibility data absent from scalar Rayleigh summaries and must
use a property stronger than bounded cap.

**Falsifier:** a response-separated code inside the proposed port, or a
four-label pinned cut-norm witness at order `n^(3/2)`.

## 6. Constrained compactness and realization

Fixed-interface response bodies have unrestricted finite recovery sequences.
Characterize which limiting response or presented-carrier states are realized
at all large sizes inside a constrained family such as linear codes, dense
CSPs, or bounded-width factor graphs.

**Success:** a Gamma-limsup/recovery theorem preserving the declared response
with vanishing normalized loss and without storing a target optimizer.

The abstract leading-rate case is now closed for unstructured landscapes:
compact microcanonical hypographs compose by supremal convolution and every
bounded usc profile has an abstract finite recovery sequence. This does not
solve constrained realization inside codes, sign quadratics, or graphs. It
also loses subexponential extremal multiplicities. The live alternatives are
a structured recovery theorem at the target scale, or a natural multi-speed
hypograph whose added resolution still closes under composition.

**Stop condition:** finite-state approximation outside the constrained model
does not answer this question.

## 7. Deterministic analogue of derivative-mass compression

The branching-random-walk benchmark proves that an exact prelimit pair
`(W,Z)` can collapse after critical renormalization to a scalar `Z` which
controls the entire unmarked limiting extremal process and composes by a
smoothing transform. This mechanism is not captured by finite-port
polyhedral response geometry.

Find a deterministic or adversarial landscape class with:

1. a declared extremal point-process query family;
2. a nontrivial centering and scaling;
3. a unique all-order limiting response law controlled by a strict state;
4. an exact limiting composition law; and
5. a falsifier showing which enriched future marks force the state to grow.

**Success:** a recovery/composition theorem in a model other than a branching
process, or a proof that one proposed deterministic class cannot have such a
state because its marked extremal measures are not tight or not unique.

**Stop condition:** subsequential point-process compactness without a unique
all-order recovery law is not a composable macroscopic state.

## 8. When does sparse semantic replacement become a congruence?

The finite-dictionary sparsification theorem gives the general upper law

```math
b\le m-\Omega\left(
\min\{m,E^2m/(V_Phi\log |X|)\}\right)                     \tag{OQ.8}
```

for one-shot all-future response, using a public list of weighted sparse
surrogates. It applies to quadratics, Littlewood polynomials, CSPs, and code
correlations, but the decoded family need not close under the model's own
composition.

Characterize dictionaries for which the masks and importance weights can be
chosen from an invariant algebra so that the same saving survives repeated
composition without a depth factor.

**Success:** a depth-stable sparsification theorem for two genuinely
different compositional models, or a response-packing lower bound matching
the `E^2/(V_Phi log|X|)` saving scale for one natural dictionary class.

**Falsifier:** a dictionary with small one-shot response cover for which
every invariant surrogate family has full coefficient information rate.

## 9. A nonlocal lower-spectral representation for growing interfaces

Finite-width discrete adversarial chains now have a complete scalar limit:
their pressure is a lower spectral radius, the ground-state density follows
by a uniform zero-temperature sandwich, and positive transfers admit a
contractive projective cavity net. This is not merely convex minimax and it
does not require a periodic optimizing disorder word.

The mechanism fails under the standard path decomposition of a dense
quadratic system. Every `n by n` sign bridge has

```math
max_(x,y)x^TBy\ge n^{3/2}/\sqrt3,
```

so the interface contribution is leading, its transfer dimension is `2^n`,
and its positive-temperature projective contraction degenerates with `n`.

Find a nonlocal semigroup, antinorm, or projective operator representation
for a growing-interface class whose dimension is strictly smaller than its
full Boolean landscape and whose multiplication defect is subleading.

The approximate-residual question is now closed at the generic level.
Terminal profiles admit a depth-uniform last-window theorem, but exact reset
maps can retain a positive directed compatibility rate, and arbitrarily
small residual shells contain scaled copies of arbitrary weighted response
algebras. Any positive scalar theorem must control the compatibility
cocycle, not only the residual radius or contraction coefficient.

**Success:** a restricted dense or rank-growing model with an exact scalar
lower-spectral limit plus a sub-landscape carrier, or a proof that one
natural proposed nonlocal product necessarily retains an extensive response
rate.

**Stop condition:** repackaging the `2^n` transfer matrix or convexifying the
discrete sign alphabet is not a new representation.

## Reconnection rule

Do not return to the motivating signing problem.  Reconnection requires a
carrier or synchronized quotient that arises naturally there, has controlled
growing-interface entropy, closes under the relevant composition, and has a
finite realization theorem.  The signed-balanced benchmark supplies all of
these only for a restricted permutation-invariant class; none is currently
proved for the motivating sign matrices.

## 10. Synchronization among amplification phases

Regular-Hadamard tensor powers and summably perturbed signed replications now
give genuine near-original limits along each prescribed hierarchy.  The
negative half is now solved sharply: every regular-Hadamard automatic prefix
hierarchy has a continuous mantissa profile, and the explicit order-four
Walsh hierarchy has

```math
L(1)=L(4)=1/2<L(3).
```

Thus convergence along every fixed phase and exact all-order realization do
not imply phase synchronization.

Logarithmic averaging is now completely classified and does not solve this
problem.  It always yields the Haar pushforward law of the phase profile, but
can identify distinct ordered profiles; every positive power-weighted mean
retains the full phase through an explicit differential inverse.

The remaining positive target is a checkable mechanism forcing this phase
profile to collapse.  One such mechanism is now proved: averaged Boolean
pullbacks inducing a Doeblin phase refresh collapse the limit when recovery
plus accumulated operator defect is little-oh of refreshed mass.  The Walsh
profile quantitatively forces this ratio to be at least `0.0605036` for a
specific full-support refresh law.

The remaining model-specific target is to find a **small-description**
certificate of this kind in a hierarchy not already uniformly operator-close
to one model, or a different inter-phase transfer estimate implying

```math
d_H(K_(d,r),K_(d',r'))\longrightarrow0
```

whenever the represented orders are asymptotically adjacent.

Uniform stationary expander refresh is no longer a small-description escape:
an observable phase excess at one-step toll `delta` requires
`Omega(1/delta)` bits, and the Walsh profile forces `Omega(sqrt N)` bits at
the natural transfer scale.  The live alternatives are nonstationary phase
laws, scrambling that vanishes with scale, or a transfer not represented by
the one-sided averaged pullback inequality.

**Success:** a pullback/refresh certificate with sub-landscape description
complexity in a dense hierarchy not already uniformly operator-close to one
fixed model, or a theorem deriving refresh from natural pseudorandomness.

**Stop condition:** a condition that simply assumes constancy of the phase
profile or vanishing Hausdorff distance is not a recovery theorem.  The
Walsh construction is now the mandatory falsifier for weaker proposals.

## 11. Tangent mass for multiscale extremal convolution

Finite lexicographic valuations compose exactly on fixed descriptor
alphabets, but Vandermonde convolution gains a logarithmic saddle-mass term
that pointwise coefficients miss.

The nondegenerate Gaussian target is now closed.  A fixed-query discrete
Laplace theorem supplies the Hessian and lattice-covolume factor; the
finite-parameter Gaussian tuple closes associatively and has finite integer
all-order recovery.  This repairs Vandermonde and multinomial tangent mass
without storing the descriptor grid.

The fixed power-exponential extension is now falsified.  Its leading
homogeneous roofs close for every `p>1`, but central and off-centre tangent
exponents differ, and the normalized power-exponential tangent family is
closed under self-convolution only for the Gaussian power `p=2`.

The remaining target is a **different finite stratified extension**:
identify a finite menu of non-Gaussian saddle types whose parameters remain
closed when saddles merge, split, or become degenerate.

**Success:** a class containing at least one genuinely nonquadratic saddle
and stable under two successive convolutions, with a finite parameter count
and integer recovery.

**Stop condition:** arbitrary functional amplitudes, the entire descriptor
measure, or a separate postulated asymptotic for every query are not strict
state reductions.  The quartic `n^(3/4)` example is mandatory.
