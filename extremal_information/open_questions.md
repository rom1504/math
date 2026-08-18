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

The residual-error half of **approximate residual/core synchronization** is
now solved on a common-law contracting carrier:

```math
\text{how do profile-covering error, cycle toll, and support leakage combine}
\quad\text{without paying the same approximation at every depth?} \tag{OQ.1}
```

Terminal error decays as `rho^t`, fresh centred reward error costs
`B/(1-rho)`, and scalar recurrent error is exactly the cycle-mean norm.  A
fixed strongly connected carrier consequently has response code size
`N_U(epsilon)N_B((1-rho)epsilon)(1+2L/epsilon)^r_G`, with a matching cycle
exponent.  The unresolved structural target is to **derive**, rather than
assume, a common transported contracting law or an approximate replacement
from natural max-plus/Bellman generators.  If supports are also claimed,
charge the greatest-core leakage separately: the response theorem does not
turn proof memory into semantic state.

**Success:** a checkable generator-level condition producing the certified
carrier for a class broader than promoted rank-one/reset systems, or a
pumpable counterexample proving that no such law can exist.

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
bounded-fan-in alternative is now completely classified: fan-in `t` permits
at most `O(t sqrt k)` total-scale hidden coordinates, and unrestricted exact
signings attain this scale with a matching pairwise packing.  Moreover the
flat alternating-form code necessarily has linear fan-in in every hidden
basis.  Concise public advice is no longer an obstruction: a small-bias label
sampler and a `6k`-wise-independent polynomial seed reduce the shared
description to `O(k log k)` while preserving flatness and the full response
packing.  The next discriminating target is therefore to retain the positive
rate under a **bounded-cap** disjoint bridge composition.  Bare exact-sign
and disjoint metric closure is now solved: a rank-one bridge plus a positive
auxiliary clique exposes every child coordinate exactly on `N=2k` vertices,
preserving the response packing.  Its common calibration is `Theta(N^2)`.
Moreover every future that universally pins one prescribed coordinate
against all sign children must have `Omega(N^2)` cap.  The live cost is thus
not exact signs alone but removing the quadratic calibration through
child-dependent optimizer switching or a genuinely joint absolute channel.

The first compiler boundary is now sharp.  One auxiliary per old edge gives
an exact sparse identity, but its quadratic order destroys the total scale.
Among independent-star futures, endpoint-local and fully dense selectors
both require quadratic order; bounded intermediate fan-in obeys
`sum_a d_a^(3/2)>=k(k-1)`.  Arbitrarily interacting auxiliaries still need
at least `0.565...k-o(k)` states for the universal complete cut shell, and a
bounded-cap linear parent cannot approximate that shell because its response
oscillation is only `O(k^(3/2))`.

These no-go results do **not** settle the flat alternating-form family.  The
precise live target is now one of the following mutually discriminating
statements:

1. construct an `O(k)`-vertex, complete-sign, jointly interacting compiler
   whose parent cap is `O(k^(3/2))` and which preserves a fixed fraction of
   the alternating-form response packing; or
2. prove that every such compiler needs linearly many selector bits with an
   incompatible slope/cap budget.

The signed exposed-set input is now proved in full generality: one orientation
of every nonzero quadratic has a constant-fraction antipodal exposed bulk.
It gives `Omega(k^(2/3))` auxiliaries without cross-block flatness and
`Omega(k)` with it.  The remaining logical gap is orientation: if the old
child stays outside the selector envelope, the argument sees the residual
query rather than the child--query difference.  This cannot be dismissed,
because an explicit flat block-clique family has a target-scale absolute
maximum but exponentially small absolute near-top mass.

The full quadratic-character equality-lock route is also closed.  Exact
pullback of every pair character forces only a signed permutation and common
gauge, and every fixed complete bridge then has an `Omega(k^(3/2))`
worst-case defect, even with bounded coordinate replication.  A bare bridge
robust against all spectrally bounded sign children is forced to rank one and
quadratic cap.  Therefore the remaining positive route may preserve only the
alternating-form query subalgebra or must let interacting auxiliaries and the
child energy choose a joint witness; transporting the full pair algebra and
locking it afterward is no longer live.

Bounded cap by itself does not force child-dependent switching.  A fixed
one-spin future has a common optimizer across `2^(Theta(k^2))` bounded-cap
children.  The response-sensitive theorem is a dichotomy: the flat Gram
packing needs an `Omega(k)` approximate witness dictionary, realized either
as common query pins or child-dependent active states.  Coordinate pins can
be much more expensive still.  For a switched spectrally flat exact-sign
pair with block-clique contrast, every fixed-fraction extremal pin library
has size `exp(Omega(k))`; reciprocal extremal mass gives the sharp exponent
up to a polynomial factor.  The unresolved question is whether one bounded-
cap **nonlinear** context can aggregate those rare coordinate witnesses
without reconstructing the full absolute landscape.

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

## 12. Reusability of bounded-cap anti-pin compilers

The regular-Hadamard anti-pin theorem gives an `exp(Omega(n))` family of
exact cap-`n^(3/2)/2` children and exact order-`n+sqrt(n)` contexts whose
projective response metric has constant distortion.  The context shore is
asymptotically minimal within its repeated-rank-one architecture.

The first multi-port boundary is now known.  At fixed port count, two Gram
matrices give a polynomial-size spherical certificate.  They do not arise by
taking the product of separate one-port states: orthogonal top eigenvectors
create a leading cross-Gram gap.  More strongly, `A` and `-A` are uniformly
one-port close but a second macroscopic Walsh shore exposes their relative
orientation by `(2-sqrt(2))n^(3/2)`.  The orientation quotient remains valid
for every flattened continuation of total width `o(n^(3/4))`.  This exponent
is now sharp for unrestricted sign futures: a biased flat bridge and clique
at width `n^(3/4)` attain a leading gap.  Futures which themselves have
`O(m^(3/2))` cap still require linear width.

The fixed-scale collective metric question is now closed for the relaxed
PSD carrier.  Trace-bounded spectral truncation and factor nets give
`exp(O_eta(p))` covers in the signed quadratic query metric, while Boolean
rank-one ports give `exp(Omega_eta(p))` packings.  At bounded total port mass
this yields a linear-rate spherical response carrier with a sharp
square-root hard-edge modulus.  Thus the `O(p^2)` exact cross-Gram table is
not the operational fixed-accuracy state.

**Success:** a lower bound exposing the true number of collective cross
variables at the correct total-system scale, a reusable approximate gluing
law for the low-rank factor carrier, or a deterministic Boolean recovery
theorem transferring its spherical response to exact old-spin caps without
a fixed leading loss.

The third item now has a strict query-local solution.  One near-global
spherical optimizer with vanishing coordinate-flatness deficit rounds with
`o(n^(3/2))` loss, and a close-pole rank-two Walsh family satisfies the
condition exactly.  Orthogonal rank-two poles retain a fixed gap, while a
uniform Boolean net for an entire multidimensional trust span is impossible.
The live problem is no longer the existence of a recovery condition; it is
whether exposure and flatness are **reusable under switching composition**
without tabulating all trust channels.

Two further boundaries now sharpen that question.  A common Boolean pole
has an exact multiplicative correlation-deficit algebra, but equal deficits
need not have equal responses.  A supplied exposed optimizer has an exact
block flatness chain rule: local nonflatness and RMS-amplitude allocation add
as separate positive resources, and fixed allocation imbalance pumps to a
maximal gap.  Uniform recovery of a whole `d`-dimensional active sphere needs
an explicit near-eigen witness cover of order
`(c/sqrt(epsilon))^(d-1)`.  The remaining target is a switching law that
selects and transports only the exposed subset without assuming its orbit
language in advance.

One narrow synchronization law is now exact.  If every onsite block is
`+-H` and every bridge is `+-H` for one common regular Hadamard factor, the
projective signed-graph class closes: `r` joining edges create exactly `r`
compatibility bits, split into relative antipodes and cycle flux.  Any
proposed general theorem should recover this fibre law without assuming the
whole scalar block graph as unexplained input.

The relaxed PSD amalgamation problem is also now classified at fixed scale.
Two marginal factors glue through a contraction; spectral truncation makes
the conditional compatibility cover port-independent, and global
multi-piece realizability is one block PSD constraint.  This does **not**
subsume signed coefficient holonomy: a rank-one global Gram alignment has
positive cycle products, whereas the negative cycle witness above lives
outside the PSD carrier.

The Boolean information boundary is now sharper.  The complete labelled
linear port table is an invertible convolution of the projective row
histogram.  Its exact dimension is `2^(p-1)-1`, but its fixed-error metric
entropy is exactly `Theta_eta(p)`, and declared union/tensor occurrence trees
have mergeable `O(p/eta^2)`-bit randomized carriers.  Uncontrolled semantic
reuse can accumulate linearly, while a uniform convolution component gives
geometric forgetting.

Pairwise Gram--Rayleigh data are definitively insufficient: equal states in
one regular-Hadamard tensor family have Boolean trust responses separated by
`rn/8`.  A positive replacement is also known.  Odd product-algebra closure
makes the row histogram an exact trust carrier, and a dense family with
`p=(log_2 n)/2+1` has an `O(sqrt(n)log n)`-bit exact state at bounded port
mass.  The live question is no longer whether a sub-landscape Boolean carrier
can exist.  It is whether a comparably generated affine pole algebra, or an
approximate selector closure with summable defect, is forced or can be
constructed for a class broad enough to interact with near-minimizers.

**Mandatory falsifiers:** the equal-`(G,R)` four-port collision, the
orthogonal-top-eigenvector cross-Gram pair, and the `A` versus `-A`
macroscopic orientation bridge.  Merely showing that one planted Walsh
witness moves is already covered by the Hadamard synchronization obstruction.

**Stop condition:** storing the complete multi-query response table, the
entire Boolean optimizer language, or one fresh `n`-bit state for every
composition node is not closure.  Full odd-product closure counts as progress
only when it has a generated sublinear presentation, as in the tensor coset;
postulating all `2^(p-1)` identities at arbitrary dense arity does not.

## 13. Near-minimizer selector synchronization

PC.3 meets the literal near-original benchmark: it has a dense sign bridge
with `Theta(n^(3/2))` edges, growing rank, leading response oscillation, and
a strict sub-landscape state.  It does **not** authorize a direct attack on
`M_n`.  The new shore has only `Theta(sqrt n)` vertices, the child is one
regular-Hadamard tensor hierarchy, and exact positive-pole closure already
fixes its child cap at `n^(3/2)/2`.

The weakest useful reconnection question is instead the following agent-
authored hypothesis.  Predeclare, independently of the signing, a family of
at most `exp(o(n))` admissible port systems.  Along a near-liminf sequence
`A_n`, can one select an `o(n)`-bit-indexed member, without enumerating its
full maximizing-spin landscape,

```math
p_n={1\over2}\log_2n+O(1),
\qquad m_n=\left\lfloor{\sqrt n\over p_n}\right\rfloor,
```

Boolean ports `W_n` and one fixed antipodally odd selector `tau` such that,
for every endpoint `epsilon`, its coordinatewise selector witness
`x_epsilon` obeys

```math
Q(A_n)-|H_(A_n)(x_epsilon)|=o(n^(3/2)),            \tag{NMSS.1}
```

while the declared query language is nonvacuous,

```math
osc_epsilon\;m_n||W_nepsilon||_1
>=c n^(3/2)                                        \tag{NMSS.2}
```

for one fixed `c>0`?

If so, every labelled response differs by `o(n^(3/2))` from

```math
Q(A_n)+m_n||W_nepsilon||_1,                        \tag{NMSS.3}
```

and Theorem 21.50 makes `(Q(A_n),mu_(W_n))` a
`O(sqrt(n)log n)`-bit leading-scale carrier.  The at most `sqrt n` new
vertices can be completed at cost `O(n)`.  This would be a strict response
reduction for actual liminf objects, but would still require a separate
scale-changing reuse theorem before saying anything about convergence.

NMSS has a strictly query-narrower conclusion than the original optimization:
it controls only
`Theta(sqrt n)` selector words and one rank-`O(log n)` continuation language,
not all `2^n` spins or any cross-order recurrence.  It becomes circular if
the port description itself encodes a maximizing-spin table, if Gram data
replace the histogram, or if compression of the completed parent is mistaken
for proof that the parent is itself near-minimizing.  More cautiously, its
conclusion is query-narrower but is not yet known to be mathematically easier.
The robust spectral product theorem approximates the roof `rn/2`; it proves
NMSS only in a class where `rn/2-Q(A_n)=o(n^(3/2))`, and therefore cannot be
silently applied if the true asymptotic cap lies below the spectral roof by a
fixed amount.

## 14. From labelled full response to an unconstrained parent

The weighted majority-tail construction proves that individual product
deficits can vanish while one joint selector has constant loss.  At its
diagonal-completed Hadamard roof, the new exact-sign order-16 seed proves the
same logical failure with marginal deficit `3/16`, and common-top tensoring
preserves it at every order `16^j`; hollowing preserves its Boolean energies
but uses a slightly different valid contraction roof.  That fixed seed did
not settle the combined asymptotic regime; the sparse-flip theorem below now
does at certificate level.

The certificate question is now answered.  For the public PC.3 family with

```math
p_N={1\over2}\log_2N+O(1),
```

choose `r_N>=||A_N||_op` with `r_N=Theta(sqrt N)` and define
`R=Z^TA_NZ/(r_NN)`, `D=G-R`.  The sparse-flip construction proves

```math
\max_S D_(SS)=o(1),
\qquad
Delta_(Maj)=\kappa+o(1).                          \tag{VEC.1}
```

The labelled **full** Boolean trust-response question is now answered.  Let
`h_j=W_jepsilon_j`.  Corollary 21.63 proves

```math
||h_j||_1=\sqrt{7/(2\mathop{\rm pi})}\,N_j\sqrt j+O(N_j).
```

Thus `m_j=floor(lambda sqrt(N_j/j))` makes
`m_j||h_j||_1=Theta(N_j^(3/2))`, while the new shore has only
`O(sqrt(N_jj))=o(N_j)` vertices.  Even its whole internal signing contributes
only `O(N_jj)=o(N_j^(3/2))` energy.
Compare

```math
\max_(y\in\{+-1\}^{N_j})
\left\{{1\over2}y^TA'_jy+m_jy^Th_j\right\}       \tag{VEC.2}
```

with the unflipped value, where `x_epsilon` pays both terms exactly before
the sparse flips.  Theorem 21.64 uses a spherical relaxation to prove a
fixed normalized gap; at `kappa=1/2,lambda=1/10` its explicit lower bound is
greater than `0.0146`.  Thus no alternative Boolean child spin repairs this
fixed field.

The realization question is now sharper still.  The direct free shore is
falsified: its all-positive endpoint has field norm at least
`N_j(1+j/2)`, hence cap `Omega(N_j^(3/2)sqrt j)` at this multiplicity.  A
shore-only completion cannot remove that bias because its entire energy is
`o(N_j^(3/2))`.  Theorem 21.66 now gives a balanced rowwise microcanonical
compiler with global cross cap `||m_jh_j||_1+o(N_j^(3/2))`.

The strengthened uniform affine estimate in Theorem 21.66 answers endpoint
stability positively, and Theorem 21.67 yields an unconstrained exact-sign
parent cap collision.  Theorem 36.1 now settles the qualitative one-hot
amplification question: a common physical query bank separates
`N_j^c` states, or `Theta(log N_j)` message bits.  It does **not** build a
`2^{k_j}` independently writable edit cube.  The remaining Phase-I question
is whether nonlinear multiplexing can beat the conditional disjoint-edge
and canonical rank-one-superposition ceilings; do not spend near-minimizer
time optimizing the crude one-hot exponent.

Two natural searches are already closed.  The PC.3 pole-conjugation twirl
has a depth-independent gap, so every diagonal switching has selector defect
at most twice its mean marginal defect.  And every monomial tensor/direct-
mixture amplification of the fixed five-port seed that reproduces its
prescribed selector blockwise keeps some odd active channel at deficit at
least `3/32`; the common-factor subclass has the sharp
`3/8` selector relation.  The successful sparse-flip construction is
genuinely non-switching and acts at leading edge scale, exactly the escape
left by those no-go theorems.

On the positive side, bounded scalar visibility was not the right necessary
target: it decays in PC.3 even while the collective twirl stays coercive.
For future positive classes, seek a natural exact-sign hypothesis that
**implies**, rather than assumes, a uniformly gapped, low-description
positive observable whose common kernel is an eigenvalue-one space for every
declared selector conjugate.  Theorem 21.60 supplies the conclusion once that
collective object is given.

## 15. Near-minimizer structural frontier

Theorem 36.2 proves the first arbitrary-order collective consequence of
near-minimality: an `epsilon`-near-minimizer has a shell of width
`O(sqrt(epsilon))n^(3/2)` carrying one signed-cut law with total edge bias
`O(sqrt(epsilon))`.  The exact-minimizer shell also has
`exp(n^(1/2-o(1)))` members at vanishing normalized width.

The smallest positive question is no longer to find any collective law.  It
is:

```text
Does the forced edit-thick shell law have joint selector coercivity
f(epsilon)+o(1) on one nonvacuous O(log n) port language,
after o(n^(3/2)) response residual?
```

An affirmative answer plugs directly into Theorems 21.50, 21.53, and 21.54
and yields `o(n^(3/2))` response error with a sublinear reusable state.  The
claim must control the Fourier-product quadratic form jointly; another
coordinatewise first-moment estimate is insufficient.

The strongest backup asks instead for an edit-thick common-law contracting
fibre whose persistent information lies in `o(n)` terminal and recurrent
cohomology coordinates.  This allows a small number of hidden coherence
bits rather than requiring them to vanish.

Counterexamples impose four mandatory qualifications:

1. the shell width must pay at least the optimality excess;
2. spectral claims must peel or otherwise quotient localized edit implants;
3. exact active multiplicity and adjacent-order deletion optimality are not
   stable;
4. full amplitude-`n` pinning may retain linear information even inside a
   near-minimizer halo, so any positive transfer must declare a balanced,
   low-cap context class rather than silently answering the full landscape.

## 16. Frontier after the multiscale affine-shell theorem

The previous selector-coercivity target must be revised twice.

First, an exact-sign vanishing spectral-roof defect already forces the
conjectural constant `1/2`; it is not demonstrably weaker than the original
problem.  Second, FB.1 first-marginal balance alone cannot choose one
spectral orientation, as the Walsh exact-shell example shows.

At the same time, Theorem 36.7 proves something stronger than the proposed
cap-relative shell algebra for **every** bounded-cap signing: an
orientation-pure affine cube of dimension `n/q`, with one-sided selector
error `8Q/q`.  Thus neither shell entropy nor coherent odd-product closure
is the remaining near-minimizer property.

The narrow unresolved alternatives are now:

1. **Joint physical cancellation.**  Can a balanced exact-sign bridge be
   chosen so its omitted channels cancel against the quadratic child before
   absolute values, beating the separately paid `n sqrt s` residual on a
   fixed-ratio interface?
2. **Cross-level congruence.**  Can favorable affine charts be selected by a
   rule whose gauges and response fibres survive composition, rather than
   being recomputed independently at every order?
3. **Low-cap incompressibility.**  Can the linear pinned-response packing in
   the genuine near-minimizer halo be compiled into bounded-cap all-spins-free
   contexts, proving that even near-minimizers retain an extensive response
   rate?

None currently qualifies as a proved strict reduction.  Projection plus
separately paid residuals is under an explicit target-scale audit; any
positive formulation must state the joint cancellation or nonlinear
congruence rather than hide it inside “contracting fibres.”

Theorem 36.9 sharpens the second alternative.  Every exact minimizer already
has a balanced atlas of `O(n^(1/3))` signing-dependent charts, each with
`Theta(n^(5/6))` affine dimensions and `O(n^(4/3))` one-block response
error.  The next theorem must therefore compress **cross-chart transition
data** or give a common low-cap physical context.  More within-chart shell
states, even with a common orientation, do not advance this question.

## 17. Fixed-projective-gap packing in deepest cut-code cosets

Theorem 36.11 removes the physical-compiler obstruction.  The selected
near-minimizer question is now purely structural:

```math
\exists\gamma>0,\quad d_n=o(n^{3/2}),\quad K_n\to\infty
```

such that every exact minimizer has `K_n` positive `d_n`-shell atoms with
pairwise absolute edge overlap at most `1-gamma`.  The preferred quantitative
form is `K_n=exp(Omega(n))`.

In cut-code language, this asks for projectively separated nearest or
near-nearest codewords in every deepest coset of the augmented coboundary
code.  It is strictly smaller data than the full coset histogram and has an
exact finite falsifier.  A positive answer compiles immediately into a
growing all-spins-free contextual packing of total cap `O(n^(3/2))`.

The mandatory obstruction is the two-cap alternative: an FB-balanced shell
may put nearly half its mass close to one projective cut and nearly half near
the opposite signed orientation.  Positivity rules this out only below
radius `Theta(n^(-1/2))`.  Either prove that deepest cut-code cosets exclude
this geometry at fixed radius, or construct it scalably.  Do not substitute
raw shell cardinality: Theorem 36.10 gives
`2^{Omega(sqrt n log n)}` thin-shell witnesses universally without a
projective packing theorem.

The exact-minimizer quantifier is indispensable.  Theorem 36.12 constructs,
for every prescribed `Delta_n=o(n^(3/2))`, a sequence with
`Q(B_n)=M_n+o(n^(3/2))` whose whole positive `Delta_n` shell lies in one
vanishing projective cap.  Thus no uniform version over all
`epsilon_n`-near-minimizers can hold when the admitted shell is narrower
than the optimality slack.  Any proof of the selected statement must use a
genuinely discontinuous consequence of exact global optimality.

Theorem 36.13 isolates the cut-specific branch more sharply.  Generic
deep-hole theory and every finite-flip cover identity allow a scalable
collapsed countermodel.  In the augmented cut code, however, a genuinely
vanishing opposite-lift cap must be supported by a complete bipartite
interface `delta(S)` with

```math
\Omega(\sqrt n)\le |S|=o(n),
\qquad
\sum_{e\in\delta(S)}a_e\ge M_n-o(n^{3/2}).
```

Can exact minimization force a third near-ground projective direction away
from both poles of every such mesoscopic interface?  This is the concrete
cut-specific form of `Cut-DH(3)`.  Replacing the interface by a separately
optimized rectangular signing merely reopens the archived joint-bridge
obligation and is not an answer.

At finer radii, Theorem 36.14 gives a second discriminating target.  A small
shell cover below `n^(-1/4)` has a local-field response roof; a large packing
above `n^(-1/4)` is physically visible.  What extra invariant controls a
shell whose metric entropy lives at the critical `n^(-1/4)` scale?  Any
answer must beat the quadratic internal-edge Taylor term or the exact-sign
compiler fluctuation without retaining the complete coefficient matrix.

## 18. Higher-order non-recycling in exact-minimizer shells

Higher-order non-recycling from an arbitrary finite anchor set is now proved
at the energy scale.  Theorem 36.19 replaces literal common-correct
intersections by a diffuse fractional reservoir and gives

```math
L_n=Theta\left({\log n\over\log\log n}\right)\longrightarrow\infty
```

words in one `o(M_n)` positive shell, pairwise separated by
`(1/4-o(1))M_n`.  Negative four-anchor holonomy is therefore a ceiling only
for literal intersections, not for jointly weighted localized flips.

This is not the strongest known global packing count.  Theorem 36.7 already
gives `exp(Omega(sqrt n))` such energy-scale-separated shell words with
`O(n)` deficit by taking `q=Theta(sqrt n)`.  Accordingly, improving the
fractional conditioning constant is a question about flexible sequential
extension, not a route to a shell-cardinality frontier which the affine cube
has already passed.

Two quantitative questions remain internal to this mechanism.  Is the
worst finite-phase conditioning constant for **actual exact-minimizer shell
patterns** substantially smaller than the arbitrary-sign-pattern bound?
And can the uniform response approximation be obtained with fewer than
`Theta(C_K^2n)` physical flips?  Either improvement would accelerate the
packing count, but neither alone changes its ambient separation scale.

The whole-shell version is already false as a naive finite heuristic:
Example 164 gives exact minimizers through orders 13 and 14 whose complete
active-shell fractional LP forces `w_e=1` on every edge.  This does not
answer the sequential growing-anchor question, but it rules out treating
fractionalization alone as a shell-wide `O(M_n)` quotient.  Any positive
asymptotic strengthening must exploit how anchors are selected, not merely
that all of them have positive response.

Even a growing packing at distance `Theta(M_n)` does not plug into the
current physical compiler at target response scale.  The fixed-gap version
of `L_projective`, or a different compiler which jointly amplifies these
energy-scale directions without raising total cap, remains the operative
frontier.

## 19. Nonlinear two-cap entry beyond deficit-scale diffusion

Theorem 36.17 proves that low dimension forces projective shell diameter
`(2-o(1))t`, where `t` is the covering-radius deficit.  Can an antipodal
`[N,k]` code with `t=Theta(sqrt(Nk))` have every `o(t)` shell contained in a
nonlinear, affine-subspace-evasive union of two Hamming caps of diameter
`o(N)`?  Direct products, large affine carriers, and a single sub-threshold
cap are now obstructed; no generic theorem controls the opposite lift.

For the augmented cut code, the corresponding question is whether exact
minimality converts repeated opposite-lift entries into a fixed-scale
projective direction.  Any proof must use more than dimension, finite-flip
covering, pairwise reservoirs, or shell cardinality.  A scalable nonlinear
countermodel would sharply delimit what generic coding theory can contribute
to `L_projective`.

## 20. Intrinsic rare-event states beyond generic spectral synchronization

Theorem 37.1 gives a strict bulk-plus-spike response algebra only for a
declared presentation and predeclared generic directions.  Is there a
nontrivial deterministic class in which a finite rare-event state is
intrinsic to the landscape and remains closed under every allowed
continuation?  What hypothesis makes rooted spectral measures functions of
a smaller unrooted state without Haar genericity?

A useful answer must survive correlated perturbations or state precisely
the extra relative-geometry marks they require.  Retaining the full
eigenbasis or the full family of rooted measures merely reconstructs the
matrix and is not a strict quotient.

## 21. Activation-cost states beyond hard marked phases

Theorem 37.2 exactly identifies the bounded-temperature pressure and replica
budget needed to observe one marked phase.  Can an intrinsic finite carrier
do the same for several competing phases under an interaction that is more
general than hard conjunction?  The carrier should compose without retaining
the complete energy histogram, while preserving both activation costs and
the geometry needed by the declared future queries.

The immediate falsifiers are labelled futures, which distinguish equal-mass
marked sets, and ordinary additive composition, which creates intermediate
phases.  A positive theorem must specify a class whose phase table closes or
whose synchronization collapses it.  Merely restating a finite-temperature
large-deviation rate function is not enough; the new state must have an
exact or quantitatively controlled update law.

Theorem 37.3 adds a sharper stopping condition.  Even a uniform fixed-rate
endpoint tail, exact variance, fixed-parameter monotonicity, and exact
centered scalar subadditivity permit oscillation of every fixed-temperature
diagonal.  Thus neither a scalar activation cost nor one endpoint entropy
deficit should be promoted again.  The next admissible target must couple
energy shells across orders or couple spin entropy to bridge-disorder
geometry in a way not reducible to the scalar pressure curve.

Theorem 37.4 further rules out using the new tail theorem as a proxy for a
bridge basin: the fixed-small-tilt conference defect survives after
conditioning on an overwhelming class where every exact-sign output has the
uniform tail.  The concrete surviving question is to determine the
large-deviation speed and support edge of low parent pressure over bridge
disorder.  A useful theorem must show an `exp(-O(n))` favorable basin or a
joint spin--disorder contraction; an `exp(-Theta(n^2))` isolated optimum
cannot be recovered by fixed tilt.

## 22. Root-to-state broadcast at the energy scale

Theorems 36.20--36.22 separate three formerly conflated tasks.  Every
`Theta(M_n)` shell pair has a bounded-cap all-spins-free **root selector**.
Every exact minimizer also has an exponential switching orbit at mutual
Boolean distance `Omega(n^(3/2))`.  But the existing AO affine map sends a
mesoscopic root into children only `o(n^(3/2))` apart, after which no common
future can amplify it.

Is there an exact-sign map from a selected shell root to a child-owned block
which has total cap `O(n^(3/2))`, broadcasts `Theta(n^(3/2))` uniform
distance, and retains strictly less than the full rooted landscape?  A
positive statement must name the compressed input and a common update law.
A negative statement should prove that every such map reconstructs a
root-dependent child optimization or pays a fixed leading loss.  More pair
queries do not address this question.

## 23. What can information heaviness contribute after `L_tail`?

The former `L_tail` question is solved.  Theorem 36.26 proves a stronger
statement: every complete signing with `Q(A)<=Cn^(3/2)` has a fixed-rate
thin tail near both one-sided extrema.  Theorems 21.8 and 36.25 therefore
give every exact minimizer an `exp(Omega(n))` switching suborbit separated
by `Omega(n^(3/2))` scalar responses under all-spins-free exact-sign,
bounded-cap contexts.  The switching label gives a matching `O(n)`-bit
description for this orbit.

The open issue is no longer whether near-minimizers are response-heavy, but
whether this lower bound is relevant to cross-order optimization.  Can an
`O(n)`-bit state, necessarily retaining the switching gauge, close under a
fixed-ratio composition with `o(n^(3/2))` error?  Or can one prove that
arbitrary dense bridges expose additional information beyond every
linear-bit state?  Either result must specify the future language and an
update law; the scalar contextual packing alone proves no recurrence.

Theorem 36.29 answers the switching-label part negatively.  When the bridge
is optimized covariantly for each parent, every child switch is pure gauge
and the entire `Theta(n)` public response packing collapses to one state.
After root-gauging both complete child classes, however, the residual exact
bridge fibre still has `2^(mn-1)` coefficient classes.  The live question is
therefore gauge-invariant: can pressure/support-edge information over this
anchored fibre be compressed, or can one prove a leading response packing
between distinct bridge holonomies or distinct child switching classes?

The internal anatomy of a hypothetical smaller fixed-level principal core
is also now sharp.  It must be positive-dominant, have an oppositely
oriented near-minimal complement, and may nevertheless be diffuse and far
from every switched clique (Theorem 36.27).  Excluding such a core is not
needed for the thin-tail theorem and should not be pursued unless a new
cross-order mechanism makes its interface observable.

## 24. What is the lower-deviation speed of conference bridge pressure?

For the deformed bipartite Rademacher pressure in Theorem 37.5, prove or
disprove that, for every fixed `delta>0`,

```math
\Pr\{L_{\epsilon,B}\le(h_\beta-\delta)r\}
\le\exp\{-r\,\omega_{\beta,\delta}(r)\},
\qquad \omega_{\beta,\delta}(r)\to\infty.
```

A speed-`r^2` theorem would match the conjectured lower-tail scale in some
mean-field spin glasses; an `exp(-O(r))` construction would instead locate
a finite disorder-temperature transition.  Either result must handle the
deterministic conference diagonal blocks, Rademacher cross disorder, and
finite-temperature pressure.  Existing Gaussian upper-deviation results do
not transfer.  The useful next object should be a concrete joint
spin--disorder statistic or support edge, not a reformulation of the same
negative moment.

## 25. What can a speed-`r` conference basin look like after regular and quartic localization?

Theorems 37.6--37.8 rule out three broad mechanisms for an
`exp(-O(r))` favorable basin at small fixed temperature:

1. it cannot lie in the operator-regular class, nor within
   `Theta(r^(3/2))` Hamming distance of that sign set;
2. on the FMW-power-regular class it cannot reach the child target without
   moving the quartic completion statistic to its quadratic-entropy edge;
3. under any `O(r)`-relative-entropy change of bridge law it cannot lower
   the frame potential or full quartic coordinate by a fixed leading
   fraction.

Determine whether the remaining deeply operator- and power-irregular sector
contains a target-reaching family of probability `exp(-O(r))`.  A useful
positive answer must exhibit a counted exact-sign family and a rigorous
pressure upper bound.  A useful negative answer must prove that every
target-reaching bridge either has a fixed quartic deficit or pays another
explicit `omega(r)` information charge.  Merely adding higher moments one
at a time is not acceptable unless they close under a uniform pressure
remainder.

Equivalently, isolate a finite or compact gauge-invariant support-edge state
whose low-pressure fibres are superexponentially small, or construct a
linear-entropy higher-cumulant cancellation showing that no such bounded
state can suffice.

## 26. Can a nonlinear diffuse bridge law make low pressure typical?

Theorems 37.9--37.10 exclude concentrated conditional information, every
`O(r)`-codimensional affine fibre with a small repair, independent weak
template tilts, exact finite-template type shells, and adaptive switching
cross-sections as conference pressure basins.  The remaining candidate must
simultaneously have:

```text
D(q_r||U_r)=O(r),
s_*(q_r)=Theta(r^2),
operator-irregular mass bounded away from zero,
and E_q f <= (h_beta-eta)r.
```

Construct such a law with a pressure certificate not defined by thresholding
`f`, or prove it cannot exist.  The first concrete test is a rowwise
sign-invariant magnitude fibre, for example one fixed-probability constraint
per row on `|<B_i,v>|`.  It has the right entropy and is not a switching
cross-section; its dependence changes one population-covariance direction
while leaving most local moments typical.

A useful negative theorem should cover a natural class of dependent diffuse
laws—perhaps through a rowwise invariance principle, a convex-order
comparison, or a uniform quenched/annealed estimate.  A bounded-difference
tail inside the conditioned class stops at speed `r` and does not answer the
question.

The row-product test is now closed.  Theorems 37.11--37.14 prove that every
centrally symmetric bounded-Renyi-two row product, and every joint law within
`o(r)` row total correlation of it, cannot lower pressure anywhere in the
strict conference interval; at smaller beta, every event below an explicit
linear entropy budget is typical regardless of its structure.  The remaining
law must carry `Omega(r)` genuine cross-row information and be irreducible
after every cheap low-rank response peel, not merely diffuse and
operator-irregular before projection.

## 27. Arbitrary constant-density row fibres: solved positively

Theorem 37.13 answers the former question and proves more.  Every row law
with uniformly bounded Renyi-two density admits a deterministic
`O(sqrt r)` mean-plus-covariance peel, an `o(r)` Frobenius coupling on the
complement, subcritical removed nuclear mass, and the sharp projected edge

```math
\|B_r(I-P_r)\|_{op}/\sqrt r\longrightarrow2.
```

The nontrivial last step follows by verifying the Strong Tail Projection
hypothesis and combining a dependent-row upper-edge theorem with the
Marchenko--Pastur bulk lower edge.  Hence no iid constant-density row fibre
can be the favorable conference basin.

## 28. Can conditionally singular cross-row information create the remaining phase?

Theorem 37.14 extends the no-gain result through arbitrary dependence whose
row total correlation is `o(r)`.  A favorable law with common
constant-density marginals must therefore satisfy, simultaneously,

```text
D(q_r||mu_r^{otimes r})=Omega(r),
s_*(q_r)=Theta(r^2),
and E_q f <= (h_beta-eta)r.
```

Theorem 37.15 additionally closes every common-latent mixture of product row
laws whose component Renyi-two constants are uniformly tight, even when its
total correlation is linear.  Construct a law with a pressure certificate
independent of thresholding `f`, or prove a response-specific inequality for
the remaining class.  A latent-mixture candidate must put nonvanishing mass
on component row laws whose Renyi-two constants escape every fixed bound;
alternatively, its dependence must remain irreducible after conditioning on
every proposed latent state.

The exact positive compactness statement is now Theorem 37.17: tightness of
the random component `D_2` is enough, and fixed-`K` convergence supplies some
nonconstructive slowly growing window.  A named growth-rate theorem would
need new quantitative pressure input; it is no longer the selected target.

Theorem 37.16 supplies the first upper obstruction to that growth window.
At

```math
\log K_{2,r}=Theta(r\log\log r/\log r)=o(r),
```

an exactly isotropic Hadamard-cluster law has a super-Bai--Yin edge after
every deterministic `o(r)`-rank peel.  Therefore no theorem based only on
subexponential component density and deterministic population peeling can
work.  The gap is now between the nonconstructive slowly growing window
guaranteed by fixed-`K` compactness and this explicit occupancy scale.  A
positive result must give a quantitative smaller window or use pressure
structure beyond the operator edge.

A mere `Theta(r)` KL calculation, a copied-row spike, or a singular mixture
without a rigorous pressure direction does not count.  The conference
projected-edge benchmark is now frozen.  The live question is whether the
exact negative-disorder Gibbs law associated with contracted-temperature
child minimizers has tight conditional component `D_2`, or instead forces a
fixed mass of non-tight conditional complexity or dependence irreducible
after every explicit latent-product/row-filtration decomposition.  This must
be a theorem about that induced law, not another surrogate row fibre.

## 29. Actual-child row complexity: filtration side solved, product side open

Theorem 37.18 settles the filtration half of the former question.  For every
actual contracted-temperature minimizing child pair and every row prefix,

```math
D_2(q_\lambda(R_i\mid R_{<i})\Vert U_n)
\le\lambda^2\beta^2n/N.
```

Thus no conditional `D_2` mass escapes at fixed `beta,lambda`.  This does
not close the phase: tight autoregressive kernels need not admit a tight
latent-iid disintegration, and the archived no-gain observable is specific
to conference children rather than the actual parent `L`.

The exact replacement question is the directed row-product target excess

```math
\Delta_N=V_{\lambda,N}^{\rm row}-T_N.
```

Theorem 37.19 proves that, conditional on the full soft bridge reaching
`T_N+o(N)`, either `Delta_N^+=o(N)` and a bounded-component row-product law
constructs a linearly rare `o(N)`-accurate basin, or
`Delta_N>=eta N` on a subsequence and

```math
\inf_{p\ {\rm row\ product}}D(p\Vert q_\lambda)=\Omega(N).
```

Determine which alternative holds uniformly for actual optimizing children.
For a convergence route, the positive alternative needs a summable rate such
as `Delta_N^+=O(N^(1-delta))`.

## 30. Find a strict nonradial optimizer statistic

The full inhomogeneous contraction box (37.57) is equivalent to exact child
minimization.  Its homogeneous and fixed-size averages are only the absolute-
energy histogram, and Example 189 proves that this radial quotient does not
determine overlap geometry or even one-vertex rank-one response inside the
actual minimizer set.

Extract a statistic `S(A)` satisfying all three requirements:

1. `S(A)` follows uniformly from contracted-temperature pressure
   minimality;
2. its information content is strictly smaller than the full sign-flip or
   bridge response landscape;
3. `S(A),S(D)` control `Delta_N` above with a power-saving error, or certify
   its fixed linear lower bound.

The coordinate best-response oracle for `V_lambda^row` reconstructs the
complete bridge table on point-mass queries, so evaluating the nonconvex
product problem directly is not such a statistic.  Generic mean-field
theorems currently give only `O(N)` error; an `o(N)` theorem must exploit the
rank-one log-partition structure and genuinely nonradial optimizer data.

Theorem 37.22 identifies the first such coordinate: the squared
fixed-orientation child overlap tensor is exactly the infinitesimal cross-row
ANOVA mass.  Determine whether a uniform cluster/cumulant estimate transports
this coordinate to physical bridge amplitude `beta/sqrt(N)` and fixed
negative-disorder temperature.  Without an order-uniform remainder this is
only a tangent formula, not a solution of Question 29.

## 31. Actual-child collision influence or high-transport dependence?

The canonical row-erased product is now explicit and child-only.  Its error
is

```math
J=\log E_{r_{\rm row}^{\otimes m}}
       e^{-\lambda(\log G-E\log G)},
```

where `G` is the collision--cavity partition (37.64).  The entire hybrid
path from this product to the full escort has bounded conditional row
`D_2`, now also relative to the canonical factor by Theorem 37.27.  Its
error splits exactly into the two weighted path masses

```math
\lambda\int_0^\lambda{\operatorname {TC}(q_s)\over s^2}\,ds,
\qquad
\lambda\int_0^\lambda{\sum_iD(q_{s,i}\Vert r_i)\over s^2}\,ds.
```

Prove both are `o(N)` from an optimizer-specific child statistic strictly
smaller than the complete external-field response, ideally with a power-
saving rate.  A proof makes the canonical row product target-accurate
whenever the full negative-disorder bridge reaches the child target.

The opposite branch must be directional.  Produce fixed-alphabet,
high-transport aggregate row maps whose image of the **actual optimized
child escort** has reverse product gap `cN-o(N)`.  Low-transport maps cannot
work: Theorem 37.26 proves that any `o(N^2)` raw-coordinate support and any
bounded number of Walsh parities per row expose only `o(N)`.

Theorem 37.27 excludes two evasions.  If the canonical error is linear, a
positive density of rows has order-one scaled marginal drift or row-versus-
rest information.  If the best product repairs a linear part of that error,
a positive density of its regular factors must retune by order one.  Thus a
converse may instead certify one of those positive-density alternatives
from a high-transport child observable without solving the global product
oracle.

Finite exact data reject worst-context projective diameter as an efficient
target and show that neither total correlation nor marginal drift alone
accounts for every actual child.  They do not establish either asymptotic
branch.  Generic rank-one and fixed-overlap examples prove that weak
coordinates, bounded row Renyi complexity, central symmetry, and every
fixed overlap order are insufficient; child optimality must enter the proof.

## 32. Sector-oriented external-field stability of actual minimizing children

The actual erased-row law is now identified exactly as a sector-biased
one-vertex extension escort.  Exact minimization gives a dimension-free
density bound for its neutral version, and the path from the canonical row
product to the full escort has uniformly bounded conditional Renyi-two
components.  If its canonical error is linear, a positive density of rows
must carry marginal drift or row-versus-rest information.

The missing implication is no longer generic row regularity.  Prove, for
contracted-temperature minimizing children, one of the following genuinely
directional statements:

1. **External-field cavity stability.**  The summed one-row curvature
   `Xi_N` in (37.82) is `O(N^(1-delta))` (or has another summable density
   bound), uniformly through the sector-oriented inverse-tilt path.  This
   makes the canonical joint interaction sublinear.
2. **High-transport persistence.**  A child statistic smaller than its full
   external-field response proves that a positive density of the regular
   rows has order-one scaled mutual information or coherent factor retuning,
   with a directional reverse-product certificate.

The statistic cannot be only pressure, entropy, homogeneous/fixed-size flip
data, or any other radial transform: exact thermal minimizers with identical
such data have different actual row-path curvature.  It also cannot be an
unlabelled separate-child overlap hierarchy for a fixed orientation.  The
first viable finite coordinate is the joint `tau`-oriented four-spin tensor
in (37.83), but its small-bridge expansion is not uniform in order.  The
smallest concrete next theorem is therefore a uniform cumulant or
synchronization result transporting an oriented overlap state through the
physical external-field cavity path without reconstructing the complete
field-response table.


The canonical orientation subproblem is now separated more sharply by
Theorem 37.31.  Its exact minimal row coordinate is the sector contrast,
and its log likelihood ratio has a polynomial-size per-row quotient whose
two orientation kernels differ by `o(1)` max divergence.  Conditional on
the exact contrast, the canonical product kernel is orientation-independent.
The open object is therefore not the canonical orientation bit itself but
the joint label dynamics and the within-label full-escort residual in
(37.88).  A valid solution must do one of the following from a child state
strictly smaller than the full bridge/external-field response:

1. prove both the contrast-label image gap and the conditional residual
   are `o(N)` with a usable rate; or
2. prove that one has a fixed positive density and identify a scalable
   high-transport child observable which certifies it.

The four-real-coordinate sector--Gram state exactly controls the first
nonzero tangent and an extensive tangent forces an aggregate sector
covariance eigenmode of order `sqrt(N)`.  What remains is a uniform
physical-scale tail/synchronization theorem.  Merely applying the KL chain
rule or storing the complete contrast/scale tables does not count as a
strict reduction.


Theorem 37.32 additionally removes unbounded scalar sector bias from the
component-support side: some orientation/filtration presentation always
has uniformly bounded canonical `D_infty`, independent of both biases.
The smallest remaining question should therefore be asked in that balanced
presentation.  Does exact child minimality force its joint canonical
interaction `J` to be `o(N)`, or does a child-only high-transport
observable certify a fixed positive density of the already classified
row dependence/retuning resource?  A proof must also check that the
balanced orientation is relevant to the desired target; zero-bridge
optimality alone does not establish this.

## 33. Balanced physical cluster tightness or extensive product phase?

In the balanced actual-child presentation, Theorem 37.33 reduces one
physical promotion step to a sector--Gram coefficient and an all-order
absolute connected cross-row cumulant tail.  Decide, uniformly for exact
contracted-temperature minimizing children, whether the tangent stays below
the explicit `N^2` threshold and the tail is `o(N)` with a summable rate.

If not, produce a directional theorem: either connect fixed normalized
sector--Gram mass or linear cluster mass to `Omega(N)` reverse-product
dependence, or certify coherent positive-density regular-factor retuning.
The certificate must use less operational information than the complete
child Gibbs/external-field landscape.  In every branch, quantify whether
the balanced orientation is target-reaching or price its loss.  Merely
renaming the full high-order Gibbs law as one exact scalar does not solve
this question.

Theorem 37.34 removes the separate orientation obligation: the physical
promotion theorem now holds directly in whichever orientation reaches the
target.  Theorem 37.36 also settles the factor-representation part of the
product alternatives.  At every fixed extensive accuracy, the full optimal
row factors may be replaced uniformly by diffuse, norm-bounded,
fixed-degree carriers, and the restricted projection robustly separates
alternatives (ii) and (iii) once `J` is known to be linear.

The remaining question is therefore the following strictly narrower robust
closure problem.

> **Macroscopic finite-degree child closure.**  For one fixed carrier degree
> and norm, derive an `O(N)`-accurate rule for the restricted carrier value
> from polynomially many actual-child observables, or construct a scalable
> pair of actual minimizing children on which every such proposed rule
> collides by a fixed density.  The rule must either establish `J=o(N)` with
> a usable rate or decide the reverse-product/coherent-retuning alternative.

Exact carrier responses do not solve this: Theorem 37.37 shows that even
degree-one clipped responses are table-complete at infinite precision.
Conversely, that transform is exponentially ill-conditioned and high row
degrees are macroscopically attenuated, so the exact no-go does not falsify
an `O(N)`-precision closure.  A pressure point oracle is acceptable for
verifying a declared retuning witness but cannot be counted as a child-side
generation theorem.

Theorem 37.39 strengthens this reduction further.  The factors can be chosen
as literal nonnegative degree-`2d` Walsh densities with a fixed collision
norm, while still recovering the complete entropy-regularized product value.
The positive-part leakage and square-root entropy blockade are therefore
closed.  The operative version of the SML is now:

> **Bounded-row-degree cross-row closure.**  At one fixed row degree, evaluate
> or bound the jointly optimized square-polynomial response from an
> actual-child state which does not store the full cross-row coefficient
> tensor, to a fixed-density or power-saving accuracy; or exhibit a scalable
> actual-minimizer obstruction.  Equivalently, synchronize or truncate the
> number of participating rows while preserving their joint cancellation.

The formal bounded-row-degree tensor still has
`(sum_(a<=2d) binom(n,a))^m` entries.  Merely listing it is the old bridge
landscape in a thinner coordinate system and does not count as closure.

Theorem 37.40 supplies a minimal falsification test for any proposed answer.
Already at degree two, the spiked products expose the actual-pressure query

```math
R_{L,y}(v)=E_{\otimes_iq_{v_i,y}}L,
\qquad v\in\{+-1\}^m,
```

with an unattenuated fully active channel.  A closure theorem must therefore
prove an optimizer-specific `o(N)` range bound for this response, derive its
minimum from a genuinely smaller child algebra, or exhibit a scalable
actual-minimizer sequence with linear range.  Merely invoking low row degree
or rank-one carrier means fails this test.

The natural-channel representation sharpens this test.  With
`dmu_y=e_y dU_n` and `a=z_y/e_y`, the response is a conditional expectation
through a binary channel whose squared contraction tends to
`0.6886409151...`.  Thus the remaining issue is not one unattenuated high
coefficient but the rare minimum created by exponentially many attenuated
coefficients.  Theorem 37.44 proves that even arbitrarily strong strict
channel contraction, bounded row mutual information, and fixed literal
degree do not control this minimum in the generic carrier class.

The previously proposed all-parameter external-pressure estimate with
subgaussian proxy `o(N)` is also closed: Theorem 37.45 proves it false for
every actual optimizing-child pair, because exponentially rare rank-one
bridges have `N^(3/2)` pressure excursions.  At zero tilt, however, the same
theorem gives the exact lower-information identity

```math
{1\over mn}\sum_aE r_a^2=o(1)
\quad\Longleftrightarrow\quad
E[\tau^1\tau^2R_XR_Z]=o(1),
```

and this would imply `Var(L)=o(N)`.

Theorem 37.46 now solves the truncation half for the **entire** recovered
square-carrier class: one universal linear cap changes every admissible
product response and both optimized variational values by
`N^(3/2)e^(-Omega(N))`.  Thus the rank-one tail and all pressures above
`CN` are rigorously irrelevant to alternatives (ii)--(iii).

For the spiked subfamily, Theorem 37.47 gives a quantitative compact-tilt
replacement.  If the normalized capped-tilt cavity overlap obeys
`rho_N(S_N)=O(N^(-alpha))` through the matching window
`S_N asymp N^(alpha/2)`, then its response range is
`O(N^(1-alpha/2))`; conversely a linear range forces positive overlap at
one bounded tilt.  This is one-sided and still ranges over exponentially many
directions, so it is not yet an operational branch selector.

Theorem 37.48 removes both qualifications at the level of response bounds.
Relative entropy to the fair bridge law charges every bounded-collision
product directly, so one carrier-independent curve

```math
\rho_N(S)={1\over mn}\sup_{|s|\le S}
 E_{\Pi_{s,L\wedge CN}}\sum_a r_a^2
```

controls the pressure range of the **complete** recovered square carrier.
Decay `rho_N(S_N)=O(N^(-alpha))` on the matching window yields
`O(N^(1-alpha/2))` range, while a linear range forces positive overlap at
one bounded tilt.  This is an analytic response certificate, not yet an
operational branch decider or a way to find the optimizing carrier.

Theorem 37.49 identifies this same curve, on compact windows, with one
secant of a fair-base Gaussian-smoothing derivative.  The value at replica
number one cancels the overlap identically; the missing information is its
mixed replica/smoothing derivative.  Growing windows require an adaptive
cap `C_N=O_beta,K(1+S_N)` with quantitative replica-moment slack (or
separate plateau control).

For the optimized product value, Theorem 37.50 is sharper still.  Applying
Donsker--Varadhan duality at the objective's own parameter `-lambda` cancels
the complete row-product entropy penalty.  On the raw actual-child path

```math
{d\widehat\Pi_s\over dU}={e^{sL}\over E_Ue^{sL}},
\qquad -\lambda\le s\le0,
```

it gives, without clipping or a tail remainder,

```math
0\le E_UL-V_\lambda^{\rm row}
\le C_{LS}\lambda t^2mn\widehat\rho_N^-(\lambda),
```

where

```math
\widehat\rho_N^-(\lambda)
={1\over\lambda mn}\int_{-\lambda}^0
 E_{\widehat\Pi_s}\sum_a r_a^2\,ds.
```

More directly,
`I^leftarrow<=C lambda^2t^2mn hat rho_N^-`.
Thus small integrated overlap selects coherent retuning whenever `J` is
linear, whereas extensive reverse dependence forces positive integrated
overlap.  At `s=-lambda` this is exactly the actual negative-disorder escort.
In the smoothing coordinates it is the integral of the ordinary **raw**
secant `(widehat A_N-1)/(s-1)` on `[-lambda,0]`; the interval avoids replica
one and negative tilts suppress the high-pressure tail, so no cap, mixed
derivative, adaptive window, or carrier optimizer is needed.

The previous SML was:

> **`L_raw-negative-overlap`.**  For actual contracted-temperature
> minimizing children in target-reaching orientations, prove
> `hat rho_N^-(lambda)=o(1)`, with a power/summable rate for the basin
> recurrence; alternatively prove that positive integrated overlap produces
> an explicit favorable reverse-product direction, or construct an actual-
> minimizer sequence showing the obstruction.

Theorem 37.52 now resolves the third alternative and closes the decay route
in the strong-channel regime.  If `mn/N^2>=gamma_0` and
`beta^2 gamma_0>2log2`, then every actual child pair satisfies

```math
\liminf_N\widehat\rho_N^-(\lambda)>0.
```

For balanced splits this begins at `beta>sqrt(8log2)`.  The obstruction is
the exact `2^(O(N))` rank-one support of the actual latent child channel in
`Theta(N^2)` bridge coordinates.  Since an unbounded sequence of fixed
`beta` values is required to squeeze the ground-state problem, raw overlap
decay cannot be the Level-6 mechanism.

Theorem 37.53 shows that this floor is not automatically directional: a
rank-one channel with one projective factor frozen has positive raw overlap
at every fixed `beta,lambda`, while every negative tilt is exactly row
product and `I^leftarrow=0`.  Hence support entropy, positive raw overlap,
and actual rank-one algebra alone cannot decide reverse dependence.

Theorem 37.54 supplies the correct inverse-order diagnostic but also fixes
its scope.  For the canonical row product `r` and actual escort `q`,

```math
J=D(r\Vert q)
=\int_0^\infty\mathfrak K_r(P_t^r(dq/dr))dt.
```

This lifetime integral vanishes on the row-factor sharpness channel and
charges pure parity at its true `Theta(N^-1)` scale.  It is nevertheless
exactly `J`, not yet a strict optimizer-specific reduction, and it cannot
by itself distinguish `I^leftarrow` from `J-I^leftarrow`.

The replacement SML is therefore:

> **`L_row-lifetime-closure`.**  From actual-child data strictly coarser than
> the complete bridge density, prove either a summable `J=o(N)` bound for the
> row-refresh lifetime, or split its extensive mass into (i) an explicit
> coherent factor-retuning direction and (ii) a cross-row component that
> lower-bounds `I^leftarrow`.  The split may not invoke the unknown optimal
> product or store the full smoothed density `P_t^r(dq/dr)`.  This resolves
> only the product phase; a Level-6 recurrence still needs a separate
> target-reach/gain theorem in the correct inequality direction.

Theorem 37.51 fixes the scope of the alternative.  At the best reverse row
product, the exact centered score is a row-ANOVA interaction residual and
`I^leftarrow` is its negative cumulant.  This centering requires the unknown
product optimizer.  Moreover, a full-parity landscape has unit normalized
centered edge mass but `I^leftarrow=Theta(N^(-1))`.  Hence a positive-overlap
converse must use actual-child structure that excludes high-row-order
concentration or controls the centered negative tail; another quadratic
centering cannot decide the branch.

The quartic coefficient already has `O(sqrt(N))` physical range, so any
obstruction must be nonperturbative or begin at higher auxiliary order.
Existing child minimality gives only one-replica `O(N)` endpoint bounds on
internal edges and cannot prove the required replicated cross-edge
derivative rigidity directly.  Another raw MGF, carrier net, or singular-
value estimate is no longer an admissible route.

### Actual-minimizer closure after the common-sign retuning no-go

Theorem 37.58 shows that rank-one latent support and bounded conditional
row Renyi complexity allow

```math
J=\Theta(N),\qquad I^\leftarrow=o(N),\qquad
J-I^\leftarrow=\Theta(N),
```

with `Omega(N)` row-lifetime mass already on a fixed bounded time interval.
The example is generic, not an actual minimizing-child sequence. The next
question must therefore use an optimizer identity:

> Do contracted-temperature minimizing children exclude the two-word
> shared-latent retuning mechanism, or can one construct from their
> one-child data an explicit product `s_N` satisfying
> `D(s_N||q)<=I^leftarrow+o(N)` with a summable rate?

Using the unknown optimal product `p^*` is not an answer. Nor is another
bound based only on support cardinality, conditional row `D_2`, or
bounded-time row-noise regularity: all three hold in the counterexample.
Even a positive answer resolves only the product branch; target reach
remains separate.

Theorems 37.55 and 37.57 also close the scalar effective-entropy variant of
the support argument.  Weighted-max transport can interpolate from
`log |supp mu|` to prior Renyi entropy, but every nonzero weight on the prior
introduces

```math
\mathcal C_{q_s}
=E_{\bar\mu_{q_s}}[-\log\mu]
=H(\bar\mu_{q_s})+D(\bar\mu_{q_s}\Vert\mu).
```

Thus an improved effective-support threshold requires a sublinear theorem
for the negative-path posterior retuning.  If the retuned posterior puts
fixed mass on an exponentially rare child set, that divergence is already
linear.  Prior Shannon entropy and positive-temperature full support do not
decide between these cases.  This is a scope closure, not a replacement for
`L_row-lifetime-closure`.

Theorem 37.61 now rules out the literal common-sign mechanism and its
low-rate narrow-cluster variants for actual children at sublinear retuning
cost.  Scalar pressure contraction bounds `t max|H_A|/m`, while the exact
quadratic Hamming-sphere identity turns this into a uniform sector
min-entropy rate

```math
H_\infty(\mu_{A,s})\ge(\eta_\beta-o(1))m,
\qquad \eta_\beta\ge e^{-\beta^2}/16.
```

Thus the actual rank-one prior has maximum atom and collision probability
`exp{-Omega_beta(N)}`.  A posterior placing fixed mass on any latent
catalogue of rate below `eta_beta` pays linear KL; the same holds for a
union of narrow rank-one overlap caps below the explicit entropy-rate
threshold in (37.220).

The strict remaining question is no longer whether one or finitely many
common-sign words can carry the row lifetime.  It is:

> Can the actual inverse-disorder posterior retune diffusely across an
> exponential-rate family and create `J-I^leftarrow=Omega(N)`?  Either
> identify an optimizer-specific low-dimensional direction of such diffuse
> retuning, or prove that its contribution to `J-I^leftarrow` is `o(N)`
> with a summable rate.

Theorem 37.62 splits this question exactly.  The rank-one prior is uniform
on each combined-energy shell
`H_A(x)+epsilon H_D(y)=e`, and there are only `O(N^2)` shells.  Every actual
latent retuning KL equals the KL of that polynomial-state shell distribution
plus the conditional entropy deficit inside the uniform shells.  Hence the
next proof should first test whether the shell KL controls a coherent
low-complexity product direction.  If it does not, the irreducible remainder
is now specifically diffuse exponential-rate **within-shell** selection,
not an arbitrary full-table retuning.

Theorem 37.63 prevents treating that split as a generic row-lifetime proxy.
A partially shared row-sign group has global exponential min-entropy,
collision decay, and narrow-cap decay, and its averaged posterior equals its
prior under every disorder tilt. Thus both shell terms vanish, while a
positive-density common block still has `J-I^leftarrow=Theta(N)`. The fixed
right factor makes this a non-actual channel. Accordingly the next lemma
must use factorwise actual-child spread or another optimizer identity; a
bound in terms of global prior spread and averaged-posterior shell KL alone
is false.

Theorem 37.64 supplies the factorwise spread that actual children really
have.  For every prescribed positive-density child block `U`, even after
conditioning on every exterior child spin, its sector law has a uniform
positive min-entropy rate.  Thus the partial-common-row obstruction is not
merely absent globally: an actual child cannot hide a frozen or low-rate
common phase on any macroscopic coordinate block.

This changes the smallest missing lemma but does not close it:

> **`L_diffuse-factor-retuning`.**  Under the exact actual-child Gibbs law
> and its uniform macroscopic conditional min-entropy bounds, prove that
> extensive canonical excess `J-I^leftarrow` either (a) is impossible up to
> a summable error or (b) yields one explicit low-information coherent
> direction.  Equivalently, decide whether an exponentially diffuse common
> phase, with no positive-density frozen coordinate block, can carry
> `Theta(N)` row-lifetime mass under the actual inverse-disorder tilt.

The Hamming-spread estimate alone is not an independence theorem: it does
not bound total correlation, overlaps at moderate radius, or negative-path
posterior retuning.  A valid proof must use an additional optimizer-specific
identity and cannot invoke the complete latent table or the unknown best
reverse product.  Even success would still require a separate target-reach
theorem before producing a Level-6 recurrence.

Theorem 37.65 shows that the word "factorwise" must retain its local
meaning.  Making both complete factors exponentially diffuse does not help:
a double-partial product group still has `J-I^leftarrow=Theta(N)` and zero
averaged-posterior retuning.  It is excluded only because each factor hides
a positive-density two-point marginal.  Consequently the sharp generic
falsifier for `L_diffuse-factor-retuning` is a genuinely noisy common phase
which obeys the conditional entropy bound on **every** macroscopic block.
Either such a channel retains linear coherent regret--proving that spread
alone is still insufficient--or its cancellation identifies the first
useful diffuse-synchronization theorem.

Theorem 37.66 resolves that falsifier in the negative: BSC smoothing of the
common phase satisfies exponential conditional min-entropy on **every**
factor subset and still has `J-I^leftarrow>=cN`.  Therefore even the exact
conclusion of Theorem 37.64 is not, by itself, a row-lifetime closure
principle.  The example also supplies the positive clue: its retuned phase
is captured exactly by the polynomial-state pair of factor magnetization
counts.

The smallest missing lemma is consequently sharper:

> **`L_actual-diffuse-phase-quotient`.**  For the actual optimizing-child
> inverse-disorder posterior, either construct an optimizer-specific
> quotient of subexponential (preferably polynomial) state complexity whose
> retuning controls `J-I^leftarrow` up to a summable error, or prove that no
> extensive diffuse phase retuning occurs.  The quotient must be defined
> from child data without the full Gibbs/bridge table and must include a
> theorem converting its displacement into one coherent row-product
> direction.

Factor magnetization counts solve the generic exchangeable example only
because it has `S_m times S_n` symmetry.  Arbitrary minimizing signings have
no such symmetry, so importing those counts verbatim is not a candidate.
The combined-energy shell is polynomial but Theorem 37.63 shows that it
does not generically control lifetime.  A valid actual quotient must use a
new optimizer-specific synchronization or rare-event/renormalization
identity.  Target reach remains separate even if this lemma is proved.

Theorem 37.60 records the information boundary.  Values on all subsets of
an `m`-edge cycle basis recover the entire augmented Gibbs table by an
invertible transform, so that route is not a compression.  The inequality
directions alone admit abstract half-atom laws.  The successful spread
input is genuinely coarser: one scalar pressure bound plus the exact
Hamming-sphere algebra of quadratic forms.  No target-reach or Level-6
recurrence follows from spread alone.

### Rooted symmetry and overlap regularization: the next exact boundary

Theorem 37.68 gives the first finite actual-child saturation test.  At the
two order-eight minimizing classes, the negative posterior uses every one
of the `19` or `22` simultaneous signed-similarity/rooted-profile cells,
while the polynomial combined-energy shell is strictly too coarse.  Hence
an all-order orbit proposal must prove a genuinely asymptotic statement:
subexponentially many cells, `o(N)` symmetry-breaking KL for an approximate
action, or response coalescence across exponentially many exact cells.
Finite symmetry by itself is not evidence for such a theorem.

Theorem 37.69 supplies a distinct, literature-backed regularization.  An
`o(N)` generic Gaussian perturbation of every deterministic child forces an
ultrametric limiting replica Gram law and changes the normalized bridge
response uniformly by `o(N)`.  But the Gram law forgets the physical
coordinate embedding, and the available synchronization theorem either
fails for two separate perturbations or destroys the product law through a
joint perturbation.

This sharpens `L_actual-diffuse-phase-quotient` to one admissible form:

> **`L_inverse-escort-coordinate-lifting`.**  For separately GG-regularized
> actual optimizing children, augment the discretized overlap laws by
> subexponential coordinate data and approximate the inverse-escort-
> typical row/cavity responses with total `o(N)` error.  The construction
> must synchronize the two child trees without using a joint perturbation
> that changes their product, and it must convert the retained displacement
> directionally to `J-I^leftarrow`.

Uniform reconstruction over every bridge is suspected to retain the whole
latent law and is not the target.  The quantifier must be the actual
negative escort, with an explicit escaping-mass bound.  Proving such a
lifting theorem, or a linear information lower bound against every such
coordinate augmentation, is the single rare-event/renormalization fork.
