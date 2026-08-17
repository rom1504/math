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
`o(N_j^(3/2))`.  Seek a balanced-endpoint compiler whose whole accessible
language remains at target scale, or prove that every sublinear compiler
exposing the periodic field necessarily admits such a super-target endpoint.

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
