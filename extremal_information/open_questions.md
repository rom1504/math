# Minimal open questions

These are theorem targets selected after the carrier-capacity checkpoint, not
a catalogue of neighboring ideas.  The unrestricted syndrome dichotomy is no
longer the lead question: a composition-stable subexponential response net and
a quadratic-bit macroscopic packing are both proved.  Optimizing that gap is
secondary to understanding the new carrier law.

## 1. Robust converse to finite tropical lumpability

The forward structural problem now has two rigorous positive criteria. Exact contextual
refinement stabilizes if and only if a finite deterministic path-realizing
quotient exists. More usefully, a finite oriented affine arrangement closed
under generator pullback is a checkable sufficient certificate; for selector
coordinate differences, finiteness of the unrestricted pullback closure is
equivalent to zero weight on every reachable unit-transport cycle. An
approximate block quotient plus uniform switching contraction gives a finite
depth-independent response simulator without enumerating active cells.

The remaining problem is a **robust converse**. Failure of one affine
certificate need not kill every quotient: formal selector branches may be
dynamically incompatible, and saturation may eliminate their infinite
normal orbit. Nor must failure produce amplitude drift: irrational rotation
has zero drift but an unbounded exact contextual partition. The compact
Cantor encoder--decoder instead has exponential approximate future memory.

For a compact rational selector-affine system and fixed response accuracy
`epsilon`, characterize which of the following is forced:

1. a finite approximate bisimulation/congruence;
2. a realizable positive response cocycle, giving extensive drift; or
3. an explicit packing of the horizon response trees whose cardinality tends
   to infinity.

The desired theorem must give a checkable witness, not merely define the
minimal contextual pseudometric. A promising restricted case is a finite
unit-transport normal graph with polyhedral enabling domains: decide whether
a nonzero syntactic holonomy cycle has a repeatable feasible subcycle, and if
not, bound the pullback-refinement entropy of the transient language.

**Success:** a decidable dichotomy/trichotomy for one nontrivial selector
class, with either a finite quotient of controlled size or a quantitative
realizable drift/response-packing witness.

**Stop condition:** common normal directions without affine offsets, pairwise
face compatibility, or an unconstrained orbit-language enumeration does not
advance the target.

In parallel, the bounded-presentation benchmark leaves a sharp scale gap.
Binary `m`-parameter max-affine grammars have an `O(m^2)` coarse entropy
upper bound, while robust facet deletion gives `Omega(m log m)`.  A robust
common-normal-fan charging theorem, or an `Omega(m^2)` separated
multi-polytope construction, would settle which exponent is real.

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

## 5. Effective-rank composition beyond fixed dimension

Shapley--Folkman bounds the loss from replacing a nonconvex reachable sum by
its convex roof by the `r` largest component diameters, where `r` is the
effective affine-difference rank.  Determine whether **exposed** or
query-dependent dimension can replace `r`.

**Success:** a bound controlled by the dimension or Gaussian width of faces
actually reached by a declared Lipschitz query, with subextensive error in a
model where ambient rank grows.

**Falsifier:** a same-roof pair with low exposed complexity for every local
component but a leading composed response gap.

## 6. Constrained compactness and realization

Fixed-interface response bodies have unrestricted finite recovery sequences.
Characterize which limiting response or presented-carrier states are realized
at all large sizes inside a constrained family such as linear codes, dense
CSPs, or bounded-width factor graphs.

**Success:** a Gamma-limsup/recovery theorem preserving the declared response
with vanishing normalized loss and without storing a target optimizer.

**Stop condition:** finite-state approximation outside the constrained model
does not answer this question.

## Reconnection rule

Do not return to the motivating signing problem.  Reconnection requires a
carrier or synchronized quotient that arises naturally there, has controlled
growing-interface entropy, closes under the relevant composition, and has a
finite realization theorem.  None is currently proved.
