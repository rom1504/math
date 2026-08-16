# Minimal open questions

These are theorem targets selected after the carrier-capacity checkpoint, not
a catalogue of neighboring ideas.  The unrestricted syndrome dichotomy is no
longer the lead question: a composition-stable subexponential response net and
a quadratic-bit macroscopic packing are both proved.  Optimizing that gap is
secondary to understanding the new carrier law.

## 1. Resolve the Hamming Grassmannian gap

For mixed-relation models, composition creates a homomorphism

```math
h:K=Z/Z_{\rm loc}\longrightarrow W
```

and a presented carrier `C_h=h(K)`.  Theorem 12.3 reduces response complexity
to carrier metric entropy.  Theorem 13.3 identifies the optimal synchronizing
quotient rank with anticode codimension, while Theorem 13.4 proves that binary
Hamming separated rank lies a linear amount below it.  Universal scale-rank
duality is therefore false.

The remaining object is the actual Grassmannian packing

```math
P_{D,k}(\Delta)=
\log_2\operatorname{Pack}
\bigl(\operatorname{Gr}_k(F_2^D),d_H,\Delta\bigr).
```

Determine whether large families of subspace carriers not contained in one
common separated host can fill the code--anticode gap.  Does the exponent
equal the common-host lower bound, the quotient-state upper bound, or a
genuinely intermediate invariant?

**Success:** matching exponential upper and lower bounds for
`P_(D,k)(Delta)` in one nontrivial linear regime, transferred through the
presentation toll to a response rate--distortion theorem.

**Useful negative result:** prove that common-host packings are exponentially
suboptimal, or that anticode-quotient summaries retain exponentially too many
states, thereby isolating a new middle invariant.

This is the single most important next theorem.  The two-scale and rank-
metric models have zero code--anticode gap; binary Hamming space is the first
rigorous setting where the general theory leaves a leading interval.

## 2. A natural strict metric synchronization

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

## 3. Query-mass-sensitive carrier capacity

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

## 4. Effective-rank composition beyond fixed dimension

Shapley--Folkman bounds the loss from replacing a nonconvex reachable sum by
its convex roof by the `r` largest component diameters, where `r` is the
effective affine-difference rank.  Determine whether **exposed** or
query-dependent dimension can replace `r`.

**Success:** a bound controlled by the dimension or Gaussian width of faces
actually reached by a declared Lipschitz query, with subextensive error in a
model where ambient rank grows.

**Falsifier:** a same-roof pair with low exposed complexity for every local
component but a leading composed response gap.

## 5. Constrained compactness and realization

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
