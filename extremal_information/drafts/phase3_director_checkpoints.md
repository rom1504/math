# Phase 3 research-director checkpoints

This file records steering judgments, not additional theorem claims.  A
checkpoint changes the next experiment; it does not end the research run.

## Checkpoint 1: unrestricted syndrome responses

**Question.** Did the macroscopic syndrome target reveal a general law, or
only another code-specific dynamic program?

It produced two positive compression mechanisms and a stronger obstruction.

1. A word-length profile `lambda_S` is one-Lipschitz in the Hamming chart of
   any basis contained in `S`.  Storing its values on a radius-`r` covering
   code determines every adversarial appended-support response to error `r`,
   because min-plus convolution and the final maximum are sup-norm
   nonexpansive.  At `r=1` this uses

   ```math
   O(2^w\log w/w+w^2)=o(2^w)
   ```

   bits.  This answers the requested full-context quotient question much
   more strongly than a fixed-relative-error net, although composing two
   independently approximated sketches can add their errors.

2. Thresholding supplies a separate, exactly closed algebra.  Supports of
   radius below `R` form an absorbing ideal under union.  Collapse that ideal
   and retain every radius-at-least-`R` support exactly.  Kneser's theorem
   implies that such a retained support has at most
   `2|F_2^w|/(R+1)` elements.  Hence the quotient uses

   ```math
   O(2^w\log R/R)
   ```

   bits, composes with no accumulated error, and has radius error at most
   `R/2-1`.  Taking `R -> infinity` with `R=o(w)` gives the stronger target:
   an `exp(o(2^w))`-state exact composition algebra with `o(w)` response
   error.  The algebraic operation is the classical Rees quotient; the new
   mathematical content is the high-radius sparsity theorem and its response-
   compression consequence.

3. These upper bounds do not make the response class low-dimensional.  A
   Grassmannian family of dense-carrier supports has pairwise future-response
   distance controlled by subspace injection distance.  Constant-dimension
   packing yields an independently audited lower bound

   ```math
   ((1/2-2\epsilon)^2-o(1))w^2
   ```

   bits at error `epsilon*w` for every fixed `epsilon<1/4` (with the displayed
   coefficient understood as a supremum over a strict packing margin).  This
   is a moving global-carrier obstruction, not the earlier fixed direct-sum
   source.  It does not approach the exponential-in-`w` upper description.

The selected dichotomy is therefore resolved positively, while its optimal
rate remains open.  The useful abstraction is not “store the dynamic-program
table.”  Two different laws appeared:

- **contractive profile sampling:** regularity on a low-covering-entropy
  interface plus nonexpansive continuation;
- **hard-core ideal compression:** collapse an absorbing easy ideal and count
  the surviving difficult objects.

Both mechanisms are candidates for a general future-response law.  Neither
by itself explains when local information becomes macroscopically exposed.

### Immediate next theorem

The next phase will test generality rather than optimize the syndrome rate.
It must do both of the following:

1. formulate the hard-core quotient as a response theorem and validate it on
   a non-code model with a strict quotient; and
2. determine whether selective neutralization admits a non-tautological
   packing theorem that explains the block and Grassmannian lower bounds and
   applies to a second model.

If these reduce only to the definitions of a Rees quotient and a metric
packing, the framework has learned a code theorem but not a general law; the
next checkpoint must pivot to query-mass-sensitive posterior geometry.

## Checkpoint 2: terminal ideals, congruence entropy, and a third model

**Question.** Did the generalization explain feature-algebra growth, or only
rename dynamic programming and semigroup quotients?

The classical algebra was not new, but three quantitative consequences passed
independent audit.

1. For any scalar response monoid with an absorbing terminal state, a future-
   response ball about that state is an ideal and hence has a Rees quotient.
   For a nonnegative antitone deficit `F`, the result is sharp: at uniform
   error `eta`, the largest possible summary cell containing the terminal
   state is exactly

   ```math
   \{x:F(x)\le2\eta\},
   ```

   and midpoint decoding makes this whole cell safe under arbitrarily many
   future compositions.  This is a genuine stopping rule for feature growth,
   although the size of the remaining hard core is model-specific.

2. Matroid residual rank supplies a genuinely different third model.  Its
   exact all-future quotient is the flat lattice, and its response metric is

   ```math
   d_F(X,Y)=\max\{r(X\vee Y)-r(X),r(X\vee Y)-r(Y)\}.
   ```

   For the binary projective matroid this becomes subspace injection distance.
   Grassmann packing proves `Theta(w^2)` response bits at every fixed
   distortion `epsilon*w`, `epsilon<1/4`, matching the basis-description upper
   bound.  Thus the framework predicts a sharp macroscopic rate rather than
   merely observing that matroid closure exists.

3. A one-shot response net need not be a reusable feature algebra.  Exact
   closed scalar summaries are characterized by monoid congruences whose
   classes have `F`-oscillation at most `2epsilon`.  On the prime cycle

   ```math
   F_p(x)=\cos(2\pi x/p),
   ```

   all translated response maps have an `O(1/epsilon)` uniform net, but for
   fixed `epsilon<1` every exact closed algebra has all `p` reachable states
   once `p` is large.  Smooth low-dimensional response geometry therefore
   does not imply low-index algebraic closure.

The bounded-multi-carrier attempt also ended decisively.  Replacing `m`
dense projective carriers by their span changes every future sparse-synthesis
response by at most `m-1`, even under an arbitrary background and future
dictionary.  Hence `o(w)` carriers collapse to one at `o(w)` normalized error;
flags and bounded-carrier variants cannot improve the quadratic packing
scale.  A stronger syndrome lower bound needs nonlinear carriers, rooted
incidence/multiplicity, or exponentially many selectively addressable
carriers whose union does not collapse to its span.

### Director judgment

The framework is becoming unified around a real distinction:

- response-metric entropy prices a stored system queried by a raw future;
- congruence entropy prices an exact summary-only algebra of unbounded depth;
- absorbing response ideals are a structural mechanism making the two
  compatible without error accumulation;
- selector/carrier geometry is a structural mechanism forcing information
  growth before approximation is paid.

This is more than convex duality, but the congruence characterization itself
is classical universal algebra/Myhill--Nerode logic.  Its value is the new
separation and the model-specific rate theorems.

### Immediate continuation

The next test is a bounded-depth interpolation: determine whether an actual
`delta`-net can be re-rounded through `ell` compositions with a universal
error law, and whether a convex or tropical state can beat that accumulation
without secretly defining a congruence.  This checkpoint does not end the
campaign.

## Checkpoint 3: algebraic reuse and structural hard-core entropy

**Question.** Can the metric/congruence distinction be sharpened into a law
that predicts when approximation survives repeated composition?

Three independently audited results answer substantial parts of that
question.

1. For a fixed-element blur `P_b(x)=x star b` in any translation-contractive
   metric monoid, the exact `m`-factor algebra defect is

   ```math
   d(b^{star m},b).
   ```

   Thus repeated local approximation is governed by one power orbit, not by
   a generic sum of local errors.  Idempotent `b` gives an exact retraction;
   an external bounded-repair metric plus a one-sided Lipschitz response then
   gives an arbitrary-depth midpoint error paid only once.  Tropical subgroup
   profiles and fixed-flat matroid contraction satisfy this certificate,
   while Minkowski uncertainty has a linearly escaping power orbit.  This is
   a constructive sufficient mechanism, not a classification of all useful
   congruences.

2. Projective residual rank has a non-Rees exact quotient.  Projection modulo
   a `d`-dimensional subspace is a join homomorphism whose every fiber has
   exact response width `d`; it uses

   ```math
   ((w-d)^2/4+o(w^2)) log_2 q
   ```

   bits.  More generally every join congruence factors canonically into such
   a linear kernel followed by a zero-separating congruence.  Exhaustive
   enumeration of all `3,616` congruences of `L(F_2^3)` verified the
   decomposition and oscillation formula.  This supplies a third strict
   composable model rather than another coding relabeling.

3. The syndrome hard core is exponentially thinner than the Kneser estimate.
   A diameter-`D` Cayley geodesic forces every affine fiber of the support to
   be a binary diameter-two anticode.  For `D>=3` each fiber has at most
   `D+1` points, and counting the possible fibers gives

   ```math
   log_2 #\{S:D(S)>=R\}=O(w^2+w2^{w-R}).
   ```

   The corresponding Rees quotient is exactly closed.  At error `epsilon w`
   it uses

   ```math
   2^{(1-2epsilon)w+O(log w)}+O(w^2)
   ```

   bits, with no accumulation under arbitrarily many unions.  Two independent
   exhaustive checks covered all spanning supports through width four; the
   `D=2` anticode exception is explicit and excluded.

The emerging general law has two independent prices.  **Algebraic reuse** is
controlled by congruence or saturation dynamics (exactly by the power orbit
for principal blurs).  **Information size** is controlled by geometric
repair and the structural entropy of the response hard core.  Neither price
alone implies the other.  Prime cycles have tiny metric entropy but rigid
congruence entropy; convex blurs have bounded one-step geometry but escaping
reuse defect; syndrome supports have a large ambient source but a thin
high-response core.

### Immediate continuation

This checkpoint is a steering event, not a stopping point.  The geodesic
proof exposes a natural finite linkage condition rather than an assumed
overlap synchronization hypothesis: every projected-zero generator cycle
must have no more transverse Hamming mass than its length.  The next theorem
will determine whether that cycle condition synchronizes all affine-fiber
labels to a linear section, with a quantitative zero-temperature error and a
sharp counterexample if the constant cannot be improved.  That test decides
whether the hard-core count is merely enumerative or also reveals a general
mechanism by which composition constrains latent information.
