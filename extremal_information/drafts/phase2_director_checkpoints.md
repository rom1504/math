# Phase 2 research-director checkpoints

This file records selection judgments, not additional theorem claims.  The
surface files are updated only after the cited claims pass independent audit.

## Checkpoint 1: response geometry, separator kernels, and outer spectra

**Question.** Has this produced a concept or theorem that explains something
we could not previously formulate, or are we merely repackaging convex
duality, dynamic programming, or information theory?

### Survives as potentially generative

The response-separation polytope `Gamma(R)` and its scalar inverse-Hamming
modulus have yielded a sharp posterior-width theorem.  Its entropy step is
classical binary rate--distortion theory, but the exact identity

```math
\inf_\pi
{\operatorname{Var}_\pi(R)\over
 \sum_e\operatorname{Var}_\pi(A_e)}
={\kappa(R)\over4}
```

turns nonlinear extremal-response geometry into a continuous information
price.  It applies with the same proof to a shifted dense Ising response and
to rooted nearest-code distance, and it is sharp for a hypercube response.
This crosses model classes and strengthens the earlier threshold/packing
argument.  The subsequent audit confirmed the constants and sharpness, so it
is **Level 3** at the project level.

### Retained as explanatory baselines

- Boundary response kernels, max-plus gluing, and separator profiles give a
  correct exact quotient, isometry, and universal `Q^2`-bit lower bound.  This
  is a useful response-theoretic unification of transfer matrices and
  variable elimination, but its current mathematical mechanism is classical.
  It is **Level 2**, not a new theory by itself.
- The outer distance polynomial is the exact quotient for the complete
  root-averaged pressure curve and has a product algebra.  It is classical
  coset-layer/generating-function structure.  Covering radius alone is a
  smaller sufficient scalar for Cartesian products, so the outer spectrum is
  not the sought minimal rooted augmentation.  It remains a **Level-2
  taxonomy example** and a strict collision with inner replica data.

### Next discriminating tests

1. Apply the posterior-width theorem to the endpoint-pinning experiment for
   universal boundary kernels.  A correct continuous `Q^2` rate curve would
   show that the information object predicts feature-algebra cost rather than
   merely redescribing a prior proof.
2. Search inside a restricted realizable kernel semigroup (trellis or code)
   for response-metric entropy smaller than the universal matrix cube.
3. Audit the proposed deterministic synchronization theorem.  It counts as
   new only if it supplies uniform zero-temperature control and a
   counterexample proves that ordinary average synchronization cannot replace
   its cross-root condition.

## Checkpoint 2: deterministic synchronization

**Question.** Is the proposed theorem a real deterministic compression
mechanism, or does its linkage hypothesis simply assume the conclusion?

The draft separates two obligations:

```math
\text{mixture ultrametricity (local no crossing)}
+\text{scalar cross-root linkage}
\Longrightarrow
\text{uniform species synchronization}.
```

The linkage condition is expressed only through the retained total-overlap
labels and path backtracking; it does not mention the hidden species values.
It is therefore logically weaker than synchronization.  The audited bound is

```math
\max_e|R_s(e)-L_s(q(e))|
\le {\tau+3D\eta\over\lambda_s},
```

with the same error transferred to every Lipschitz zero-temperature coupling
query.  A perfect-matching construction satisfies PSD,
exchangeability, ultrametricity of every nonnegative mixture, and vanishing
conditional variance, yet retains a fixed conditioned extremal gap.  Both
the theorem and falsifier survived adversarial audit.  This is **Level 3 on a
limited deterministic hierarchical class**: it explains when a fibre
response body collapses to a scalar and why averaged replica criteria can
miss the obstruction.  It is not yet a deterministic Parisi theory.

The route stops if the linkage condition is found equivalent to inspecting
the hidden profile or if its only closing class assumes species homogeneity
in disguise.

## Checkpoint 3: composition stress test and known-model validation

**Question.** Have the new objects generated mathematics beyond the convex
roof, or are we accumulating alternate descriptions of dynamic programming?

### Strict advances confirmed by final audits

1. **Syndrome-rooted code algebra.**  For binary parity-check fragments over
   a fixed syndrome group, the coset-leader profile composes by min-plus group
   convolution and, because equal column types cancel in pairs, reduces to
   set union of distinct nonzero syndrome types.  Special future fragments
   expose every support bit.  This gives a strict quotient of the code with a
   closed repeated-composition algebra and sharp `Theta(2^w)` response-bit
   complexity.  It is the requested second nontrivial model beyond
   fixed-rank mean field.  An equal-outer-spectrum pair separated by one fixed
   environment proves that the group labels, not only the distance histogram,
   are necessary.
2. **Robust tropical response rank.**  A four-cell crossing gap `G` forces
   approximate min-plus factor rank at least the fooling-set size under
   uniform error below `G/4`.  Sheshadri's code transversal has `G=2`, so the
   every sub-half-unit uniform approximant still needs at least `2^s`
   channels, and the exact table attains that count.  This is a
   lattice-scale stable feature-growth theorem inside a
   structured code class.  The exact rank theorem was already in the ledger;
   the robust corollary is the new deduction.
3. **Fixed-dimensional amortized nonconvexity.**  Shapley--Folkman implies
   that the Minkowski sum of component response sets lies within the sum of
   the `p` largest component diameters of its convex roof, independent of the
   number of factors.  Hence the composable convex response body answers every
   Lipschitz nonlinear aggregate query with nonaccumulating error in fixed
   feature dimension.  A same-zonotope vector-balancing pair has discrepancy
   gap `2d` when the dimension grows, locating the sharp structural boundary.

These are mechanisms, not merely vocabulary: algebraic cancellation,
four-cell tropical separation, and amortized nonconvexity each prove a new
operational conclusion.  The underlying coding, tropical, and
Shapley--Folkman tools are classical or recently imported and must be cited
as such.

### Concepts narrowed or rejected

- `Gamma(R)` tensorizes exactly only when response channels are retained
  orthogonally.  Identical child polytopes can yield either full separation or
  zero separation after same-space addition, and universal max-plus gluing
  makes `Gamma` collapse after two factors.  It remains a sharp
  rate--distortion certificate for a fixed response embedding, **not** a
  compositional state.
- Exact tropical rank alone gives no average-error theorem.  The zero-diagonal,
  one-off-diagonal matrix has rank `r` and uniform robustness below `1/2`, yet
  a rank-one matrix has normalized mean-square error `1/r`.  Any averaged
  theorem needs a mass/exposure condition.
- REM/GREM covariance is an ensemble response-law parameter, not a quotient
  of a fixed quenched realization.  Gaussian regeneration, not deterministic
  landscape compression, explains its small state.

### Director selection

Subject to independent proof audits, the framework has become generative in
a limited but genuine sense.  The strongest common next theorem is not a
larger roof.  It is a **query-mass-sensitive response complexity theorem**
linking posterior width and tropical witness geometry.  It must interpolate
between uniform robustness and average distortion, and Proposition TR.3
shows that an explicit anti-rare-exposure hypothesis is unavoidable.

The deterministic-synchronization route remains a ranked alternative: derive
cross-root linkage from finite replica identities plus an anti-rare-face
condition.  No attempt should reconnect to the original signing problem
during this theory campaign.

## Checkpoint 4: weighted exposure and severe selection

**Question.** Did the query-mass-sensitive tropical theorem create a
macroscopic response-complexity law, or only diagnose why the uniform theorem
does not transfer to average loss?

The patched weighted four-cell theorem is rigorous when the distinguished
anchor rows and columns are separately injective.  It gives a clean finite
lower bound by the maximum exposed monochromatic matching left by every
`k`-channel coloring.  Removing injectivity is false; repeated witness
rectangles double-count the same error cells.

Its canonical code application is negative.  A Sheshadri transversal of size
`q=2^s` has exposure at most `1/(2q)` under uniform state-pair queries.  For
the graph code `C_t={(z,z)}`, the exact one-channel certificate is `1/(8q)`,
and the rank-one constant approximation to the normalized distance table has
mean-square error `1/(16t)`.  Hence exponential exact/lattice-robust tropical
rank does not imply macroscopic diffuse rate--distortion.  More pairwise
witness graphs on this transversal are stopped.

The severe audit therefore classifies the program as **Level 3 locally and
Level 2 globally**.  The posterior-width inequality, syndrome-support future-
query quotient, and finite monotone-linkage synchronization theorem are the
three core project-level results.  They are mutually compatible but do not
share one compositional state.  Shapley--Folkman is strong imported
validation; weighted exposure remains a Level-2 certificate.

The selected next theorem is the macroscopic response-metric entropy of
syndrome supports at additive distortion `epsilon*w`.  It has a strict
two-checkpoint test: either find a composition-stable subexponential quotient
with `o(w)` error, or an `exp(Omega(2^w))` realizable packing separated by
`Omega(w)`.  Stop if available packings remain only lattice-separated and
approximate quotients accumulate linear error.  The project does not
reconnect to the motivating signing problem at this checkpoint.

## Checkpoint 5: first macroscopic syndrome-response theorem

**Question.** Did the selected normalized target produce a new information
law, or merely repeat the individual support-bit packing at a larger scale?

It produced a new joint-query mechanism.  Split `F_2^w` into fixed
`L`-dimensional summands.  In each block a latent bit selects either a basis
(radius `L`) or all nonzero generators (radius one).  A legal future fragment
selects an arbitrary subset of blocks on which those choices remain visible.
The composite covering radius is exactly

```math
q+(L-1)|\{j\in P:a_j=0\}|.
```

Thus one response query aggregates a chosen set of unit support effects before
the approximation error is paid.  Hamming packing yields an
`Omega_epsilon(w)` deterministic response-information lower bound at additive
error `epsilon*w` for every fixed `epsilon<1/8`; the latent block state gives a
matching `O(w)` upper bound on this restricted source family while answering
all unrestricted appended-fragment queries.  State and environment lengths
are only `Theta_epsilon(w)`.  The theorem and constants passed an independent
audit and small exhaustive verification.

This is primary progress: unlike the old one-bit environments, whose complete
response range is only `{1,2}`, the block query creates order-`w` separation.
It makes the framework more generative by predicting and then constructing a
joint exposure operation.  It does **not** determine the normalized metric
entropy of arbitrary syndrome supports or prove an exponential-in-`2^w`
rate.

The next single target is now the unrestricted-support dichotomy: construct a
subexponential-in-`2^w` net at error `epsilon*w`, or a superlinear/macroscopic
packing not forced by a supplied direct-sum decomposition.  The signing
problem remains outside this campaign.
