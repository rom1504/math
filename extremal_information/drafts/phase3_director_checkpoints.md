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

## Checkpoint 4: deterministic synchronization from cycle contraction

**Question.** Did the geodesic condition produce a strict response quotient,
or merely restate shortest-path optimality?

The cycle criterion itself is an exact reformulation, but its consequences
are genuinely compressive and have survived two independent audits.

1. For a complete quotient projection, one representative per affine fibre
   has joint additive-triangle defect at most three.  Vector-valued BLR
   stability supplies one linear section at uniformly bounded Hamming
   distance.  Replacing the whole raw support by its linear graph changes
   every rooted word metric under every appended raw support by at most
   eleven.  The state uses `O(w^2)` bits rather than a fibre truth table.

2. For an incomplete projection, quotient Cayley diameter `h` is the exact
   worst-case scale of hidden fibre information.  Stripping to one lift per
   quotient generator costs at most `2h` in every future context, and an
   explicit selector family attains `2h`.  A linear-section state gives the
   audited bound `10h+1`; hence `h=o(w)` is a genuine submacroscopic
   synchronization regime.

3. Synchronization does not make the response algebra exact.  Replacing
   finitely many raw sources telescopes, but union may create projected-zero
   cycles that occur in no source separately.  Those *mixed cycles* are the
   first unaccounted composition feature.

This reaches the requested generative threshold: a natural finite extremal
hypothesis forces a strict all-context quotient, with a sharp parameter
measuring when omitted rooted information becomes macroscopic.  It is not a
universal deterministic Parisi object; it is a theorem for binary extension
word landscapes.

### Immediate continuation

The next result must quantify mixed-cycle creation, not optimize the constants
eleven or `10h+1`.  We will test whether a closed affine/holonomy state
controls unions, identify the exact composition defect, and construct a
scalable equal-state counterexample if bounded algebraic data fail.  This
checkpoint therefore redirects the work and does not end it.

## Checkpoint 5: composition-created holonomy

**Question.** Is mixed-cycle language only another presentation of syndrome
decoding, or does it give a quantitative information-growth law?

The elimination formula is classical in substance, but the audited package
now goes further.

1. Kernel offsets on quotient-independent columns are removable by a linear
   shear.  Their complete invariant is the holonomy map on the quotient cycle
   space.  When labeled fragments are glued, the newly available invariant is

   ```math
   \operatorname{Hom}(Z/Z_{\rm loc},W),
   \qquad
   \kappa=\sum_j\dim U_j-\dim\sum_jU_j.
   ```

   Thus composition creates exactly `D kappa` binary degrees of labeled
   gluing freedom.  This is an exact gauge count, not automatically a scalar
   response lower bound.

2. The response content is nevertheless macroscopic.  One new mixed channel
   has an exponential family of rooted endpoint profiles separated by
   `Theta(D)`, so an `epsilon D` decoder with `epsilon<1/8` needs `Omega(D)`
   new gluing bits even though each fragment separately is shear-trivial.
   At every arity `r`, there are `r` fragments whose proper unions are all
   trivial but whose full union changes the antipodal distance by `D-r`.

3. Circuit defects give the robust amplification law

   ```math
   D-\ell(t)\le\nu\max_C d(C),
   ```

   and disjoint circuit blocks attain the nullity factor.  Composition turns
   microscopic defect into macroscopic response precisely when the product
   of per-circuit defect and new cycle rank is macroscopic.

4. On the positive side, exact linear graph sources have a closed affine-
   subspace feature state.  Its response error is at most `(r+1)/2`, where
   `r` is affine map rank rather than source count, and the order is sharp.
   This does not remove the separate synchronization cost for nonlinear raw
   transversals.  Complete unrooted pair-distance data also fail: tensor
   powers with identical pair laws retain a linear radius gap.

### Director judgment

This is a unified Level-3 result for binary extension landscapes.  The
positive quotient, exact gluing space, robust defect amplification, and
response packing are four faces of one relation-space mechanism.  It is not
a universal theory of arbitrary landscapes, and neither “holonomy” nor the
underlying exact sequence is new terminology.  The generative contribution
is the theorem that **relative gauge information is born on mixed relations,
with cycle rank controlling both algebraic growth and extremal error**.

### Immediate continuation

The one-channel lower bound does not yet show how much of the exact
`D kappa` gluing freedom remains operational at macroscopic distortion.  The
next theorem will study `kappa` independent mixed channels jointly.  It will
either construct `2^{Omega(D kappa)}` rooted-response profiles separated by
`Theta(D)`, or prove that unlabeled joint queries collapse most of those
degrees.  A code/Hausdorff packing gives a concrete falsifiable route.  Work
continues immediately; this checkpoint is not a stopping event.

## Checkpoint 6: full macroscopic rate of mixed holonomy

**Question.** Do `D kappa` labeled gluing degrees survive after channel labels
are removed and response is observed only to additive accuracy `epsilon D`?

Yes, in the Hamming extension model and throughout a nontrivial linear
regime.  The result has three independently audited layers.

1. For `k` parallel mixed channels with holonomy columns `V`, the unlabelled
   kernel-endpoint profile is the weighted Cayley norm

   ```math
   F_V(u)=\min_z\{2\operatorname{wt}(z)+\operatorname{wt}(u+Vz)\}.
   ```

   Exactly, it remembers the set of distinct generator columns of weight at
   least three; in particular a general change of channel basis is not an
   exact symmetry.  Macroscopically it is within `2 rank(V)` of distance to
   the image code.

2. Put the image codes inside one asymptotically good `[D,r,d]` host.  Its
   `k`-subspaces yield at least `q^(k(r-k))` profiles, and any two are
   separated in sup response distance by at least `d-2k`.  With
   `r=floor(D/4)`, `d>D/8`, and `k<=D/32`, there are at least
   `q^(3Dk/16)` profiles separated by more than `D/16`.  Therefore accuracy
   below `D/32` costs at least `(3/16)Dk log_2 q` deterministic bits, with the
   corresponding Fano mutual-information bound.

3. Both child fragments remain individually shear-trivial and every query is
   only an endpoint in the common kernel.  The information is genuinely born
   at composition.  The proof works over every finite field, so it is not a
   binary parity accident.

### Director judgment

This closes the promised dichotomy on the negative side and upgrades the
relation-space mechanism from gauge counting to operational response
complexity.  The genuinely reusable statement is not “dimension always
forces information.”  It is:

> Mixed relations create a carrier-valued compatibility state.  A
> presentation with submacroscopic penalty transfers Hausdorff packing of its
> zero sets into response packing; good carrier codes can then expose a
> constant fraction of the full compatibility information.

The qualifiers are essential.  Low-diameter carriers, redundant
presentations, synchronized sections, or a weak query family can collapse the
same algebraic gluing space.  Relation rank is the supply of potential
information; separated zero sets and cheap access are what make it visible.

### Immediate continuation

The next theorem will isolate this carrier-capacity principle abstractly and
try it outside finite-field Hamming geometry.  In parallel, an adversarial
construction will test whether any lower law depending only on carrier size
and relation rank is false.  This checkpoint redirects rather than ends the
campaign.

## Checkpoint 7: carrier capacity replaces parameter dimension

**Question.** Has the continuation produced a general composition law, or
only repackaged metric entropy?

The base identity is elementary: distance-to-set functions embed the
Hausdorff hyperspace isometrically into `ell_infinity`.  The project-level
result is nevertheless generative because composition supplies a weighted
carrier that was absent from either child.

1. If

   ```math
   F_\theta(x)=\min_{c\in C_\theta}
   \{d(x,c)+\pi_\theta(c)\},
   \qquad 0\le\pi_\theta\le p,
   ```

   then response distance and carrier Hausdorff distance differ by at most
   `p`, sharply.  Packing, covering, deterministic rate--distortion, and
   Fano bounds transfer with that one presentation-radius loss.

2. The same law is query-relative.  Under any query distribution `mu`, the
   `L^s(mu)` response metric is within `p` of the corresponding
   distance-transform metric.  A Hausdorff witness of size `Delta` exposes
   at least

   ```math
   (\Delta-2t-p)_+\,\mu(B_t)^{1/s}.
   ```

   Thus uniform hardness and diffuse-query hardness differ exactly by
   witness-neighborhood mass; the theorem does not silently turn a rare
   endpoint into average information.

3. The law produces full-rate packings in Hamming, Lee, flag-ultrametric,
   and rank-metric Cayley realizations.  The rank-metric construction uses an
   equilateral multiplication host, so it proves portability to a
   non-Hamming ambient model but not an intrinsically rank-geometric entropy
   theorem.

4. Raw relation dimension is decisively falsified as a universal proxy.
   Surjective maps can have identical image carriers; diameter-scale
   presentation cost can erase arbitrary carrier geometry; and a two-scale
   carrier of diameter `Theta(D)` can collapse `q^(Dk)` gauges to a quotient
   of `q^(rk)` states, or to only the subspaces of `F_q^r` for endpoint
   queries.

### Director judgment

The surviving object is not “holonomy dimension.”  It is the **presented
carrier response**: a composition-created subset in the query metric,
together with the cost of accessing its points.  Its macroscopic complexity
is the exposed metric entropy of that carrier class after subtracting the
presentation scale.  This unifies the positive and negative finite-field
examples and gives a second nontrivial model through rank-metric shortcuts.
It remains a theorem for distance-transform/min-plus landscapes, not a
universal theory of arbitrary extremal responses.

### Immediate continuation

The two-scale counterexample suggests the complementary positive compression
theorem.  We next test whether a coarse metric quotient with small fibres and
controlled lifting makes every presented carrier response a function of the
quotient carrier, with nonaccumulating error under future min-plus contexts.
This is a natural deterministic synchronization hypothesis rather than an
assumed linkage condition.  Work continues immediately.

## Checkpoint 8: deterministic metric synchronization

**Question.** Does the two-scale collapse expose a natural compression
theorem, or only an engineered example?

It exposes a general, independently audited sufficient condition.  If an
onto one-Lipschitz map of query metrics has fibre diameter `a`, lift defect
`b`, and a carrier has presentation radius `p`, then its complete endpoint
profile is uniformly within `a+b+p` of distance to the projected carrier.
The three terms are jointly sharp under the stated axioms.  Hausdorff nets of
projected carriers therefore transfer directly to response nets.

The error survives every fixed min-plus continuation without amplification.
For additive carriers and a homomorphic quotient, the projected carrier
updates exactly; translation-invariant metrics give the corresponding exact
response infimal convolution.  The maintained state includes one scalar
presentation-radius certificate.

The audit identified the precise strictness conditions.  This is compression
only if projected carrier entropy is actually smaller, the error is subscale,
the declared composition descends, and the radius certificate remains
controlled.  Rank-row projection and puncturing prove small-error factors but
do not automatically prove an entropy reduction.  The two-scale model meets
all four conditions.

### Immediate continuation

The lower carrier theorem and upper quotient theorem currently use unrelated
geometric certificates.  The next theorem will define their scale-dependent
linear ranks and test whether fibre geometry forces an inequality between
them.  Hamming puncturing suggests a generalized Singleton argument.  Work
continues immediately.

## Checkpoint 9: scale-rank duality and an intrinsic rank model

**Question.** Can the lower and upper laws be related without enumerating the
full response class?

Yes.  Let `s_W(Delta)` be the largest dimension of a linear carrier host with
minimum nonzero metric weight greater than `Delta`.

1. Its `k`-subspaces give `q^(k(s-k))` response profiles separated by more
   than `Delta-2k`.
2. A dimension-`r` metric synchronization quotient decodes every profile
   from a projected subspace, using at most
   `sum_(j<=k){r bracket j}_q` states and error `a+b+2k`.
3. Fibre diameter `a` forces

   ```math
   s_W(a)\le r.
   ```

   A separated host must inject into the quotient.  This is exactly the
   puncturing proof of Singleton in Hamming space and its row-projection
   analogue in rank metric.

The two-scale model attains the inequality at every scale: its separated-rank
curve is `D`, then `r`, then zero.  More substantially, a self-contained
Gabidulin host of dimension `rD` and minimum rank `D-r+1` turns the
rank-metric Singleton geometry into response information.  At
`r=floor(D/2)` and `k<=D/16`, it produces at least `q^(kD^2/3)` profiles
separated by more than `3D/8`, forcing `(1/3)kD^2 log_2 q` bits below error
`3D/16`.  This removes the earlier equilateral-host limitation.

### Director judgment

The program has moved from local examples to one genuine two-sided law for
linear presented carriers.  Algebraic relation rank creates a potential
state; metric carrier separation turns it into information; metric quotient
fibres bound how much can be compressed; and Singleton duality makes those
certificates compatible.  The same statement explains Hamming, two-scale,
Lee/flag, and intrinsic rank-metric behavior.

It is still not a universal extremal-information theory.  The carrier
representation is special to distance-transform/min-plus responses, and the
best separated rank need not equal the best synchronizing quotient rank.

### Next theorem

Determine whether those two scale ranks are asymptotically dual in a natural
class, or exhibit a leading gap and identify the invariant between them.  A
mere restatement as covering dimension will not suffice; the result must
predict a response exponent or a strict composable quotient.

## Checkpoint 10: anticode exactness and a real duality gap

**Question.** Is synchronizing quotient rank another optimization with no
closed characterization?

No.  For an `N`-dimensional translation-invariant linear carrier, let
`A_W(a)` be the largest dimension of a linear anticode of diameter at most
`a`.  The least dimension of an exact `(a,0)` synchronization quotient is

```math
N-A_W(a).
```

The lower bound is kernel dimension.  The upper bound uses the canonical
coset metric on `W/K`, which is one-Lipschitz and lifts every quotient
distance exactly.  Consequently

```math
s_W(a)+A_W(a)\le N.
```

This is the precise code--anticode boundary between carrier packing and
carrier compression.  It is tight for the two-scale metric and for rank
metric, where Gabidulin codes and row-supported anticodes give
`s_W(a)=D(D-a)` and `A_W(a)=Da`.

It is not universally tight.  In binary Hamming space,
`A_W(floor(delta D))=floor(delta D)`, while sphere packing gives

```math
\liminf {N-A_W-s_W\over D}
\ge H_2(\delta/2)-\delta>0.
```

Thus the session ends with a decisive falsification of universal scale-rank
duality, not an unresolved slogan.

### Director judgment and continuation target

The strongest surviving abstraction is now a three-layer law:

1. mixed relations create a presented carrier;
2. carrier packing and quotient anticodes bound its response information;
3. the code--anticode gap measures what neither certificate resolves.

The next theorem is the Hamming Grassmannian problem: determine whether
families of subspace carriers not contained in one common separated host can
fill the leading code--anticode gap.  This checkpoint is followed by full
independent audit, surface promotion, and reproducibility checks; it is not
being used to return to the original signing problem.

## Checkpoint 11: the Hamming gap contains coding theory and rooted lifts

**Question.** Is the missing Hamming Grassmannian exponent another quantity
that should have a simple formula between separated rank and anticode
codimension?

Not in unrestricted form.

1. For binary lines, the packing number is between
   `A_2(D,t+1)-1` and `A_2(D,t+1)`.  Its exponent is the classical nonlinear
   coding exponent.  An actual probabilistic line cover has exponent at most
   `1-H_2(delta)`, proving that puncturing is exponentially nonminimal without
   confusing packing and covering.
2. Directed Grassmannian balls have an exact formula in terms of linear flats
   inside the quotient coset-leader ball.  This sparse-flat spectrum is useful
   but not two-sided sufficient: direct-sum carriers can have isometric
   quotient norms and a linear rooted Hausdorff gap.
3. A general injection-distance construction uses the low words in `C+C'`,
   but its asymptotic exponent is always bounded by the common-host Gilbert
   exponent.  A four-carrier finite witness and a scalable seven-plane
   alphabet show that common hosts are nevertheless not exactly complete.
4. Systematic charts reduce to `2^k`-ary column codes.  At `k=Theta(D)` that
   alphabet admits MDS codes meeting Singleton, so ordinary column coding
   cannot improve the quotient exponent; same-input recoupling is the missing
   constraint.

### Director judgment

The sparse-flat spectrum is retained as an exact one-sided counting tool and
rejected as a response state.  The actual middle geometry contains nonlinear
coding, kernel/lift rooting, and orientation.  A universal closed Hamming
exponent is not a sensible immediate objective.

### Immediate continuation

The seven-plane packing exposed a useful compositional ingredient—both
directed gaps add exactly—but its presentation toll is too large.  Work
continues with the minimal question: when does a finite two-sided carrier
alphabet retain a positive response margin after presentation?

## Checkpoint 12: exact directed-response algebra

**Question.** Can multiple channels be evaluated jointly so their useful
orientation adds before the final absolute value, without paying matching
channels separately?

Yes.  For local responses `f_a`, the finite directed table

```math
r(a,b)=\sup_x(f_a(x)-f_b(x))
```

determines every pairwise uniform distance after arbitrary direct-product
composition:

```math
\|F_{\boldsymbol a}-F_{\boldsymbol b}\|_\infty
=\max\left\{
 \sum_i r(a_i,b_i),
 \sum_i r(b_i,a_i)
 \right\}.
```

This is exact.  If local carriers have both directed gaps at least `d` and
presentation radius `p`, each differing channel contributes at least `d-p`
in either chosen orientation.  An outer code of relative distance `rho`
therefore gives gap `(d-p)rho m`; matching channels pay zero.  This improves
the coarse global-carrier estimate by `p(1-rho)m`.

Two independent model tests are proved and exhaustively checked.

* Seven binary simplex lines have `d=4,p=2`; a relative-`3/4` outer code gives
  response gap `3m/2`.
* Seven `F_8` multiplication lines in rank metric have `d=3,p=2`; the same
  outer code gives gap `3m/4`.

Both families have positive response rate
`(0.0573549...-o(1))m` bits.  The Hamming family is not contained in any
common growing-distance host.

### Director judgment

This is the strongest result of the checkpoint.  The `q x q` directed table
is a strict query-relative quotient of the full local functions, has an exact
composition algebra, and produces a theorem in two metric models.  The
framework is genuinely generative at Level 3 for product-composed presented
responses, while remaining Level 2 globally.

It is related to concatenated coding but not merely a relabeling: the new
content is the exact signed response algebra and the proof that presentation
cost is charged only on differing channels.  It directly explains when
composition converts microscopic distinctions into macroscopic extremal
information.

### Next theorem

Determine when a nontrivial min-plus continuation preserves the directed
table up to sublinear cumulative loss.  Generic continuation is only
nonexpansive.  A successful theorem must derive near-equality from exposed
minimizers or a strict synchronization quotient without storing the full
kernel rows.  The original signing problem remains out of scope.
