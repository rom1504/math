# Phase 2 checkpoint 1: adversarial verification

## Verdict convention

- **ACCEPT**: correct as stated, including its declared scope.
- **CORRECT AFTER EDIT**: the mathematical core is correct, but a stated
  interpretation, scope, or complexity claim must be narrowed.
- **REJECT**: false in the declared scope.

No theorem below is rejected.  The transfer-kernel results are correct but
mostly a sharp response-theoretic packaging of max-plus dynamic programming.
The outer-spectrum theorem is correct but is only minimal for the deliberately
declared all-temperature pressure experiment.  It is not the minimal rooted
augmentation for covering radius or Cartesian-product composition.

## Summary table

| Result | Verdict | Required action |
|---|---|---|
| FG.1(a), exact endpoint quotient | **ACCEPT** | None. |
| FG.1(b), isometry and response minimality | **ACCEPT** | Keep the arbitrary, unbounded endpoint-field hypothesis prominent.  Say “coarsest response-equivalence quotient,” not “smallest encoding.” |
| FG.1(c), optimal-transport roof | **ACCEPT** | Keep the nonempty-fibre hypothesis; otherwise restrict transport plans to feasible endpoint pairs and allow `-infinity`. |
| FG.2, max-plus closure and error propagation | **ACCEPT** | None. |
| FG.3, bounded-separator closure and bit count | **CORRECT AFTER EDIT** | Replace “number of full state energies” by “number of entries in the full landscape table”; state additive factor allocation and treatment of infeasible boundary assignments explicitly. |
| FG.4, universal separator packing | **ACCEPT** | Its sharpness is only for the universal arbitrary-kernel class.  Do not transfer it to a fixed local constraint language without a realizability packing. |
| OS.1(1), pressure quotient | **ACCEPT** | Explicitly call it the quotient for the *complete pressure curve*. |
| OS.1(2), zero-temperature limit | **ACCEPT** | None. |
| OS.1(3), product algebra | **ACCEPT** | None. |
| OS Section 3, strict nonreconstruction example | **ACCEPT** | Exact census independently verified. |
| “minimal rooted augmentation/intermediate state” interpretation | **REJECT** | Replace by “strict root-averaged quotient for the pressure experiment.”  Radius alone is sufficient and minimal for radius plus Cartesian products. |

## 1. Transfer-kernel audit

### FG.1(a): exact quotient — ACCEPT

Conditioning on `(ell,r)=(a,b)` gives (FG.4) immediately.  The conclusion is
exactly scoped to unary endpoint rewards.  It does not cover a future query
touching an internal state, and the draft correctly says so.

### FG.1(b): max-plus isometry and minimality — ACCEPT

The reverse inequality in (FG.5) is valid despite the fact that the query is
separable as `f(a)+g(b)`.  For a selected pair `(a0,b0)`, set both selected
coordinates to zero and all other coordinates to `-T`.  If

```math
T>max\{\operatorname{osc}K,\operatorname{osc}K'\},
```

then `(a0,b0)` is the unique maximizing endpoint pair for both kernels.  This
recovers the difference of that entry.  Choosing an entry of maximum
absolute difference proves the reverse inequality.  The usual maximum
inequality proves the forward inequality.

Thus equality of complete endpoint-response functions is equivalent to
equality of kernels.  This proves operational minimality as an equivalence
quotient.  It does **not** prove a canonical minimum number of real
coordinates or bits for an unrestricted real encoding; the finite packing in
FG.4 is what supplies a bit lower bound after an accuracy and kernel class are
declared.

The use of arbitrary fields is essential.  For example, under a bounded
field experiment, take a kernel whose `(1,1)` entry is so far below the other
three entries that it is never exposed.  Changing only that entry changes
the sup norm but no allowed response.  The draft already records this scope
boundary; it should remain adjacent to every use of the isometry.

### FG.1(c): response roof — ACCEPT

A probability distribution on states induces a coupling `pi` of the two
endpoint marginals.  Its conditional mean energy in fibre `(a,b)` is at most
`K(a,b)`.  Conversely, mixing one fibre maximizer for each `(a,b)` realizes
every coupling when all fibres are nonempty.  Hence (FG.6) is exact and its
values at pairs of simplex vertices recover all kernel entries.

There is one necessary boundary convention.  If a later version permits
empty fibres, (FG.6) must maximize only over couplings supported on feasible
pairs, equivalently use `K(a,b)=-infinity` on empty fibres.  Without that
change the displayed formula is false.  The present theorem assumes every
fibre nonempty, so it is correct as written.

### FG.2: max-plus composition — ACCEPT

For fixed exposed endpoints `(a,c)`, maximizing over the shared endpoint `b`
gives exactly `K odot L`.  Associativity and (FG.9) follow by expansion and
the elementary Lipschitz property of a maximum.  Iteration gives (FG.10).
No hidden independence or uniqueness assumption is used.

This is standard transfer-matrix/max-product algebra.  The response isometry
adds a clean operational metric interpretation, but neither the composition
law nor its additive error estimate is a new mechanism beyond max-plus
dynamic programming.

### FG.3: separator corollary — CORRECT AFTER EDIT

The one-sided profile is the coarsest exact quotient for arbitrary separator
potentials by the same pinning argument.  A two-sided width-`w` strip has
`q^(2w)` entries, and if each entry is an integer in `[-NW,NW]`, its raw exact
storage is

```math
O(q^{2w}\log(1+NW))
```

bits, plus a constant-size feasibility flag per entry if boundary fibres may
be empty.  The analogous one-sided count is correct.

Two exact edits are required:

1. Replace “The number of full state energies can be exponential in `N`” by
   “The full landscape table can have exponentially many state entries.”
   Under the stated integer and magnitude hypotheses there are at most
   `2NW+1` distinct numerical energy values, so the existing phrase can be
   read as false.
2. In the tree/path-decomposition sentence, state that every local factor is
   allocated exactly once and that the split is additive conditional on the
   separator.  Otherwise factor double counting can invalidate the displayed
   gluing rule.

Calling this “deterministic synchronization” is permissible only if it is
defined here to mean exact conditional screening by a separator.  It is the
Markov/separation property underlying variable elimination, not an analogue
of Ghirlanda--Guerra/Panchenko synchronization.  The latter identification
would be an overclaim.

### FG.4: separator information lower bound — ACCEPT, universal scope only

The class realizes every binary `Q by Q` kernel with one state per endpoint
pair.  By FG.1(b), distinct binary matrices are response-distance one.  If a
single deterministic message decoded both response functions to error
`epsilon<1/2`, their distance would be at most `2 epsilon<1`, a contradiction.
Thus `Q^2` bits are necessary and sufficient on this binary class.

The continuous cube has sup-metric packing and covering logarithms

```math
Theta(Q^2 log(B0/epsilon))
```

in the stated nontrivial regime.  To make this completely self-contained,
say explicitly that the one-state-per-fibre construction also realizes every
real kernel in the cube, not just binary kernels.

The persistence calculation is exact: for `C>1`, every off-diagonal term in
`(A odot I_C)(a,c)` is at most `1-C<0`, whereas the diagonal choice equals
`A(a,c)>=0`.  Hence `A odot I_C=A` at every chain length.

The scope caveat is substantive.  This construction uses an arbitrary factor
on the entire left/right boundary pair.  It does not establish a
`q^(2w)`-bit lower bound for kernels generated by a fixed bounded-arity CSP
language, a fixed trellis, planar factors, linear-code constraints, or any
other restricted semigroup.  Such a claim needs an equally large separated
packing *inside the realizable kernel class*.  The draft states this at the
end; the executive verdict should preserve the word “universal.”

## 2. Is feature growth more than dynamic programming?

The mathematical statements are sound, but the director checkpoint currently
leans too strongly toward novelty.

- Conditional maxima on a separator, max-plus multiplication, associativity,
  and reusable endpoint tables are precisely standard variable elimination,
  Viterbi/max-product, and tropical transfer matrices.
- The roof formula is standard convex duality plus an optimal-transport
  coupling polytope.
- The isometry is an elementary consequence of allowing unbounded unary
  endpoint fields.
- The binary lower bound is the full-cube packing forced by the choice of a
  model class that realizes an arbitrary table.

What is useful to this program is the **unified operational statement**:
separator tables are exact query quotients, their response distortion is sup
distance, and their universal metric entropy is exponential in separator
width.  This is a clean Level-2 explanatory theorem and a valuable baseline
for asking about restricted realizable semigroups.  It is not yet a new
Level-3 mechanism.  The proposed next problem—showing that a structured
realizable kernel semigroup has unexpectedly low response-metric entropy—is
the point at which genuinely new theory could begin.

## 3. Outer-spectrum audit

### OS.1(1): quotient for the pressure curve — ACCEPT

For integer distances in `{0,...,n}`, `O_A` is a finite polynomial with
nonnegative integer coefficients.  Knowledge of it gives `P_A(beta)`.
Conversely, equality of `P_A` on any nonempty real interval gives equality of
`O_A(e^beta)` there, hence equality of the two polynomials.  Therefore `O_A`
is exactly the response-equivalence quotient for the complete pressure-curve
experiment.

This is not a result about a single bounded temperature, a finite temperature
list independent of `n`, or covering radius alone.  The phrase “complete
pressure experiment” must remain in every minimality claim.

### OS.1(2): radius and zero temperature — ACCEPT

The largest occupied degree is the covering radius.  The largest-term
log-sum-exp sandwich proves (OS.3)--(OS.4), with no missing normalization.

### OS.1(3): product algebra — ACCEPT

For the `l_1` product metric,

```math
delta_(A times B)(x,y)=delta_A(x)+delta_B(y),
```

so the outer polynomials multiply, pressures add, and degrees/radii add.
All three statements are exact.

### Exact `Q_4` example — ACCEPT

Independent enumeration gives

```text
A outer: 0:4, 1:8, 2:4
A inner: 0:4, 1:8, 2:4
B outer: 0:4, 1:8, 2:4
B inner: 0:4, 1:6, 2:4, 3:2
```

Thus (OS.7)--(OS.8) are exact.  Different inner enumerators prove the sets are
not Hamming-isometric.  Polynomial multiplication preserves the common outer
spectrum under Cartesian powers, while the distinct inner polynomials remain
distinct after taking positive powers, so the nonreconstruction claim is
valid at every power.

The claimed scope failure under named-coordinate puncturing is also visible
in this same example.  Puncturing coordinate 2 (numbered from the left
starting at 1) gives outer counts `4+4z` for `A` but `3+4z+z^2` for `B`.

### Bit complexity — ACCEPT with wording discipline

In `Q_n`, the `n+1` coefficients are integers at most `2^n`, so the raw
coefficient vector uses `O(n^2)` bits.  The full raw labeled map has `2^n`
entries.  State this as raw/worst-case representation complexity, not as the
description length of every individual succinct code.

## 4. The “minimal rooted augmentation” claim is overstated

The outer spectrum is not minimal for the zero-temperature target under the
declared Cartesian product.  The scalar covering radius already satisfies

```math
rho(A times B)=rho(A)+rho(B),
```

and answers the radius query exactly.  It is therefore the coarsest exact
quotient for that much smaller experiment.  For an explicit strict
separation, in `Q_2`

```math
A={00,01},       O_A(z)=2+2z,
B={00,01,10},    O_B(z)=3+z,
```

yet both radii equal one and all Cartesian-power radii are respectively the
same additive scalar.  Thus the full outer spectrum contains unnecessary
information for covering radius plus product composition.

Nor is `O_A` literally a rooted state: it averages away root labels.  It is a
histogram of rooted distances, or a **root-averaged** state.  It also does not
sit between inner replica data and the labeled distance map in an information
partial order.  The present `Q_4` example has the same outer spectrum and
different inner data, while the earlier `C,D` example has the same inner data
and different radii, hence different outer spectra.  Inner and outer spectra
are incomparable quotients, although both are quotients of sufficiently rich
labeled information.

Required replacement for the checkpoint claim:

> The outer spectrum is a strict, nonreconstructive root-averaged quotient
> that is exactly sufficient and minimal for the complete uniform-root
> pressure curve, and it closes under Cartesian products.

Do not call it “the minimal rooted augmentation that restores compositional
sufficiency” without specifying that exact pressure query.  For the harder
rooted/coupled experiments posed by the program it has not restored
sufficiency.

## 5. Is the outer-spectrum result more than log-sum-exp?

No, not mathematically.  It is the moment-generating polynomial of the outer
distance distribution; product multiplication is convolution under an
additive product metric, and the radius limit is the elementary largest-term
log-sum-exp estimate.  Its value is taxonomic: it supplies a clean known-model
example of a strict query quotient and shows exactly how changing the query
family changes the minimal state.  It should not count as a substantive new
theorem or as evidence that richer code couplings have been compressed.

## 6. Director recommendation

Retain FG.1--FG.4 after the two wording fixes because they give a rigorous
baseline and one genuine external model with sub-landscape closure.  Count
them as one explanatory checkpoint, not as several independent discoveries.
Retain OS.1 as a compact known-model validation example, but downgrade the
minimal-rooted interpretation as above and do not count it as a substantive
result.

The next theorem must leave these taut regimes.  A qualifying advance would
show either:

1. a restricted local/trellis/code kernel semigroup has response-metric
   entropy asymptotically below its ambient `q^(2w)` table size while remaining
   closed under composition; or
2. a strictly coarser rooted statistic than the labeled distance map closes
   under a composition richer than Cartesian product and answers a query not
   already reducible to a single additive scalar.

Until one of these is proved, checkpoint 1 is a correct consolidation of
max-plus dynamic programming and distance-enumerator algebra, not yet a new
generative theory.
