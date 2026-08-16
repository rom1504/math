# Adversarial audit: metric-quotient synchronization

**Object audited.**
[`phase3_metric_quotient_synchronization.md`](phase3_metric_quotient_synchronization.md).

**Verdict.** MQ.1, its constants, the Hausdorff-net consequence, and the
three concrete metric synchronizations are correct.  Min-plus
nonamplification is also correct after making the extended-real convention
explicit.  MQ.4 is an exact algebra of carrier **sets**, but two scope repairs
are needed: translation invariance is required if it is also meant as an
exact min-plus algebra of the displayed distance responses, and growth of the
presentation-radius bound is sufficient rather than necessary for loss of
accuracy.  The theorem is a genuine strict quotient only after a separate
entropy comparison for the projected carrier family.

## 1. MQ.1 and the constant `a+b+p`

The proof is valid under the stated finite hypotheses.  More explicitly, for
`y in varpi C` nearest to `varpi x`, choose `c in C` over `y`.  The lifting
condition produces `z` over the same `y` with

```math
d_X(x,z)\le d_Y(\varpi x,y)+b,
```

and the fibre bound gives `d_X(z,c)<=a`.  Thus

```math
d_X(x,C)\le d_Y(\varpi x,\varpi C)+a+b.
```

The presentation contributes at most `p`.  In the other direction,
one-Lipschitzness alone gives

```math
d_Y(\varpi x,\varpi C)\le d_X(x,C)\le F_{C,\alpha}(x).
```

There is no hidden need for a section, a group structure, or compatibility of
`C` with fibres.  The product example attains `a+p` when `b=0`, so those two
terms cannot be removed uniformly.  The lifting clause is genuinely needed:
small fibres and a one-Lipschitz quotient do not bound how much shorter an
arbitrarily declared quotient metric can be.

Minor hypothesis clarification: declare `a,b,p>=0`.  Finiteness guarantees
all nearest points; for infinite spaces the statements would need infima and
an attainment or approximation argument.

## 2. Quotient nets

MQ.2 is correct.  If quotient carriers `A,A'` are within Hausdorff distance
`eta`, then

```math
\|d_Y(\mathord\cdot,A)-d_Y(\mathord\cdot,A')\|_\infty=\eta.
```

Pullback along `varpi` cannot increase this error, and MQ.1 adds
`a+b+p`.  Surjectivity even makes the pullback norm equal to the norm on
`Y`, although only the upper bound is used.

The net statement counts projected **carriers**, not merely points or vector
dimensions.  It is a strict response quotient only if

```math
\operatorname{Cov}(\{\varpi C\},d_H^Y,\eta)
```

is demonstrably smaller than the corresponding original response/carrier
entropy.  The identity map satisfies MQ.1 with `a=b=0` and gives no
compression.  Likewise, a low-dimensional target does not by itself help for
arbitrary subsets: its hyperspace can still have `2^{|Y|}` elements.

## 3. Min-plus continuation

MQ.3 is the standard order argument and is correct for any fixed kernel.
The clean hypothesis is either

* real-valued `f,g` with `||f-g||_infty<=epsilon` and
  `K:X times X->R union {+infinity}` chosen so every output is finite;
  or, more safely,
* extended-real functions satisfying the order inequalities
  `g-epsilon<=f<=g+epsilon`, with no indeterminate `(+infinity)-(+infinity)`
  norm notation.

Using `inf` rather than `min` makes the statement valid without attainment.
For the present finite applications, real-valued profiles and kernels in
`R union {+infinity}` suffice.

The nonamplification conclusion applies when the future kernel is fixed
independently of which of `f,g` is supplied.  It also applies to a sequence
of kernels, including kernels between changing finite domains.  It does not
say that a quotient carrier can be updated from itself under an arbitrary
future operation; that is a separate closure requirement, correctly noted in
the draft.

## 4. Additive carrier algebra

The set identity

```math
\varpi(C+D)=\varpi C+\varpi D
```

is exact for a homomorphism, and infimal convolution of the presentations is
associative with the certified bound

```math
0\le\alpha\square\beta\le p_C+p_D.
```

Two scope points matter.

First, if MQ.4 is intended only as an algebra of carrier sets, its hypotheses
are sufficient as written.  If it is intended to identify composition of the
responses with a min-plus context, require `d_X` to be translation invariant.
Then one has the exact formula

```math
F_{C+D,\alpha\square\beta}(x)
=\min_{d\in D}\{F_{C,\alpha}(x-d)+\beta(d)\}.                \tag{A.1}
```

Without translation invariance, (A.1) need not hold even though the carrier
set identity remains true.  Translation invariance of `d_Y` is natural if
the projected distance profiles are also to carry the same additive
interpretation, though it is not needed for the bare set identity.

Second, `p_C+p_D` is an upper certificate, not necessarily the exact parent
presentation radius.  Infimal convolution can reduce it drastically, and a
large certified radius does not imply a large actual decoding error.  Replace
the phrase “precisely while the total presentation radius ... is
submacroscopic” by:

> the theorem guarantees submacroscopic error whenever the actual or tracked
> presentation-radius bound, together with `a+b`, is submacroscopic.

For variable-depth composition, the maintained feature is really
`(varpi C, p_certificate)` rather than `varpi C` alone.  This adds only a
scalar bookkeeping state but should be stated.

## 5. Model checks

### Two-scale finite-field metric

MQ.5 is correct for `s>=0`.  The displayed formula is a metric, `varpi` is
one-Lipschitz, and each fibre has diameter one.  Given `x` and a target `y`,
choose a linear splitting of the surjection and change only the quotient
component; the distance is exactly `(s+1)1_{varpi x ne y}`.  Thus `b=0`.

The count of projected subspaces is independent of `D` only with `q,r`
fixed.  State this explicitly.  The sharper `q^{rk}` labeled-map claim is not
a consequence of MQ.1; it relies on the separate aggregate-fibre argument in
the cited carrier-relation note, as the draft already indicates.

### Rank-metric row projection

MQ.6 is correct.  Deleting rows cannot increase rank.  A fibre difference is
supported on `D-r` rows and has rank at most `D-r`.  To lift a target top
block, replace the top rows and leave the bottom rows unchanged; the rank of
the full difference equals the rank of the top-block difference, so `b=0`.

When `D-r=o(D)` and `p=o(D)`, the response error is submacroscopic.  This is
not automatically a strict information quotient: projected carriers can
still have essentially the full entropy of the chosen carrier class.  The
example establishes a small-error non-Hamming factor; strict compression
requires the MQ.2 covering number to be smaller.

### Hamming puncturing

MQ.7 is correct.  Puncturing is one-Lipschitz, fibres have diameter `h`, and
a target prefix is lifted by retaining the omitted coordinates, so `b=0`.
The `h+p` error is stable under subsequent fixed min-plus contexts.  Updating
only the punctured carrier under composition additionally requires an
operation for which puncturing is a homomorphism, such as Minkowski addition;
arbitrary future operations do not follow from MQ.3.

## 6. Is this more than a rephrasing?

MQ.1 is an elementary approximate-submetry factorization of distance-to-set
responses.  It should not be advertised as a new general synchronization
theorem in metric geometry.  Within the extremal-information program it is
non-tautological in a useful but limited sense: one quotient works uniformly
for every carrier, its hypotheses are checked on the ambient query metric
rather than on the full response landscape, and MQ.4 can provide an exact
feature update.

The strict-content test is therefore:

1. `a+b+p` is below the target response scale;
2. the projected carrier family has a proved smaller Hausdorff covering rate;
3. the declared composition descends exactly to projected carriers; and
4. the presentation-radius certificate remains controlled.

Without all four, the construction is a correct approximation identity but
not yet a strict composable quotient.

## Required scope fixes

1. Add `a,b,p>=0` and a clean finite/extended-real convention.
2. In MQ.4, distinguish the always-valid carrier-set algebra from the
   response min-plus algebra, which needs translation-invariant metrics.
3. Replace “precisely while” by the one-way guarantee described above and
   include the scalar presentation-radius certificate in the maintained
   state.
4. In MQ.5 say “for fixed `q,r`.”
5. Qualify MQ.6 and MQ.7 as small-error quotients; they become strict
   information quotients only after a projected-carrier entropy bound.
