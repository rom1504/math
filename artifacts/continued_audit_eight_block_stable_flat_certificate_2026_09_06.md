# Independent audit: four- and eight-point stable flat-profile codes

Date: 2026-09-06. Full-read independent reconstruction of
`continued_director_stable_quad_profile_count_2026_09_06.md` and
`continued_feedback_eight_block_stable_flat_certificate_2026_09_06.md`.
Both class-specific theorems pass. Neither is a full signing cap bound.

## 1. The available local information is exactly sufficient

For a Walsh supported spin xi of weight k=pm, epsilon-flatness gives

```math
\|\xi\|_{U^2}^4\le8p^2(m^{-1}+\epsilon^4),
\qquad |\mathbb E\xi|\le\sqrt p(m^{-1/2}+\epsilon).
```

These follow directly by writing the Fourier coefficients as a Boolean
sign times sqrt(k), plus the magnitude error. In the fourth-power sum,
the flat contribution is p^2/m and the error contribution is at most
p^2 epsilon^4. The mean bound uses the zero-frequency coordinate.

A high-probability selector event ensures `||1_T-p||_{U^2}=o(1)`.
The zero/positive/negative indicators are linear combinations of this
centered selector and xi, so they have the required small centered U^2
norms and near-product symbol frequencies.

Any at most four distinct vertices of an affine three-cube are either
affinely independent, or (only at size four) a parallelogram. Uniform
injective affine images of the first kind differ from independent
sampling by O(1/m). For a parallelogram, the exact Fourier expansion
forces a common frequency and Holder bounds by the product of the four
U^2 norms. Thus every marginal of size at most four is uniformly within
`O(epsilon)+o(1)` of the ternary product law. No eight-wise independence
and no U^3 estimate enter the proof.

## 2. The four-point code and its limited numerical gain

The three retained characters of a four-point plane project orthogonally
away from its fourth character. Translation within a plane changes the
omitted character only by a sign, so the product-prior conditional cost
is independent of the coset representative. Averaging the conditional
cost over directions selects one pair of directions at o(m) description
cost. Given the exact projections, product conditional probabilities
sum to one; the number of preimages with code cost at most R is at most
exp(R), without assuming independent actual blocks.

The retained 3m/4 Fourier signs give projected-coordinate error at most
epsilon^2 k. True projected coordinates are quarter-integers in
[-3/2,3/2], a 13-symbol alphabet. A rounding error costs at least 1/64,
so at most 64 epsilon^2 k corrections are required. These reconstruct
the claimed stable count and the denominator binomial(m,k).

The director's rational script was rerun to the dedicated scratch output
`tmp/continued_audit_quad_certificate_replay_2026_09_06.json`. It verifies
at p=15/16 and t=33/32 the exact-class tilted exponent
`-0.000033507010118643... < -3/100000`. The source correctly notes that
the typical Gaussian profile has positive pressure at this tilt. Thus
this negative class exponent cannot be promoted to a full bound.

## 3. Independent exact reconstruction of the eight-point dual

The original optimizer-free checker was read completely and replayed.
Additionally, the separate audit checker

`computations/continued_audit_eight_block_dual_replay_2026_09_06.py`

imports only the proposed 71 integer dual coefficients and rebuilds the
entire verification. It does not use the original sparse matrices,
orbit-union algorithm, floating marginal table, optimizer, or log-bound
function. Its own arithmetic uses integer tuples and rational fractions.

The audit checker enumerates all 6561 ternary patterns and their exact
product masses, obtaining 1537 outputs at frequencies {0,1,2,4}. Each
permutation of this affine frequency basis determines one affine map;
together with global sign reversal these produce 48 signed coordinate
actions and 255 pattern orbits. Product mass is preserved. The output
mass is checked exactly constant on each orbit.

All 1697 degree-at-most-four indicator monomials are rebuilt in the
specified ordering. On every orbit, the averaged integer polynomial
majorizes an independently computed rational upper enclosure of
`log P_Y(Y)`. The minimum certified slack is greater than
`6.0027e-11`; no constant correction or coefficient adjustment is used.
Since the product law is invariant under every signed coordinate action,
the exact expectation of the averaged polynomial is its original
product-law expectation. The result is precisely

```math
A=-65899554490727779/10485760000000000.
```

Symmetrizing an arbitrary feasible four-wise law preserves both its
constraints and the objective, so this orbit certificate applies to all
feasible laws, not just symmetric ones. For approximately correct
marginals, the same fixed finite polynomial gives a uniform error
`O(epsilon)+o(1)`. There is no limiting optimization or bounded-dual
assumption to justify: the coefficients are explicit fixed numbers.

The original checker stores rhs entries as floating values, but these
particular entries `(15/32)^r`, r<=4, are exactly binary-representable.
Its conversion to Fraction is therefore exact. The independent checker
instead constructs them as fractions from the outset.

## 4. Encoding normalization and strict finite-width margin

The chosen four characters of each eight-point block correspond to
exactly m/2 global Fourier coefficients. After their signs are stored
and their magnitudes replaced by sqrt(k), quotient Walsh Parseval gives
total squared error in the block sums at most `8 epsilon^2 k`.
Block sums are integers in [-8,8]. A wrong rounding costs at least 1/4,
so at most `32 epsilon^2 k` of the m/2 sums require correction. The
relative error fraction is therefore exactly `64p epsilon^2`, and
the correction entropy is one HALF of
`h(64p epsilon^2)+64p epsilon^2 log17` per input coordinate.

For an actual supported spin the negative log product-prior probability
is exactly `m[h(p)+p log2]`, irrespective of sign imbalance. The
conditional recovery cost adds the sum of `log P_Y(Y)` over the m/8
blocks. Averaging over affine partitions and using the fixed dual picks
one partition with the claimed total cost. Encoding the directions and
a complementary basis costs only O((log m)^2), not a linear entropy.

This gives the zero-width count rate

```math
c_8=(p+1/2)\log2+A/8=0.2108151647996559\ldots.
```

The already audited near-flat permanent perturbation contributes
`(1/2)log[(1+exp(-4t))/2]+4t epsilon+o(1)` uniformly over every removed
coordinate. At p=15/16 and t=4, adding the weave tilt gives

```math
p\log2+A/8+\tfrac12\log(1+e^{-16})+4-\sqrt{15}
 <-0.00874170918.
```

The independent checker reproduces this strict sign using rational
logarithm bounds, an integer-squared lower bound on sqrt(15), and the
first 50 positive exponential-series terms. All approximation errors
vanish as epsilon decreases after m increases. Hence a fixed sufficiently
small epsilon_0>0, and then a fixed negative margin, are genuinely
licensed. The theorem is not restricted to widths shrinking with m.

## 5. Scope boundary

This removes a stable near-flat class for Walsh bases on the declared
high-probability selector events at the common tilt t=4. The
complementary profile partition sum remains an independent obligation.
Class-specific tilts cannot simply be combined in one Finner estimate.

The separate randomized recursive-Hadamard ensemble need not preserve
Walsh affine-coordinate identities. Thus the eight-point theorem cannot
be imported unchanged into that ensemble without a new transfer or an
appropriate terminal-level coding argument. This audit makes no such
transfer and no claim of a new universal cap or convergence result.

## 6. Diffuse finite bins and the sharp scope of the counting relaxation

The complete companion
`continued_feedback_diffuse_bin_count_and_limit_2026_09_06.md` was
independently read. Its finite-bin extension and all-tilt relaxation
obstruction both pass, with their stated diffuse and Walsh hypotheses.

Uniform U2 diffuseness of the supported spins, together with selector
quasirandomness, supplies exactly the same four-wise posterior recovery
cost. It does not assert eight-wise independence. Quantizing the selected
half of the spectrum replaces its sign-only code by the product symbol
code of rate `[H(nu)+(1-nu_0)log2]/2`. The zero bin requires no sign bit:
replacing its possibly nonzero coefficients by zero is already included
in the squared quantization error. Empirical nonempty bins have mass at
least 1/m, so the exceptional DC coefficient costs only O(log m). The
partition is selected to minimize the SUM of this symbol cost and the
conditional recovery cost, which is justified by their joint average.

Quantization error d^2 m leads to at most 32d^2 k corrected integer
symbols among m/2 outputs. Thus the correction fraction is 64pd^2 and
the half-rate correction term is `[h(e)+e log17]/2`. Its use for an
at-most-e fraction is monotone precisely in the declared range e<17/18.
Conditioning, fixed-bin profile enumeration, and the trivial 2^k cap
therefore give the claimed sufficient upper expression.

The equal-mass trial amplitudes a=1/2 and b=sqrt(7)/2 have second moment
one. Every self-coupling has equal diagonal masses r/2; optimizing its
single parameter gives exactly

```math
\Psi_K=\log\bigl[(\sqrt{K_{aa}K_{bb}}+K_{ab})/2\bigr].
```

Using K_aa,K_bb>=1/2 and K_ab>=exp(-gamma t)/2 cancels the symbol code's
log2 term and yields the lower test
`C0+log(1+exp(-gamma t))/2+Delta t` for the relaxed upper expression.
Here gamma=2-sqrt(7)/2 and Delta=1-sqrt(15)/4. Its minimum is
`C0+h(2Delta/gamma)/2`, with
`2Delta/gamma=(4-sqrt(15))/(4-sqrt(7))` strictly between .09 and 1/2.
Since C0>-.136 and h(.09)/2>.151, the expression is >.015 for every
positive tilt. This is an obstruction to this PARTICULAR counting and
permanent relaxation, not a lower count or feasibility assertion for
the trial inverse-Walsh profile. Arbitrarily fine bins do not remove
the obstruction because the trial itself has only two atoms.
