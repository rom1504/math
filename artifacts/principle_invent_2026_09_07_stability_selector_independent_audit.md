# Independent reconstruction of the stability-plus-selector signing improvement

2026-09-07, approximately 17:48 UTC. **PASS.** I reconstructed the two
director proofs and their composition, rather than relying on their status
labels. The conclusion is a strict actual all-order upper improvement by
a positive symbolic amount; it does not change a reported decimal or prove
convergence.

Audited sources:

- `principle_director_stability_orthant_penalty_2026_09_07.md`;
- `principle_director_selector_entropy_stability_2026_09_07.md`;
- the exact finite normalization in
  `principle_invent_2026_09_07_stability_aware_weave.md`;
- the deterministic transfer and global split in
  `principle_invent_2026_09_07_stability_global_dichotomy.md`.

## 1. Imported concentration theorem checked at its primary source

I checked Theorem 4.1.1, Lemma 4.1.2, and the coordinate induction in
[Talagrand's primary paper, Section 4.1](https://arxiv.org/pdf/math/9406212).
The exponential moment has factor 1/4 exactly. The proof conditions on
one coordinate and integrates its marginal; the same induction applies
without requiring identical marginal measures. Thus the actual biased
independent edge signs are covered.

For a nonnegative convex K-Lipschitz function L, a subgradient at z gives
L(z)<=2K times the convex distance from z to {L=0}. Jensen consequently
gives P(L=0)<=exp[-(EL)^2/(16K^2)]. This uses the convex distance to the
zero set directly, not an invalid concentration assertion for -L.

## 2. Conditional row orthant bound reconstructed

Condition arbitrarily on all incoming coordinates of magnitude greater
than V. The conditional mean input has coefficient magnitudes at most
the original incoming magnitudes, so orthogonal contraction bounds its
physical l2 norm by sqrt(Cm). The deterministic physical diagonal adds
at most D sqrt(k). At least k/2 physical coordinates therefore have mean
bounded in absolute value by the stated M.

The retained Hadamard entries all have magnitude one, so every physical
coordinate has EXACTLY the same remaining variance in [kappa,C]. For
the light signs, each centered summand is bounded by 2V/sqrt(m), and
the sum of third absolute moments, including matched Gaussians, is
O(VC/sqrt(m)). Smoothing the hinge with fixed smoothing parameter and
then applying second-order independent Lindeberg telescoping supplies
a uniform expectation comparison. The smoothing parameter is fixed
before m tends to infinity; it may depend on the small positive g.

The Gaussian negative-part expectation decreases with the mean and
increases with the standard deviation, so the stated g is a valid
strictly positive lower bound. Summing over the k/2 selected physical
coordinates gives EL>=gk/4. The map from the independent light signs
to all fields has operator norm at most V, hence L has Lipschitz norm
at most V sqrt(k). The probability exponent g^2k/(256V^2) follows.
All these estimates are uniform over the conditioned heavy signs.

The weave substitution has d_a equal to the stated full within-fibre
rank-one contribution minus the physical diagonal term. The bound
|u_i(i)|<=sqrt(k) gives |d_a|<=2. This remains true for maximally
concentrated diagonal spectra. Deleting column i only decreases the
contraction norm. The edge tilt is tanh(2t sigma u_i(j)u_j(i)), with no
missing factor two.

## 3. Selector count reconstructed without a probabilistic spectral claim

A concentrated spectrum lies within sqrt(epsilon m) of a coordinate
subspace of dimension at most m/V^2. Nets on the radius-sqrt(m) ball
at scale sqrt(epsilon m) have size (1+2/sqrt(epsilon))^d. There is no
extra log m in this dimension-normalized covering number.

After pulling a net point back through the arbitrary orthogonal basis,
nearest-ternary rounding is well-defined with arbitrary deterministic
ties. Every coordinate disagreement from a target constant-type word
costs at least 1/(4p_m) squared distance. The total discrepancy is at
most 16p_m epsilon m, and replacing p_m by one only enlarges the
Hamming-ball bound. The number of alternative symbols is two per
changed coordinate. Requiring exactly k nonzeros can only decrease the
count.

Uniform selector averaging divides this constant-type count by precisely
binom(m,k); there is no missing 2^k, because spin summation is already
included in the count of ternary words. Each deleted-row orbit factor
is at most one. Thus the row exponent is -h(p)+xi, uniformly over the
entire basis, including arbitrary terminal orthogonals.

## 4. Global composition reconstructed

Diffuse row mass forces a positive fraction of incoming columns with
bounded energy and positive light tilted variance by the deterministic
reverse-heavy and high-column truncations. Their stability factors can
be pulled out BEFORE any output-column permutation average. All other
terms are then handled by the old graph-CS/orbit theorem unchanged.

In the complementary sector, full-row concentrated indicators are
invariant under output-coordinate permutations. They survive the same
conditional permutation/diagonal argument, after which independent
bases and selectors give a product of restricted row partitions. The
sum over assignments of the two row classes costs only exp(O(m)).

The exact full physical diagonal combines with the folded kernel into
a factor at most exp(tm^2+tm), so it contributes no new m^2 loss. This
is stronger than assuming a coordinate is diffuse and dropping it.

For finite depth, use the COMMON unrestricted row upper a_cert+zeta,
where zeta is chosen smaller than the already reserved strict saving.
Then choose one fixed depth realizing that upper. This avoids assuming
the limiting exponent a_cert is attained at finite depth. No lower
bound on E_t is needed for improvement over a certified upper: only
a_cert>-h(p) is required. To improve over the exact E_t expression,
the analogous exact gap must be positive.

## 5. Independent scalar check

I independently evaluated mpmath interval expressions at 80 digits, not
by reading the stored numerical outputs alone. For p=24/25,
epsilon=1/100000, V=100, the outward intervals imply

    h(p)+p log 2-5151/6250 > .009,
    xi(epsilon,V) < .004,
    2 xi(epsilon,V) < h(p)+p log 2-5151/6250.

Their approximate values are .00920544107172045138 and
.00333549264439380467. Exact rational checks also give
U^2=C=10^10 >=4V^2/(delta epsilon), for delta=1/2,
theta=1/16000000000, and
kappa=sech(97000000)^2/1600000.

The positivity of the final improvement is analytic. Its fantastically
small magnitude is not estimated by underflowing floating arithmetic.
There is no claim of a changed decimal endpoint. The valid theorem is
strict improvement below the old exact certified expression by some
fixed positive amount, using actual full signings at all orders.

## 6. Stronger no-bad-row conditioning principle: independently PASS

The director subsequently observed a stronger consequence of the same
selector count. For any fixed p in (0,1), choose epsilon,V with xi<h(p).
For EVERY fixed orthogonal basis O, the expected NUMBER of concentrated
row-spin words under a uniform selector T is at most exp(-c m), for a
fixed c>0 and all sufficiently large m. This is integer-valued, so

    P_T(no concentrated word exists | O)>=1-exp(-c m),

uniformly over O. Condition each independently sampled pair (O,T) on
this event. Its density relative to the unconditioned joint pair law
is at most 1/(1-exp(-cm)); consequently every nonnegative row partition
expectation grows by at most that factor. Its logarithm is negligible
even after multiplication over m fibres.

The event is invariant under arbitrary output signed-coordinate
permutations: it quantifies every row word, and concentration depends
only on spectrum magnitudes. Hence the output permutations remain
conditionally uniform, and independent per-fibre conditioning preserves
the product structure needed for graph-CS and the row theorem. There
is no assumption that a selector fixed before the basis works uniformly;
the actual construction is allowed to choose it from this conditioned
joint distribution.

Under the conditioned construction EVERY Boolean configuration has ALL
rows diffuse. Apply the deterministic transfer with delta=1 and use the
orthant/Finner loss uniformly. The bad-sector argument is no longer
needed, and neither is a+h(p)>0. Thus every finite all-vector row moment
certificate of the audited form is strictly improvable for every FIXED
p<1 and finite t>0. The old numerical endpoint is a special case.
The p=1 endpoint is excluded: h(1)=0 gives no selector entropy.

This strengthening is an actual finite-alphabet uncertainty principle,
not just a pressure-weighted distinction. Its proof uses the very same
uniform sparse-spectrum net count already audited above.
