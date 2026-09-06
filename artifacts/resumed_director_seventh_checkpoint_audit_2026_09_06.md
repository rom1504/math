# Seventh checkpoint: original-signing stability and the genuine feedback boundary

Date: 2026-09-06, approximately 06:45 UTC. This checkpoint does not end the
authorized campaign. The universal interval is still

    .4333221116640807 <= liminf M_n/n^(3/2)
                      <= limsup M_n/n^(3/2) <= .5.

Convergence and nonconvergence remain open. The statements below distinguish
original signings, restricted algorithmic limits, and numerical evidence.

## 1. A robust original-signing consequence of the Ramsey mechanism

Suppose a symmetric additive-convolution matrix A_q over an odd finite field
has all but o(q) Fourier eigenvalues of magnitude at least
(1-eta)sqrt(q), for every fixed eta>0. Every hollow symmetric B_q with
||B_q-A_q||_F=o(q) then has liminf Q(B_q)/q^(3/2)>=1/2.
No operator bound on either matrix is needed for this statement.

The director reconstructed the proof in
`resumed_bound_audit_cayley_hamming_stability_2026_09_06.md`:

- Fix the approximation accuracy before q and use its finite Ramsey host.
- Dilations avoiding deficient modes have proportion tending to one.
- Average the selected polynomial witnesses over those dilations and all
  translations. Their oriented covariance M has ||M||_F=O_epsilon(sqrt(q)).
  This remains true when the selected configuration depends on the dilation.
- Frobenius perturbations therefore cost o_epsilon(q^(3/2)).
- To pass from the polynomial witness to a Boolean spin for arbitrary B,
  use coordinatewise L2 control under translations and real Grothendieck.
  Hollow symmetry gives beta(B)<=4Q(B); no unproved spectral bound on B
  is inserted. The resulting factor tends to one as epsilon tends to zero.

For sign matrices an edit count d gives ||B-A||_F^2=8d. Consequently
o(q^2) edits cannot turn a spectrally near-flat additive-Cayley family into
a fixed sub-1/2 construction. This says nothing about all other signings.

## 2. Exact radial symmetry and a quantitative Paley exclusion radius

For a Paley core, averaging affine images of one oriented ground state gives
exactly

    E[sigma xx^T]=[2Q(A)/(q(q-1))] A.

Therefore every signing B differing on d edges satisfies

    Q(B)>=Q(A) |1-4d/(q(q-1))|.

The director checked every step, read the integer verifier, and replayed it
at q=5,13,17. The exact caps 4,20,32 and all affine covariance identities
pass. The asymptotic all-order Paley theorem then gives the explicit radius:
coefficient at most 1/2-delta requires at least (delta/2-o(1))q^2 edits
from every switched/relabelled Paley core and its negative. Border deletion
gives the same leading statement for the Paley conferences.

The director separately derived the sufficient symmetry behind the identity:
an unordered-edge-transitive signed-permutation group preserving A up to a
global sign gives the same off-diagonal radial law. The diagonal need not
vanish and is explicitly irrelevant for hollow targets. The independent
radial-symmetry audit passes. Averaging arbitrary minimizers over different
matrices does not create this law for one fixed matrix.

Proofs: `resumed_convergence_paley_radial_hamming_law_2026_09_06.md`,
`resumed_director_radial_symmetry_principle_2026_09_06.md`, and
`resumed_bound_audit_radial_symmetry_2026_09_06.md`.

## 3. The common endpoint energy has now been checked, not assumed

On normalized symmetric Hadamards and the explicitly stable hollow class,
every fixed canonical odd-tree terminal rule has

    Q_B(F)/n -> 0,   Q_B(H sign(BF))/n -> 0,
    Q_B(H sign(BF) +/- F)/n -> +/- J(F,H).

All limits are in L2, with the rule fixed before n. Thus the <.45 ceiling
on J really does limit this terminal architecture on these original
signings. It is not a ceiling on their Boolean cap or on all algorithms.

The director read the complete endpoint proof, reconstructed its normalized-
L2 tree-to-GFOM recursion, and checked the Haar traffic theorem directly in
the primary source. Doubled-tree counting gives the required old Gaussian
and marked identities. For universality, the additional fresh-noise
nonsingularity justification uses Fan's Assumption 4.2(d,e) and Theorem 4.3;
the director read those primary statements. The earlier module was updated
to avoid assuming a next-step state evolution in order to justify itself.

Exact matrix-sign parity and the deterministic W2 energy limit then force
both self-energies to vanish. The hard-sign limit uses the separately proved
positive higher-chaos variance and its small-ball bound. Mere Frobenius-near
involution structure, arbitrary growing depth, and arbitrary eigenbases
are outside the argument.

Proofs: `resumed_response_flat_involution_terminal_endpoint_2026_09_06.md`
and `resumed_bound_audit_flat_endpoint_primary_2026_09_06.md`.

## 4. A genuine departure from the old response state

Two inertial updates u_next=sign(Bu+alpha u), with alpha=1/2, already
produce a charge-odd half D_2 whose squared L2 distance from the COMPLETE
old canonical Gaussian sigma-field exceeds 1/102400. This is a statement
about the actual coupled response at B and -B, not a comparison of unrelated
one-root marginals. The director reconstructed the conditional Bernoulli
gate and exactly replayed every rational Gaussian enclosure.

The finite-depth Haar posterior gives a paired-query recursion: each fresh
query direction exposes one positive and one negative eigendirection, and
its transverse response becomes a fresh Gaussian. The finite longitudinal
beta component is retained before being shown to vanish. All three
researchers and the director independently reconstructed this mechanism.

It yields actual retained energy, not another final terminal certificate.
For every nondegenerate Boolean query with y=k+sigma N, its stability gap
g=E(|y|-u y)>0. A damped feasible update gains at least g^2/8; independent
Boolean rounding preserves normalized energy on the stated hollow class.
Neither the gain nor the innovation is uniform over growing depth.

The exact 40-step Sylvester replay at n=2^18 gives
30324543/67108864=.4518709033727645874... . The director read the source
and reproduced it using integer arithmetic. The asymptotic value near
.4516 is only floating Gaussian integration, not a certified new lower
bound. A finite value above .45 alone would not refute an asymptotic ceiling;
the separate two-step measurability theorem proves the architectural escape.

There is now a quantitative depth hierarchy, not merely a frame-specific
stationarity statement. Peeling the newest causal Gaussian pair shows that
the unexposed squared norm and the old-history gradient gap cannot both be
below an explicit eta_m>0 depending only on m pairs. Either an in-history
damped update improves actual energy, or rounding and one actual new query
does. Thus C_(m+1)>=C_m+Delta_m>C_m for the restricted fixed-query optima.
The director checked the conditional projections after deletion, the Gaussian
tail inverse modulus, and the explicit nested-exponential threshold. Delta_m
is not uniform in depth. A separate fresh-gated construction produces
arbitrarily small stability gaps in O(epsilon^-2) depth while remaining below
the broader fixed-GFOM ceiling, so deficit from 1/2 alone cannot force a gain.

Proofs: `resumed_response_mixed_charge_involution_feedback_2026_09_06.md`,
`resumed_convergence_haar_adaptive_query_kernel_2026_09_06.md`, and
`resumed_bound_audit_mixed_charge_feedback_2026_09_06.md`. The uniform
hierarchy is in `resumed_bound_audit_uniform_query_depth_hierarchy_2026_09_06.md`.

## 5. Direct spectral tests that retain all aliased harmonics

For any real symmetric finite-abelian convolution matrix, phase-averaged
character square waves give

    g(a)=(8/pi^2) sum_(m odd>=1) lambda(ma)/m^2.

Absolutely convergent Mobius inversion, including annihilated/zero modes,
gives ||lambda||_infty<=(3/2)||g||_infty and Q(A)>=n||A||op/3.
The director checked the inverse Euler product and the weighted-symbol
sharpness example. That sharpness is not a sign-kernel counterexample.

Combined with the finite Ramsey host, an actual sub-(1/2-delta) Cayley
signing must have op norm at most (3/2-3delta)sqrt(q) AND a uniformly
positive fraction of modes below (1-delta/2)sqrt(q). The latter proof
counts all dilations and has no hidden subsubsequence quantifier.

In characteristic three the exact formula is g(a)=(8lambda(a)+lambda(0))/9.
Parseval for the ACTUAL sign kernel gives the finite identity

    sum_(a!=0) g(a)^2
      =[64q(q-1)+(q-81)lambda(0)^2]/81.

Hence Q(A)>=(4/9)q^(3/2) for q>=81, without a row-sum or flat-spectrum
assumption. The director and convergence researcher then proved an explicit
improvement for this entire nonflat target class:

    liminf Q(A)/q^(3/2) >= 2080/(9 sqrt(269441))
                        = .4452346798594428... .

A projective plane's thirteen four-point lines cannot all have color split
3/1, since every point lies in four lines. A small blocking-set argument
strengthens this to 3N_mono+N_balanced>=3. Exact pair counting turns it into
asymptotic monochromatic-line density at least 1/26. Averaged Boolean nine-point
tests on those lines force a linear spectral-square deficit. Completing
one scalar square gives the constant. The director reconstructed the
proof and replayed all 8192 finite colorings and all Boolean nine-point
profiles. No spectral independence, large Ramsey threshold, or row-sum
assumption is used. This does not change the unrestricted .433322 bound.

Proof: `resumed_convergence_cayley_character_mobius_bound_2026_09_06.md`.
The incidence proof and exact program are indexed in
`resumed_convergence_characteristic_three_incidence_gain_2026_09_06.md`.

## Continuing work

The original interval has improved, and the nonlocal and feedback tracks
now have sharp, different scope boundaries. There is still no uniform
comparable-order recurrence, matching unrestricted upper construction,
or proof of two separated infinite optimal-order subsequences. Continue
the projective-incidence theorem and the quantitative finite-query escape
audit through the authorized endpoint, then state the remaining gap honestly.
