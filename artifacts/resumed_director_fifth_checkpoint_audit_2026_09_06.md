# Fifth checkpoint: causal-shape duality, critical scalar recovery, and nonlocal scope

Date: 2026-09-06, approximately 05:00 UTC. This is a continuing checkpoint,
not completion of the six-hour campaign. The director reconstructed the
proofs below; numerical experiments are distinguished from certificates.

## 1. New original-problem lower bound

The same actual 21-anchor core permits any even unit feedforward inverse
u(V), with covariance <E[g|V],u>. This is not a scalar feedback equation
for the conditional inverse, whose norm and derivative energy differ from
those of g. The three exactly orthonormal shapes are a central indicator,
its centered conditional inverse, and the tail conditional inverse.
The chosen coefficients satisfy unit norm identically, with rational
nonleading coordinates and one exactly enclosed square root.

I separately checked the Gaussian first-chaos coefficients, the inverse
identity K=lambda*g+gamma*u, the residual variance, and conditional Jensen.
I read the entire new certificate and replayed it against the rational
target 539339233892013/1250000000000000. The result is

    liminf M_n/n^(3/2) >= .4314713871136104.                 (1)

The exact lower endpoint is
.431471387113610406780573650606255124749081686370408082604706.
All central second moments use a finite alternating Gaussian series with
an outward remainder; all first moments use Hermite endpoint identities.
The final 64-bin Jensen step is a lower bound, not numerical quadrature.
The inherited finite-approximation-before-n and spectral-deletion-last
orders are unchanged. The upper remains .5; convergence is unresolved.

Files:

* `resumed_response_causal_inverse_shape_dual_2026_09_06.md`;
* `computations/resumed_response_feedforward_inverse_certificate_2026_09_06.py`;
* its exact JSON in `computations/results/`.

## 2. An exact optimization reduction, not another shape guess

At fixed mask and output covariance, the ENTIRE unit causal inverse class
has an exact two-scalar dual. The director derived it independently using
Gaussian absolute-value conjugacy, a Hilbert-sphere support function, and
the atomless three-moment closure. Every maximizing triple has strictly
positive residual moment, so the supporting dual attains finite parameters.
The verifier independently reconstructed the proof and the exact change
of variables to the response researcher's energy/covariance dual.

This is a genuine finite-dimensional optimization reduction for the
declared conditional response class. It does not optimize the rich core,
the covariance, or the mask and does not remove conditional Jensen loss.
It supplies neither an upper bound on actual signings nor convergence.

* `resumed_director_feedforward_shape_duality_2026_09_06.md`;
* `resumed_bound_audit_feedforward_shape_dual_2026_09_06.md`.

The exact all-shape upper certificate at the displayed fixed mask/covariance
uses strong-concavity residual bounds and convex chords in the conditional
polynomial. I read and replayed it independently: its upper endpoint is
.43147220696185809, giving a certified gap below 8.2e-7 from (1) for this
declared slice. The Mehler second-derivative bound, strong-concavity residual
bound, and chord integration were checked separately. This is NOT a global
upper bound on the response class or actual signings.

## 3. Whole scalar class: attained optimum, criticality, and uniform escape

I reconstructed Sections 1--10 of
`resumed_response_scalar_sobolev_full_center_variation_2026_09_06.md`.
Even unit inverse functions with Gaussian derivative energy at most one
form a strongly compact class; the exact scalar height theorem realizes
every member, including the critical boundary. The full-center objective
attains its maximum and normalized degree truncations approximate it at
O(N^(-1/2)) in objective on the compact high-value region.

The nonzero-covariance maximizing inverse necessarily has active derivative
constraint and a positive OU multiplier. The zero-covariance branch has
value below .409 and cannot contain the known maximizer. Thus the true
scalar optimum satisfies a genuine critical free-boundary variational
equation, not an unproved subcritical stationary condition.

Height iterates satisfy deficit at most 2/(h+2). For two inverse shapes,
their actual Gaussian coordinates satisfy

    ||V_g-V_h||_2 <= sqrt(2 ||g-h||_2).

This gives one common finite approximation bank for the compact scalar
class. Combined with full-gradient escape it proves a uniform strict
improvement above the ENTIRE attained scalar full-center supremum. That
increment is not evaluated, and it is not uniform over unrestricted
growing matrix-valued cores.

The scalar degree-320 and one-birth numerical tests remain diagnostics;
they are below (1), so extra degree tuning was stopped. Their scripts and
saved outputs make that negative experimental decision reproducible.

## 4. Nonlocal tests: positive realization and precise obstructions

The enlarged model allowing arbitrary Fourier-index permutations attains
the PSD-majorant norm T(B), by empirical Gaussian-cloud matching. This
does NOT assert that the allowed outer signings realize those permutations.
The director reconstructed its covariance/pseudoinverse and singular
covariance cases. A separate deterministic theorem shows that a single
dense bilinear compatibility condition on index pairs still permits
asymptotically exact cloud matching. Thus that local algebraic condition
cannot honestly be used as a broad impossibility claim.

For quadratic affine bent families, however, the entire index graph is
defined by polynomially many bilinear coefficients, even after arbitrary
physical ancillas. The resulting index alphabet is subexponential in the
number of frequencies. An independent-cloud concentration theorem shows
that such alphabets have positive Gaussian matching error, even after
arbitrary common frequency signs. This is not a statement about correlated
profiles, every nonlinear bent family, or actual signing minima.

Actual affine MM constructions provide Toffoli-type and higher polynomial
index examples. Exact finite cubic phase modules and a cubic controlled
index are also reconstructed and integer-checked. A separate simultaneous
tensor quotient already realizes all binary Fourier-sign multipliers as
limiting tests. None of these individual statements permits arbitrary
overlapping gate composition. A new phase-to-index selector mechanism is
being tested after this checkpoint; it is not a banked conclusion here.

Key files:

* `resumed_convergence_fourier_permutation_majorant_2026_09_06.md`;
* `resumed_director_bilinear_constrained_matching_2026_09_06.md`;
* `resumed_bound_audit_quadratic_index_information_obstruction_2026_09_06.md`;
* `resumed_director_phase_robust_cloud_obstruction_2026_09_06.md`;
* `resumed_convergence_affine_bent_index_modules_2026_09_06.md`;
* `resumed_convergence_cubic_affine_mm_modules_2026_09_06.md`.

## 5. Actual signings: a useful rigidity statement with narrow hypotheses

If a symmetric full Hadamard H has a complete orthogonal Boolean
eigenbasis, any hollow signing A differing on d unordered edges obeys

    Q(A)/n^(3/2) >= 1/2-1/(2n)-2d/n^2.

I reconstructed this by averaging its oriented eigenvectors. Consequently
a genuinely sub-.5 minimizing sequence cannot be o(n^2)-edit close to
this class. Tensor powers of J4-2I supply an infinite example class.
No such eigenbasis is asserted for every conference/Hadamard matrix.
The eigenvector law has an O(n), not O(sqrt(n)), cap window, so it does
not close the critical insertion-law obligation.

`resumed_bound_audit_boolean_eigenbasis_hamming_rigidity_2026_09_06.md`.

## 6. Continuing judgment

The response method continues to produce actual lower certificates and
now admits a global conditional-shape dual. Its fixed-rule ceiling still
prevents treating constant improvement as a proof of sharpness. The
independent nonlocal work must establish actual compatible realization,
not merely the availability of separate gates. No convergence proof,
nonconvergence proof, or comparable-order recurrence is claimed.

Next discriminating work: a correctly conditioned rich-core response
birth; the phase-to-index realization mechanism and its composition scope;
and independent auditing of any resulting original-problem implication.
Research continues through the remaining campaign time.
