# Third-checkpoint independent reconstruction

Date: 2026-09-06. These are proof checks, not votes by previous auditors.
The original convergence question remains open.

## 1. Full Gaussian functional: uniform quantitative separation

I reconstructed Section 6 of
`resumed_response_full_center_gain_and_purification_2026_09_06.md`.
For an original mask, purification preserves its mass mu, first projection,
adjoint K, and signed mean J, while the residual variance becomes
1-mu-||P1F||^2. The Gaussian rearrangement bound on that projection uses
the actual direction P1F/||P1F||, without assuming independence from the
support. Thus J<=sqrt(mu) 2phi(a), mu=2Phi(a)-1. The scalar cap is strictly
unimodal, and the exact endpoint evaluations confine every J>=.43 to
.47<a<.88. Signed Jensen and the monotonicities of the Gaussian absolute
gain give the stated binwise lower bound.

I read and independently ran the 128-bin rational program. All bins have
positive residual variance and gain above 1/40000; the lowest rigorous
endpoint is .000026760544759910564903901171555280747476109377036807622545.
Taking a sequence approaching the old marked-mask supremum proves

    liminf M_n/n^(3/2) >= C_mark+1/40000.

This is uniform over that certificate class, not a repeated fresh-spin
theorem. It does not assert a positive fixed improvement over every member
of the enlarged full-response class.

I also independently replayed the revised rational candidate with
alpha=3623/5000, resolvent=3479/1000 and the same 21 anchors. Its derivative
energy is below .999904932505320246<1. The full signed-Jensen bound is
strictly above .431124492502981658479186277525917256229651077680764381299370.
The smaller contraction margin remains strictly positive and is fixed
BEFORE finite approximation and the matrix limit. It is not a uniform
rate as these numerical parameters approach a contraction boundary.

## 2. Finite canonical coordinate birth

I reconstructed the exact update in
`resumed_response_coordinate_birth_escape_2026_09_06.md`. For a finite
ancestor-closed bank, a nonzero bounded function H Gamma'(K) cannot have
only finitely many Hermite coefficients unless it is constant. Its
vanishing on the positive-mass occupied set forces that constant to zero.
Above the pure-noise value sqrt(8/(27pi)), the function is nonzero.
Hence an unused creation coordinate with old-measurable inverse feature
and nonzero response correlation exists.

Filling its far Gaussian tail inside H preserves all old first-chaos
coefficients. The new one is mu z; the variance increment is
mu p-mu^2 z^2>=0 by Cauchy--Schwarz. Since K' remains old-measurable,
the new mask factors exactly, and convexity gives

    C'-C >= p[(1-p)mu(z/p)|c|-C].

A sufficiently distant but FIXED tail makes this positive. No matrix
response was declared to be an independent input. The infinite sequence
has nested occupied sets and hence an actual L2 limit, but its increments
have no uniform lower bound. The no-finite-canonical-maximizer statement
does not by itself cover arbitrary finite cyclic Gaussian frames. The
separate residual criterion in Sections 6--7 DOES cover the actual cyclic
polynomial core plus its nonpolynomial H feature: a polynomial vanishing
on the exterior open set is zero, and a nonconstant polynomial g cannot
be constant throughout the interior strip. This proves a further strict
but unevaluated escape above that precise full-response construction.

## 3. Causal-frame stationarity, including a necessary correction

I independently derived the first variation and the distinction between
an old-frame gradient and its independent innovation. I then reconstructed
all cases of `resumed_bound_audit_triangular_stationarity_tie_obstruction_2026_09_06.md`.
The largest nonzero projection index and largest gradient index must be
treated separately. Conditional Gaussian threshold differentiation handles
the latter index at least as large as the former; earlier measurability
handles the other case. This excludes a nonzero stationary gradient in
a finite strictly causal frame at positive value.

The initial conjecture that zero gradient implies the pure-noise cap was
WRONG: ties at K=0 allow a larger value. The corrected proof first shows
B=Psi(0,t). In particular, if B>Psi(0,t), the zero-gradient equation and
the defining equation for B would give an expectation of
2z phi(z)/(2Phi(z)-1)<1 equal to one. On nonzero K the support is then
empty and |K| is constant. The resulting scalar equations, together with

    (2Phi(z)-1)^2 <= (4/pi)(1-exp(-z^2/2)),

give p<=1/(1+4t) and C<=phi(0). I checked the integral derivative proving
this inequality and the rational polynomial comparison. Thus no finite
strictly causal frame can be stationary above 1/sqrt(2pi). Cyclic
fixed-point frames are explicitly outside this claim.

## 4. Original-class local-algorithm ceiling

I reread both
`resumed_convergence_involution_local_ceiling_2026_09_06.md` and
`resumed_bound_audit_gfom_ceiling_handoff_2026_09_06.md` against the previously
inspected primary WZF hypotheses. Exact flat involutions pass the uniform
matrix-power tests. Signed-permutation equivariance removes the auxiliary
random conjugation, not arbitrary coordinatewise symmetry. A fresh Gaussian
input before each multiplication repairs finite covariance degeneracy;
the zero first spectral cumulant avoids a current-input Onsager term.
Fixed Lipschitz recursions are stable as that perturbation vanishes.

The finite-tree handoff first proves the Boolean own-input-free Wick error,
then truncates parent nonlinearities and chooses child accuracies to match
their fixed Lipschitz constants. It does not propagate an L2 error through
an unbounded polynomial without control. The final sign is replaced in
ENERGY by a soft sign, not in L2 across a possible atom at zero.

The Haar beta tail and union bound yield sqrt(15)/8 for the normalized
half-energy. Removing the small diagonal handles cube-valued outputs.
Consequently the entire full-response variational supremum, including its
countable closure, is at most sqrt(15)/8. This ceiling applies to the
audited fixed-rule class; it is NOT an upper bound on actual Hadamard
ground states, arbitrary signings, or algorithms whose complexity grows
with matrix order.

## 5. Nonlocal three-cell realization

I read `resumed_convergence_h144_three_cell_2026_09_06.md`, reconstructed the
transpose/gauge equivalence to the prescribed H144, and independently ran
its integer verifier. The first seed has quadratic energy 9792; the second
has bilinear energy 9792 and R4-lift energy 78336. Both recover R=T=17/6.
The simplex-code quotient is a separate scalable theorem: every nonzero
character occupies 3^(r-1) factors, an odd number, so the product action on
the quotient is the SINGLE averaging reflection. Grouping its equal atoms
then approximates arbitrary finite probability weights. Tensor products of
these tests are valid. Noncommuting composition on the same interface, and
arbitrary polar-operator realization, are not consequences.

## 6. Conditional projection and the improved exact decimal

I independently derived the projection of the finite-resolvent inverse
feature onto V before receiving the completed certificate. Multinomial
Hermite contraction gives the scalar coefficients g_d in
`resumed_response_conditional_v_full_center_2026_09_06.md`; deleting the
anchor coefficients subtracts B_d/a and B_d/a^2 in the first- and
second-resolvent sums, respectively. The first projection of F gives
K=lambda g+gamma H. Conditional Jensen, first given V and then on
64 symmetric rational bins inside the central strip, is in the correct
LOWER-bound direction.

I read the entire rational verifier and reran it independently. The
normalized Hermite integral identity uses exact endpoint terms, so there
is no numerical quadrature error. The actual inverse feature has degree
200, so there is no unbounded Hermite tail to discard. The implementation
forms large normalized Hermite coefficients by taking square roots of
exact rational squares before interval evaluation; it does not invert an
interval for 1/sqrt(200!) which contains zero at the working grid.

The replay verifies the derivative energy below .999904932505320246,
all positive bin masses, the total first moment identity, and the bound

    liminf M_n/n^(3/2)
      >= .431460392823700540631522256404469656378004801252582084430582.

The safely rounded public endpoint is .4314603928237005. Conditioning
retains genuine center variation which the previous single Jensen step
discarded. No high-dimensional sampling or unexplained decimal is used.

## 7. Explicit ascent rather than an assumed full update

I checked the first variation and the second-derivative estimate (13) in
`resumed_bound_audit_variational_stationarity_dual_2026_09_06.md`.
On the feasible mixture, t_theta>=t/2. The two support-derivative terms,
three Hessian terms, and t_theta''=-2|a_new-a|^2 give exactly the displayed
M bound using only L2 inverse-feature norms. Hence a nonzero gradient
innovation gives a certified positive line-search step. Taking the full
step is NOT licensed; the saved numerical example in fact decreases.
This is a generative finite-stage theorem, not a uniform convergence rate
for the response optimizer or for the original signing sequence.
