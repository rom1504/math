# Scalar regression and the second-query gain, with ordered cutoffs

Date: 2026-09-06. Proposed bounded-response assembly, submitted for
independent audit. The fixed-polynomial contractions and scalar Stein
regression have been separately reconstructed; the ordered approximation
passage is written explicitly below. No positive innovation bound or new
universal constant is asserted.

## 1. Actual first-history field, not a surrogate

Use the smooth bounded scope of
`continued_feedback_first_marked_history_energy_projection_2026_09_06.md`.
Thus B is an actual normalized hollow signing with bounded operator
norm, and

`W=(S,G,Y,QS,QD)`, `G=BS`, `D=S h2(G)`, `Y=BD`,
`f(G,Y)=b0 G+b1 Y+r(G,Y)`,
`V=b0 QS+b1 QD`, `Z=B r(G,Y)`,
`C=H(G,Y) psi(V+Z)`.

The source f is bounded, odd, and Gaussian-a.e.-continuous as stipulated
there. The mask H is bounded, even, and has the stipulated actual-marginal
L2 approximation. The odd response psi is bounded by one, with bounded
first and second derivatives. Assume `0<=H<=1-|f|`, so the actual first
endpoint `X=f+C` satisfies `|X|<=1`.
Let c0 and the deterministic a be the actual coherent conditional mean
and expected derivative from that theorem, using its ideal old-noise
variance. Define EXACTLY

`L=B c0+B D_a Z`, `eta=B C-L`.                               (1)

Keep this definition of L. The actual endpoint return is
`BX=V+Z+L+eta`; the already retained old shift `V+Z` is not absorbed
into L or into its covariance with eta.

No law is assigned to L. All these fields are odd and centered. Operator
norm bounds and the earlier averaged L2 estimates give a constant M with

`(1/n) sum_i E[L_i^2+eta_i^2] <= M+o(1)`.                    (2)

The first-history proof supplies finite polynomial approximants to the
triple `(C,L,eta)` with arbitrarily small averaged L2 error, in the order
`fixed approximation -> n limit -> approximation limit`. Their exact
squarefree main obeys the quantitative collision estimates already used
in the marked energy theorem.

## 2. The finite-polynomial scalar Stein statement

At one fixed polynomial stage write `eta=sum_p eta_p`. Its surviving
Walsh degrees satisfy p>=5, and all proper cuts are small. The now-proved
full-contraction criterion against every higher-degree component of L
is in `continued_audit_boundary_graph_and_next_return_2026_09_06.md`.
The same source argument applies to old Z, with its corresponding bounded
outer transport. The equal-degree eta/old-Z covariances have vanishing
AVERAGE ABSOLUTE value, by the first-history source pairing bounds with
arbitrary bounded diagonal tests. They are used here, not assumed zero.

Let `U=sum_p eta_p/p`, the inverse-generator field. Exact Boolean
integration by parts gives

`E[eta F]=sum_a E[Delta_a U Delta_a F]`.

Every unequal-degree derivative-product coefficient is controlled by
a proper contraction or the new full-contraction criterion. The scalar
constants are exactly

`E Gamma(eta,eta)=E eta^2=sigma^2`,
`E Gamma(eta,L)=E eta L=rho`.

Hence their centered fluctuations tend to zero in averaged L2. The
maximum noise influence tends to zero as well. Fixed-degree
hypercontractivity and the exponential finite-difference remainder
therefore extend the self-contained proof in
`continued_audit_boolean_stable_noise_stein_2026_09_06.md` to this
inhomogeneous scalar eta and the literal L.

On rows with `sigma^2>=epsilon` and `E L^2+sigma^2<=K`, set

`lambda=rho/sigma^2`, `R=L-lambda eta`.

Here `|lambda|<=sqrt(K/epsilon)` and `E R^2<=E L^2`. The mixed Stein
covariance with R vanishes. Thus eta is stably compared to an independent
`N(0,sigma^2)`, while keeping the ACTUAL joint law of `(R,W,old Z)`.
This is after averaging rows; no unsupported uniform row statement is
needed. In particular old X, a function of `(W,old Z)`, is retained.

The independent Gaussian term in `BX=V+Z+L+eta` consequently has coefficient

`t_i=sigma_i+rho_i/sigma_i`.                                 (3)

This coefficient may have either sign. Its absolute value is the
innovation standard deviation in the scalar-regression comparison.
Aggregating is legitimate; it merely gives a weaker Jensen innovation
floor than regressing separately on every homogeneous noise component.
The comparison before Jensen keeps the other Gaussian directions inside
the literal regression residual R.

## 3. Bounded-response passage at fixed variance cutoffs

All quantities in this section refer to the actual bounded field (1).
Define deterministic row statistics

`sigma_i^2=E eta_i^2`, `rho_i=E eta_i L_i`, `v_i=E L_i^2`,

and row sets

`G_(epsilon,K)={i: sigma_i^2>=epsilon, v_i+sigma_i^2<=K}`,
`A_(epsilon,K)={i: sigma_i^2<epsilon, v_i+sigma_i^2<=K}`.

On G define R and t by (3). Put `h_i=1-|X_i|` and retain the old shift
`T_i=V_i+Z_i` literally.
For every fixed epsilon>0 and finite K, the proposed comparison is

`(1/n) sum_(i in G) E[h_i |(BX)_i|]`
` - (1/n) sum_(i in G) E[h_i E_N |T_i+R_i+t_i N|] -> 0`.      (4)

Here and below E_N acts only on a fresh independent scalar Gaussian.

The passage from the fixed polynomial result to (4) uses these specific
orders and bounds.

1. Fix epsilon and K before selecting an approximation. Averaged L2
   approximation of L and eta implies averaged L1 approximation of
   their variance and covariance statistics, by Cauchy--Schwarz.
2. Discard the deterministic rows where any of these statistic errors
   exceeds a fixed small threshold. Their fraction tends to zero with
   the approximation error. On retained good rows the approximating
   variance is at least epsilon/2 and its second moments are bounded
   by K+1. The rational maps defining lambda and t are uniformly
   continuous and bounded there.
3. This converts coefficient errors into weighted L2 errors of R:
   first use a uniform bound for lambda on the retained rows, then the
   averaged L2 field error; for coefficient differences use their uniform
   smallness times the bounded row second moments. No merely averaged
   coefficient convergence is multiplied by an uncontrolled raw square.
4. The factor `h=1-|X|` is bounded by one and is Lipschitz in the bounded
   X coordinate. For polynomial approximants it may be replaced by
   `(1-|X_polynomial|)_+`, another bounded Lipschitz function. Its error
   is controlled by the earlier averaged L2 approximation.
5. Truncate the absolute-value test first. Stable characteristic-function
   comparison applies to the fixed finite old polynomial coordinates and
   R, with deterministic coefficients in a compact set. Their bounded
   second moments give uniform integrability for tests of linear growth.
   The discarded-row contribution is controlled by Cauchy--Schwarz.

Take the n limit at the fixed polynomial stage, then remove the
approximation and statistic-error thresholds. Only after (4) is obtained
does one send epsilon down to zero or K to infinity. No positive
pointwise variance floor is used.

## 4. Full gain formula with harmless low-variance rows kept literally

Let

`J_n=(1/n) E sum_i (1-|X_i|) |(BX)_i|`.

The low-variance set A obeys

`(1/n) sum_(i in A) E h_i ||T_i+L_i+eta_i|-|T_i+L_i|| <= sqrt(epsilon)`.

The rows with `v_i+sigma_i^2>K` have fraction at most `(M+o(1))/K`.
The average squared norm of BX is at most `||B||^2`, because `|X|<=1`.
Their gain contribution is therefore at most
`||B|| sqrt((M+o(1))/K)`. Combining this with (4) gives

`J_n = (1/n) sum_(i in G) E[h_i E_N|T_i+R_i+t_i N|]`
`      + (1/n) sum_(i in A) E[h_i |T_i+L_i|]`
`      + O(sqrt(epsilon)+sqrt(M/K)) + o_n(1)`.                (5)

Equation (5) is a cutoff formula, not an assertion that the raw regression
coefficient is bounded on all rows. It includes the case of identically
zero new variance: those rows simply keep T+L. No hard-threshold limit is
taken inside the OLD construction of C.

## 5. Operational extremal-energy consequence

For each seed realization the two vectors

`x_plus=X+(1-|X|) sign(BX)`,
`x_minus=X-(1-|X|) sign(BX)`

lie in the continuous cube. Because B is hollow, multilinear rounding
shows that its Boolean maximum dominates their absolute energies. The
difference of their two energies, divided by two, is exactly

`(1/n) sum_i (1-|X_i|)|(BX)_i|`.

Thus the finite normalized extremal value for the actual signing,
`max_x |x^T Bx|/(2n)=Qabs(A)/(n sqrt(n-1))`, is at least J_n.
In particular, after the stated cutoff limits, (5) is an actual
second-query lower certificate, retaining the coherent return.

A simpler but weaker consequence of (4), since
`E_N|T+R+tN|>=sqrt(2/pi)|t|`, is

`liminf J_n >= liminf (sqrt(2/pi)/n)
     sum_(i in G_(epsilon,K)) E(1-|X_i|)
                       |sigma_i+rho_i/sigma_i|`.             (6)

No claim that the right side improves 0.4333221116640807 is made.
It can vanish, and the new theorem does not provide a lower bound on
its innovation or on the remaining slack. The original convergence
question is unchanged.
