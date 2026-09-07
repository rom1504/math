# Independent audit: actual optimizer Gibbs covariance

Date: 2026-09-07. Verdict: PASS for the full canonical theorem
`decisive_director_actual_optimizer_covariance_2026_09_07.md`.
I reconstructed the score calculation before reading the final proof,
then checked all its quantifiers and constants against the saved version.

## 1. Exact row optimality, including the quenched field

Integrating x_i gives `Z=2 Z_cavity C_i(g)`, with
`C_i(g)=E_cavity cosh(sum_j J_ij x_j)`. The cavity law is independent
of all signs subsequently replaced on row i. Global optimality of
`E_g log Z`, not a pointwise-in-g optimum, therefore gives

`E_g log C_i(g)<=E_g E_row log C_i^row(g)
 <=(n-1)log cosh(beta/sqrt(n))<=beta^2/2`.

At each g, full conditional variance is
`v_i(g)=E_cavity sech(h_i)/C_i(g)>=C_i(g)^(-2)`.
Jensen on `exp(-2 log C_i)` then gives `E_g v_i(g)>=exp(-beta^2)`.
This step remains valid even when row fields are atypically large on
some disorder samples. It does not assert a pointwise uniform bound.

## 2. Score covariance and Schur complement

For `eta_i=x_i-s tanh(h_i)`, conditioning on x_i gives
`E eta_i x_j=delta_ij v_i` and `E eta_i^2=v_i`.
For i!=j, with `h_j=b+J_ij x_i`,

`E eta_i eta_j=-E[s sech^2(h_i) Delta_i tanh(h_j)]`.

The central-difference integral for Delta_i tanh, and
`|(sech^2)'|<=2`, give the explicit pointwise error

`|Delta_i tanh(b+J x_i)-J sech^2(b+J x_i)|<=2J^2`.

Hence the error matrix is symmetric, hollow, and has row absolute
sums at most `2 sum_j J_ij^2`. Gershgorin/Schur bounds it above by
`2 beta^2 I`. The leading term is
`-E[s V J V]<=E[V D V]<=D`, where V is a bounded diagonal random
matrix and the SAME deterministic D majorizes both J and -J.
No independence of V and J, or between coordinates of V, is needed.
Thus `E eta eta^T<=I+D+2 beta^2 I` exactly as stated.

The joint second-moment block matrix of x and eta is PSD. Increasing
its lower-right block to the positive diagonal upper bound preserves
PSD; its Schur complement yields the claimed diagonal covariance
lower bound. Global spin reversal gives E x=0, so second moment and
covariance agree. The Grothendieck majorant's factor four was checked
from the polarization and simultaneous SDP proof, not inferred from
separate one-sided SDPs. Jensen gives the diagonal certificate mass
in (2).

## 3. Flat bridges and the exact limitation

For independent optimized child laws, PSD order applied successively
to `tr(B Sigma_y B^T Sigma_x)` gives the weighted sum
`sum_ij B_ij^2 c_i c'_j`. For a flat sign bridge this factors into
the two certificate masses. This is a genuine non-annihilation result.
The nontrivial trace is the mass of the DIAGONAL PSD lower certificate;
`tr Cov(x)=n` by itself would be tautological.

Independence of the child laws is used in the trace formula. The common
orientation of an actual two-sided block endpoint is not silently
replaced by independent child orientations. The canonical note retains
that scope. A variance of order mn alone gives no order-(m+n) lower
bound on log moment-generating function at inverse scale sqrt(m+n),
as its final caveat correctly states.

## 4. A direct weighted-profile version

The SAME proof works at a global optimizer for any fixed nonnegative
edge magnitudes lambda_ij, with `J_ij=lambda_ij A_ij`. Define
`r_i=sum_j lambda_ij^2`. Row replacement gives
`v_i>=exp(-r_i)`. For a simultaneous majorant D>=+/-J, the symmetric
score error is bounded above by `diag(2r_i)`. Therefore

`Cov(x)>=diag(exp(-2r_i)/(1+d_i+2r_i))`.

This includes bounded-row-square weighted interpolation profiles.
It is still a LINEAR-spin covariance assertion; no quadratic-chaos
variance lower bound under a nonzero bridge tilt is proved here.

## 5. Finite reproduction

`computations/decisive_audit_gibbs_residual_covariance_2026_09_07.py`
enumerates all switching classes at orders 3 through 6 at beta
0.3, 0.9, and 1.8. It checks the score and Schur matrix inequalities
on every signing, and the row-optimality/variance bound on every
minimizer identified by exhaustive partition evaluation. It also
checks rectangular bridge variances and quenched-field versions at
order 4 by 64-node Gaussian quadrature.

These are finite floating-arithmetic diagnostics, not interval
certificates for transcendental partition minima. The theorem rests
on the exact inequalities above. The all-positive order-6 signing at
beta=1.8 supplies a checked NON-optimizer whose row variance violates
the optimized lower floor, while its score inequality still passes.
Thus the finite tests distinguish the optimizer-specific and general
parts of the statement rather than conflating them.
