# Independent audit: quenched Gaussian-sign universality

2026-09-07. PASS of the complete canonical proof in `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`. This is a genuine quenched comparison, unlike the earlier deterministic-source MGF argument. It does not solve the resulting Gaussian parent optimization.

## Conditional replacement and smoothness

Conditional on W in G=W+sqrt(epsilon)z, the sign coordinates are independent. Their replacements F_i=f(W_i)+sigma(W_i)xi_i match both conditional moments exactly, and every third coordinate derivative of the Gibbs log partition is bounded by 8|lambda|³. Centered conditional third moments are uniformly bounded, including conditional means near +/-1. Therefore the Lindeberg cost C|lambda|³d is valid for arbitrary deterministic offsets.

The replacements have exactly the ACTUAL sign covariance, not merely its first chaos: distinct-coordinate covariance is E f(W_i)f(W_j), and the diagonal is one.

Global sigma-derivative bounds hold. An elementary alternative to the source's Mills estimates is to integrate the normal density over [z,z+1/(1+z)] to get Q(z)>=c phi(z)/(1+z), while Q(z)<=exp(-z²/2)/2 follows by translating the positive half-line density. Since sigma(z)=2sqrt(Phi(z)Q(z)), its logarithmic derivative is O(1+z), its logarithmic second derivative is O(1+z²), and its Gaussian decay dominates these factors. Thus after rescaling, |sigma'|<=C epsilon^(-1/2) and |sigma''|<=C epsilon^(-1). No hidden bounded-auxiliary-noise assumption occurs.

## Exact Stein algebra and the adaptive Gibbs issue

The source's two-dimensional block covariance T=diag(R-epsilon I,I) has bounded operator norm. Its local resolvent gradient a_i and ordinary gradient b_i depend ONLY on block i. The kernel Gamma_ij=a_i^T T_ij b_j has the correct orientation for integration by parts in F_i, and E Gamma_ij is the exact sign covariance.

Differentiating Gamma is supported only on blocks i,j. Applying the Gaussian OU covariance identity to Gamma_ij and the pressure Hessian gives two stars, with covariance factors T_ij T_ik or T_ij T_jk. Singular covariance in the W component is harmless via the Gaussian square-root representation.

I also reconstructed the same expansion in independent latent coordinates. Writing a_i=f'(W_i)+sigma'(W_i)xi_i and c_i=sigma(W_i), the scalar-covariance version has exactly five types of terms:

    S_ij S_ik, S_ij S_jk, delta_ij S_ik,
    S_ij delta_ik, S_ij delta_jk,

where S=R-epsilon I. Their scalar factors remain separated by i,j,k even when the third factor is evaluated at a DIFFERENT OU environment. This independently confirms the block-index proof and its treatment of diagonal auxiliary terms.

The third Gibbs cumulant expands into five expectation products with total absolute coefficient sum six. Each product is represented by at most three replicas and has coordinate form v_i w_j z_k, with all coordinates bounded by one. For fixed environments AND fixed replicas, every star obeys

    |sum_i d_i v_i (S_1 diag(e)w)_i (S_2 diag(f)z)_i|
       <=d ||S_1||op||S_2||op ||d||infinity||e||infinity||f||infinity.

This bound is pointwise, so adaptive dependence of the replicas on the environments does not invalidate it. This is precisely the missing ingredient that could not be inferred from a deterministic-source MGF alone.

There are only finitely many block-coordinate choices. The coefficient product is bounded by a constant times the cube of one plus the two auxiliary Gaussian maxima; its expectation is O(log(d+1)^(3/2)), uniformly in the OU interpolation time. Integrating both OU time and Gaussian comparison time is finite; the only comparison-time factor is sqrt(t), not an endpoint singularity.

Consequently the canonical error C_(epsilon,K)|lambda|³d log(d+1)^(3/2) checks in full. Differentiation is justified by the displayed Gaussian-polynomial moment bounds. No pressure Hessian is assumed small.

## Quantitative constants and cap transfer

The smooth local derivative products above give a permissible explicit dependence C K² epsilon^(-2) for 0<epsilon<=1, K>=1. The worst star contains one second derivative O(epsilon^(-1)) and two first derivatives O(epsilon^(-1/2)) each; terms involving auxiliary-coordinate derivatives are smaller. The conditional replacement has an absolute constant.

For d=n², lambda=beta/sqrt(n), and exp(O(n)) configurations including an absolute-polarity bit, the normalized expected-cap error is

    O(1/beta)+O_(epsilon,K)(beta² n^(-1/2)log(n+1)^(3/2)).

The canonical beta=n^(1/12) is valid. Optimizing the stated bound at fixed epsilon,K permits beta=n^(1/6)/sqrt(log(n+1)), giving O_(epsilon,K)(n^(-1/6)sqrt(log(n+1))). Arbitrary actual-child deterministic energy offsets are unchanged throughout.

For R=I-rho A tensor D/p with p=||A||op||D||op and fixed rho<1, this comparison applies to EVERY pair of actual sign children. Its Gaussian comparator has exactly covariance I-[(2/pi)arcsin(rho/p)]A tensor D. No spectral-flatness assumption is needed. A favorable value for this Gaussian parent is still unproved; expected-cap equivalence alone is not the requested recurrence.
