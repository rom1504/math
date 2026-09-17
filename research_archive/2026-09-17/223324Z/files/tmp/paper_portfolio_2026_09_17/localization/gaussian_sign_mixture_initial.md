# Frozen mechanism before archive comparison: 2026-09-17

Actual full signing A, bounded spectral norm ||A||<=L sqrt(n).
For 0<t<=1/L, K_s=I+s t A/sqrt(n) are Gaussian correlation matrices.
Coordinate sign rounding gives EXACT covariance

    E h h^T=I+s rho_n A, rho_n=(2/pi)arcsin(t/sqrt(n)).

The equal mixture of the two laws is EXACTLY isotropic. It can nevertheless
have smaller absolute response on both absolute near-level polarities:
if a uniform Gaussian sign-sum CLT is valid, then

    E_mixture |h dot x|/sqrt(n)
      =kappa/2 [sqrt(1+4t H_A(x)/(pi n^(3/2)))
                 +sqrt(1-4t H_A(x)/(pi n^(3/2)))]+o(1).

The right side is strictly smaller than kappa away from energy zero.
This uses latent variance heterogeneity, not covariance reduction after
mixing. It is a bridge-column law, not same-map Krivine rounding.

At c=1/2,L=t=1 the scalar slope is
sqrt((1+sqrt(1-4/pi^2))/pi), slightly above 3/4. Thus a sharp negative
threshold may accompany the positive uniform-response theorem.

Main obligations: prove the CLT uniformly over x without a false weak-
dependence assumption; show adequate column concentration or use exact
bounded Gaussian density away from spectral singularity; quantify parent
value with all old energy levels, all new spins, and child offsets retained.
No optimizer spectral bound is assumed or claimed.
