# Independent director audit: high-value zero-strip counterexample

Date: 2026-09-06. This audit reconstructs the mathematical steps rather than
using the response researcher's positive conclusion as a premise.

The audited claim is existence of a genuine marked-creation mask H with
J(H)>.4301875 and H=1 on |UH|<=.001. It does not improve the banked lower
bound and does not concern actual optimizing signings or the original limit.

## Analytic reconstruction

1. The inverse image g of the banked unit Gaussian V has the exact decomposition
   g=P+A R_a chi_D, where the innovation resolvent multiplies innovation degree
   l by 1/(a+l). All removed anchor features have innovation degree zero;
   hence their correction coefficient is rho_T-(A/a)c_T. The code's grouped
   squared coefficients omit their constant and degree-200 terms correctly.
2. The full Ornstein--Uhlenbeck resolvent is positive. Adding and subtracting
   its untruncated value bounds ||g_-||_2 by the negative part of P and the
   L2 truncation remainder; positivity of the finite polynomial is not used.
3. For P=delta_0+R with ER=0 and delta_0>0, the pointwise bound
   (-delta_0-R)_+^2<=R^4/(16 delta_0^2) follows by optimizing one scalar ratio.
   The normalized Hermite multiplication formula in the executable is correct;
   the 21 features use only the first nine anchor variables, so its nine-index
   coefficient tuples omit no variables. Parseval computes E R^4 from R^2.
4. At each total degree d the resolvent tail has binomial innovation-degree
   weights. Coupling Bin(d+1,s^2)=Bin(d,s^2)+Bernoulli(s^2) proves that its
   inverse-square multiplier decreases in d. Thus the first omitted even
   degree 202 gives a valid bound on the entire remaining tail.
5. The affine ball intersected with the orthogonality plane and covariance
   halfspace is closed, convex, and complete. Since H0 is binary and the
   feedback K is supported on 1-H0, its image increment UK is orthogonal to
   W0. The one-sided bound <g,K>>=-||g_-||_2||K||_2 gives the required invariant
   halfspace once the norm bound closes. No arbitrary covariance realization
   is substituted for the actual creation map.
6. For u=(alpha-cw/v)/sqrt(1-c^2/v), direct differentiation gives the signs
   asserted in the proof: its c derivative has sign alpha*c-w and its v
   derivative has sign w(2v-c^2)-alpha*v*c. The analogous plus threshold has
   the required signs without the latter restriction. The checked rational
   inequalities apply throughout the stated parameter box.
7. The contraction bound includes every Gaussian direction. Decomposing a
   direction along W, its V residual, and their independent complement gives
   a positive 2-by-2 weighted second-moment matrix and a separate scalar term.
   The trace dominates the largest eigenvalue, and t2>=t0 dominates the
   independent term. The cross term is not assumed zero.
8. The worst-parameter tail and truncated second moment increase with w on
   the interval because u>=sqrt(2); the Gaussian density and triangular q
   decrease. Endpoint products are therefore genuine Darboux upper sums.
   Segment integration on the convex domain turns the all-direction estimate
   into a contraction. The piecewise linear corners have Gaussian measure zero.
9. Finally B(v,c)=E|W|1{|V|>alpha} increases with v at fixed c, and
   0<=partial_c B<=2phi(alpha). If necessary decrease c first to c0-er (which
   is below sqrt(p0)) before comparing v with p0; this avoids evaluating the
   formula outside its covariance domain. The resulting value loss is at most
   2phi(alpha)er plus the separately bounded removed absolute energy.

The scalar covariance can alternatively be kept uniformly nonsingular using
Delta perpendicular to W0 and its norm bound; singular conditional laws, if
included by an enlarged envelope, are handled by continuous limits.

## Verdict and scope

The analytic proof passes this independent reconstruction. The director's
exact replay passed; the full report is saved at
`computations/results/resumed_response_strip_certificate_2026_09_06.json`.
This invalidates high-value zero-strip slack at the stated threshold, but not
slack at true maximizers of the marked variational functional. It motivates
testing the distinct feasible center-best-response mechanism; the stronger
nonlinear transport theorem needed by that mechanism is still open.
