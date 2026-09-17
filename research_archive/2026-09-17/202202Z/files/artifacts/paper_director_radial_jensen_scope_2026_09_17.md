# What covariance-only radial certificates can and cannot certify

2026-09-17. Director derivation. **Independently reconstructed by the
discrepancy researcher.** This is a limitation of a particular upper certificate,
not a lower bound on actual sign-column response. It is motivated by the
new exact nonlocal covariance laws and finite correlation-polytope tests.
No external novelty claim is made.

Let A be a hollow full signing of order n, d=binom(n,2), and
Q(A)=c_n n^(3/2). Suppose a family of genuine sign-vector laws nu_s has

    E_(nu_s) h h^T = I+s A/sqrt(n).

Means need not vanish for the statements below. Symmetrizing h by an
independent global sign makes them vanish without changing anything used.

## 1. Actual energy constrains every radial component

For EVERY realizable s,

    E_(nu_s) H_A(h)=s d/sqrt(n),
    |s| <= Q(A)sqrt(n)/d = 2c_n n/(n-1).             (1)

This uses complete sign support: every squared edge coefficient is one.
It is distinct from positive-semidefinite feasibility. The sharper one-sided
inequalities replace Q by max H_A or max(-H_A), respectively.

Consider any mixture of these laws with mixing variable S and ES=0. It is
exactly isotropic. For a fixed spin x with H_A(x)=e n^(3/2), the
componentwise Cauchy--Schwarz/Jensen certificate is

    E|h dot x|/sqrt(n) <= E sqrt(1+2eS).              (2)

Equation (2) is an UPPER certificate, potentially very loose.

## 2. A sharp floor on this certificate, not on the response

Let c_n -> c, where 0<c<1/2, and take an absolute ground word, so
|e|=c_n. By (1), S is in [-a_n,a_n], a_n=2c_n n/(n-1).
For all sufficiently large n this entire interval lies in the domain of
sqrt(1+2e s). The concave square-root graph lies ABOVE its endpoint chord.
Taking expectations and ES=0 proves

    E sqrt(1+2eS)
      >= .5[sqrt(1-2c_n a_n)+sqrt(1+2c_n a_n)]
      -> J(c):=.5[sqrt(1-4c^2)+sqrt(1+4c^2)].       (3)

Allowing every possible radial endpoint and every mean-zero mixing law
does not improve this bound. Equality in this relaxed certificate problem
uses the two extreme endpoints equally; those endpoint laws need not be
physically realizable. Thus using only their second moments and applying
(2) cannot certify a coefficient below J(c).

Solving J(c)<3c/2 on 0<c<=1/2 gives exactly

    c > 6/sqrt(145) = .49827... .                    (4)

Indeed J(c)^2=(1+sqrt(1-16c^4))/2. The desired inequality requires
4.5c^2-1>0; squaring then gives c^2>36/145, which also implies the
required sign condition. At equality the two coefficients agree.

Consequently the already reported strict upper bound below .494 is in
the regime where covariance-only radial Jensen certificates are too
expensive for infinitesimal near-optimal extension. This does NOT exclude
these actual column laws: their true absolute means may be substantially
smaller than their standard deviations. The finite exact n=8 witnesses
with identical covariances but different mean responses illustrate exactly
that missing information. Gaussian absolute-response comparison or a
different nonlocal calculation could recover it.

## 3. Contrast with Gaussian response and with Gaussian-sign rounding

If, as an EXTRA theorem, each component had Gaussian absolute mean at
its own variance, the coefficient would instead be kappa J(c),
kappa=sqrt(2/pi). Then kappa J(c)<3c/2 precisely when

    c^2 > 72pi/(256+81pi^2).                         (5)

The threshold is about .463, below the current reported upper bound but
above the current reported lower bound. The algebra uses
sqrt(1-16c^4)<(9pi/4)c^2-1 and checks the right-hand side is positive.
Neither endpoint feasibility nor that Gaussian-response hypothesis is
proved for general near-minimizing A. Covariance matching alone supplies
neither. This calculation is a transfer diagnostic, not a convergence
criterion or an improved signing bound.

Gaussian-sign rounding supplies a different, provable response law, but
its radial covariance loses 2/pi and has the independently reconstructed
coefficient barrier .750855... . Nonlocal exact sign laws can evade that
particular barrier. Equations (3)--(5) identify what their subsequent
response calculation must retain rather than prematurely replacing it
by a second-moment estimate.
