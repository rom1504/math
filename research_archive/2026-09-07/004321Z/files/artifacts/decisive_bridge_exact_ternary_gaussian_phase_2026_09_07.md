# Exact Gaussian phase for E at p=31/32, t=4

Date: 2026-09-07. Status: exact analytic reduction plus directed-interval
rectangle exclusion. The director supplied the three-atom reduction idea;
the bridge agent independently proved its monotonicity and ran the exclusion.
This is a phase theorem, not an approximate optimization of E.

Let nu_p=(1-p)delta_0+(p/2)(delta_(-1/sqrt(p))+delta_(1/sqrt(p))). Then

`E_4(nu_(31/32)) = g_4(1)`.

The trivial, uninformative channel attains this value. All genuinely
informative channels have strictly smaller value. Consequently a proposed
Paley/nonlinear-mean mechanism cannot improve on E at this particular
parameter: the universal Gaussian/mask-blind sector already attains E.
This does not assert that the Gaussian phase holds for all p or t.

## 1. Exact three-point reproduction support theorem

For fixed lambda>0, rate-distortion duality gives

`J_lambda(nu)= -max_pi sum_x nu(x) log integral exp[-lambda(x-y)²] dpi(y)`.

One may restrict reproduction y to the convex hull of the input alphabet
(posterior means minimize squared error), so the maximizing probability
measure exists. Symmetrizing it preserves feasibility and improves the
concave objective. For the symmetric ternary source, its KKT witness has form

`q(y)=exp(-lambda y²)[a+b cosh(2 lambda y/sqrt(p))] <=1`,

with a,b>0, and equality on the reproduction support. Put r=a/b and z>0.
The positive critical points of q correspond to horizontal-level crossings
of `h(z)=sinh(z)/[z(r+cosh(z))]`. The sign of h' is the sign of r-R(z), where

`R(z)=[sinh(z)cosh(z)-z]/[z cosh(z)-sinh(z)]`.

The numerator and denominator have positive power series beginning at z³,
with respective coefficients `4^j/(2j+1)!` and `2j/(2j+1)!`, j>=1. Their
coefficient ratios 4^j/(2j) are strictly increasing. Expanding the derivative
of the quotient pairwise, or using weighted covariance of j with this ratio,
proves R strictly increasing, from 2 to infinity. Thus h is either decreasing
or increases once and then decreases once. There is at most one positive
LOCAL MAXIMUM of q. Since its support lies at global maxima q=1, the symmetric
optimal reproduction measure is supported on `{0,+a0,-a0}` for some a0.
Boundary y>1/sqrt(p) cannot maximize q since every Gaussian summand decreases
there. This proves the support reduction without assuming it numerically.

## 2. Eliminate the reproduction mixture weight exactly

Write pi=w delta_0+(1-w)(delta_a0+delta_-a0)/2 and put z=a0 sqrt(p), so
0<=z<=1. Define

`A=exp(lambda z²/p)-1`, `B=cosh(2 lambda z/p)-1`,
`v=(1-w)exp(-lambda a0²)/[w+(1-w)exp(-lambda a0²)]`.

Then

`J_lambda=lambda-max_(z,v) G(A,B,v)`,
`G=p log(1+vB)-log(1+vA)`, 0<=v<=1.

The maximizing v is the clipping to [0,1] of
`(pB-A)/[(1-p)AB]`, with the zero cases interpreted continuously. In the
interior, the maximum simplifies to

`p log p+(1-p)log(1-p)-p log A+log B-(1-p)log(B-A)`.

The two endpoint values are 0 and `p log(1+B)-log(1+A)`. Therefore

`E_t(nu_p)=max_(0<lambda<=t, 0<=z<=1)
                [c_t(lambda)-lambda+max_v G(A,B,v)]`.

## 3. Low precision: exact subgaussian argument

For p>=1/3 the normalized ternary source is unit subgaussian. Indeed
`E X^(2j)=p^(1-j)<=(2j-1)!!` for all j>=1: the induction step multiplies
the left side by 1/p<=3 and the right side by 2j+1>=3. Its exponential
moment series is therefore coefficientwise dominated by exp(s²/2).
Entropy duality implies `I(X;L)>=E[E(X|L)²]/2`.

Thus J_lambda=lambda for lambda<=1/2. The maximum of c_4(lambda)-lambda
in this interval is exactly g_4(1), attained at

`lambda_star=4(1-rho)`, `rho=(sqrt(257)-1)/16`,

which lies strictly below1/2. The exact rational log/square-root routines
already used by the campaign certify

`g_4(1) >= -622136276211/800000000000 > -777671/1000000`.

## 4. High precision: a finite covering certificate

The companion script covers the CLOSED rectangle
`lambda in [1/2,4], z in [0,1]` with dyadic rectangles. On a rectangle
`[l,u] x [a,b]`, G is increasing in B and decreasing in A, so replace

`A by exp(l a²/p)-1` and `B by cosh(2u b/p)-1`.

Also c_4(lambda)-lambda is decreasing for lambda>=1/2: c_4'(1/2)=7/15<1
and c_4 is strictly concave. Thus the interval upper bound is

`c_4(l)-l + max_v G(A_min,B_max,v)`.

An adaptive subdivision covers the entire domain with 212505 accepted
rectangles (425009 tree nodes). Every accepted upper bound is STRICTLY
less than `-777671/1000000`, itself strictly below the Gaussian lower bound.
The largest floating diagnostic upper was -0.7776710006464511; that number
is not used as an exact constant.

Replay:

`.venv/bin/python computations/decisive_bridge_gaussian_phase_rectangles_2026_09_07.py --verify`

The floating pass only chooses subdivisions. Every accepted rectangle is
re-evaluated using mpmath1.3.0 directed intervals at60 decimal digits. The
installed interval exp/log implementations call their underlying arbitrary-
precision functions separately with round_floor and round_ceiling (inspected
directly). The interval endpoints are dyadic; the cases selecting the
optimized-v branch are checked with Python EXACT rational arithmetic.
All rectangle endpoints are dyadic, so the split tree covers the domain
exactly. An incorrect floating prune therefore fails its interval assertion
instead of silently producing a certificate. The Gaussian comparison is
separately replayed by the campaign's exact rational routines.

## 5. Phase consequence and uniqueness

Combining the two precision regions proves the displayed equality. In the
low-precision proof, lambda_star<1/2 forces zero posterior-mean variance at
equality, and then zero mutual information. High precisions have a strictly
positive certified gap. Hence the sole optimal posterior law is the source
itself: redundant labels independent of X are harmless, informative labels
are not optimal. This is a local Gaussian-phase statement for the recursive
Bellman certificate after H=E, not a universal identification of actual
minimax values or arbitrary ensemble pressures.
