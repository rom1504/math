# Exact shape of the ternary posterior entropy envelope

Date: 2026-09-06. Independent full-read audit and sharpening of
`continued_director_two_orbit_latent_envelope_2026_09_06.md`, using the
previously audited latent lower theorem. This computes a LOWER-envelope
obstruction class, not an upper bound on the Bellman operator.
The shape proof and the exact rational certificate in Section 5 have
both passed independent director/auditor reconstruction; the auditor
replayed the rational verifier and obtained the same upper endpoint.

## 1. The posterior and two-orbit reductions pass

The posterior coordinates z and b are exact. Reflection preserves its
entropy and squared posterior mean, while enforcing zero signed mean.
The only remaining prior constraint is `E z=p`; every such finite
mixture defines a valid channel by Bayes' rule.

The objective depends only on `E z`, `E z^2 s^2`, and expected posterior
entropy. Four posterior points reproduce those moments. With their
locations fixed, feasible weights form a polytope with two affine
constraints. Its extreme points have at most two supported points.
The Gaussian potential is convex, being a supremum of affine functions
of its variance parameter. Thus at least one of those extreme weight
vectors has objective no smaller. Compactness gives attainment. This
establishes the two-reflection-orbit reduction without a fixed-distortion
or minimax assumption.

Exposing the Gaussian scalar parameter introduces another SUPREMUM,
so interchanging the two optimizations is legitimate. The least concave
majorant of the resulting z function is exactly the mean-z mixture
optimization in the director's equation (7).

## 2. A single possible linear segment

Write

`g_k(z)=h(z)+z max_(0<=s<=1)[h((1+s)/2)+k z s^2]`,

for k>=0. For `2kz<=1` the maximizer is s=0. Otherwise put
`a=atanh(s)>0`; the unique stationary maximum obeys

`a=2k z tanh(a)`.                                          (1)

On this ordered branch the envelope derivative simplifies to

`g_k'(z)=log((1-z)/z)+log(2 cosh a)`.

Differentiating (1), with
`D=1-2kz(1-s^2)>0`, gives

`g_k''(z)=[2kz(1-zs^2)-1]/[z(1-z)D]`.                      (2)

Use `z=a/(2ks)` in the numerator. Its sign is exactly the sign of

`k-K(a)`, where `K(a)=a^2/[2(a coth a-1)]`.                  (3)

The function K increases strictly from `3/2` to infinity. Indeed its
derivative has the sign of

`a sinh(a) cosh(a)-2 sinh^2(a)+a^2`
`=(a/2)sinh(2a)-cosh(2a)+1+a^2`.

The coefficients of a^2 and a^4 vanish, and the coefficient of a^(2j)
for j>=3 is

`2^(2j-1)(j-2)/(2j)! > 0`.

Consequently:

- For `k<=3/2`, g_k is concave on its whole interval [0,1].
- For `k>3/2`, it is strictly concave up to `z_c=1/(2k)`, then
  strictly convex until the unique ordered solution of `K(a)=k`,
  then strictly concave up to one.

The derivative is continuous at z_c. The final inflection occurs
strictly before z=1: at the endpoint, the numerator of (2) is
`-D<0`. Hence no further convexity intervals are hidden.

## 3. Unique common tangent and an explicit scalar equation

For k>3/2 the least concave majorant replaces exactly one interval
`[z_minus,z_plus]`, with `z_minus<z_c<z_inflection<z_plus`, by its
common tangent. Both endpoints are strictly inside (0,1).

One elementary uniqueness proof uses the line slope c. On the left
concave branch and on the right concave branch, maximize `g_k(z)-cz`.
The RIGHT maximum value minus the LEFT maximum value is strictly decreasing in
c, because its derivative is `z_left-z_right<0`. At the two endpoints
of the overlapping derivative ranges the maximizing branch changes.
Thus exactly one slope gives equal maxima. Its line lies above the
entire graph and agrees with it at exactly the two tangent endpoints.

The ordered endpoint has field a>0 from (1). Put

`E=4k sinh(a)/a-2cosh(a)`.

At a stationary point with slope c, `E=exp(c)`. The unordered endpoint
is `z_minus=2/(E+2)`, while

`z_plus=a/(2k tanh(a))`.

The grand-canonical height of the ordered stationary point is
`log(1+2exp(-c)cosh(a))-a^2/(4k)`; that of the unordered point is
`log(1+2exp(-c))`. Equality therefore gives the SINGLE scalar equation

`log[(E+2cosh(a))/(E+2)] = a^2/(4k)`.                        (4)

The desired nonzero root lies on the ordered concave branch, beyond
the solution of `K(a)=k`. It is unique there by the preceding tangent
argument. It lies before the positive solution of `a=2k tanh(a)`
(the endpoint z=1), so E is strictly positive.

The exact concave-hull value at any p is consequently

`cav(g_k)(p)=g_k(p)` outside `[z_minus,z_plus]`,

`cav(g_k)(p)=log(1+2/E)+p log(E)` inside that interval.        (5)

For k<=3/2, simply use `cav(g_k)=g_k`. Thus an arbitrary grid-based
concave-envelope operation is unnecessary. There is one outer k
variable, and at most one uniquely specified common-tangent root.

## 4. Exact outer formula and numerical role

Putting `k=t(1-r)/p` in the director's formula gives

`Envelope(p,t) = -h(p)+t(1-sqrt(p))`
` + sup_(0<k<=t/p) {-p k+(1/4)log[1-(1-pk/t)^2]
                             +cav(g_k)(p)}`.                (6)

For k<=1/2, `g_k(z)=h(z)+z log2`; that part of (6) is exactly the
Gaussian branch. The derivative bound for the outer expression can
also be read directly: an active posterior mixture contributes
`E z^2 s^2` in [0,p], while the logarithmic term has derivative between
zero and `1/(4k)`. Wherever derivatives exist, the outer derivative
therefore lies in `[-p,1/(4k)]`. This supplies a simple Lipschitz
bound on every interval bounded away from k=0, including hull switches.

The companion code uses (1)--(5), not a spectral-profile grid. At
`p=31/32,t=4` it finds, diagnostically,

- Gaussian maximum `-0.04318795110195317` at k approximately .25000786;
- competing ordered/hull local maximum `-0.0468259502773587` at k
  approximately 3.15376757.

These floating extrema are NOT upper certificates for the latent
envelope, and a certified negative latent envelope would still not
upper-bound the Bellman value. The exact shape theorem makes a direct
one-dimensional interval check of this obstruction class feasible.

Files:

- `computations/continued_feedback_ternary_concave_envelope_2026_09_06.py`
- `computations/results/continued_feedback_ternary_concave_envelope_p31over32_t4_2026_09_06.json`

## 5. A fully rational negative certificate

The separate verifier
`computations/continued_feedback_ternary_latent_interval_certificate_2026_09_06.py`
now proves, with exact rational arithmetic,

`Envelope(31/32,4) <= -19678127864847/800000000000000 < 0`.     (7)

Its full output is saved in
`computations/results/continued_feedback_ternary_latent_interval_certificate_p31over32_t4_2026_09_06.json`.
No floating transcendental value, optimizer, or floating comparison enters
the enclosure. Floating values in the JSON are displays only. The verifier
uses a discrete concave hull rather than the analytic tangent formula,
providing a separate enclosure that does not need to locate coexistence.

Here are the complete error bounds and finite operations.

### 5.1 Rational logarithms, entropy, and square roots

For a positive rational x, first write x=2^j y with 1<=y<2. With
`z=(y-1)/(y+1)`, use

`log y = 2 sum_(l=0)^24 z^(2l+1)/(2l+1) + R`,
`0 <= R <= 2z^51/[51(1-z^2)]`.

The same formula at y=2 encloses log 2. Multiplication by the possibly
negative integer j uses the appropriate interval endpoint. Final endpoints
are rounded OUTWARD to multiples of 10^-12 using integer division. Binary
entropy intervals are weighted sums of these logarithm intervals. For a
positive rational x, `isqrt(floor(10^24 x))` provides square-root endpoints
at spacing 10^-12; their squared inequalities are checked exactly.

### 5.2 Inner entropy optimization

For alpha<=1/2 the inner maximum is exactly log 2. For alpha>1/2 the
unique positive maximum s solves

`T(s)=atanh(s)/(2s)=alpha`.

T increases strictly from 1/2 to infinity, as is immediate from its
positive power series. Its values on the grid s=j/100 are enclosed using
the rational logarithm routine. Binary search on the lower and upper
threshold arrays locates rational l<=s<=u; the exact inequalities
`T(l)<=alpha<=T(u)` are verified using outward interval endpoints. The
last endpoint u=1 needs no logarithm because T(1)=infinity.

Write `H(s)=h((1+s)/2)` and m=(l+u)/2. Concavity gives

`H(s)+alpha s^2 <= H(m)+H'(m)(s-m)+alpha s^2`.

The right side is a convex quadratic, whose maximum on [l,u] occurs at
an endpoint. The verifier uses an upper interval endpoint for H(m), and
the lower/upper derivative endpoints for the negative/positive displacements
respectively. The larger of these two rational values bounds the ENTIRE
inner maximum. In the saved run all brackets have width at most 1/100.
Thus a numerical root finder or a smoothness estimate near the transition
is unnecessary.

### 5.3 Concave hull with preserved mean

For each k, compute upper values for g_k on z=j/1000. Their least concave
majorant at p is found by an exact rational monotone-slope stack and linear
interpolation. This equals the maximum average grid-upper value over grid
laws with mean exactly p.

For delta=1/1000, uniformly in s and in |z-z'|<=delta,

`|g_k(z)-g_k(z')| <= omega_k`,
`omega_k=h(delta)+delta log2+2k delta`.                       (8)

Indeed binary entropy has modulus h(delta); the term z H(s) has modulus
delta log2; and k z^2 s^2 has modulus 2k delta. Taking the maximum over s
preserves the bound. Round any z randomly to its adjacent grid endpoints
with probabilities preserving its expectation. Applying this rounding to
any posterior law preserves `E z=p` EXACTLY, and loses at most omega_k in
its objective. Hence

`cav(g_k)(p) <= cav(grid_upper)(p)+omega_k`.

The verifier computes omega_k outward using the same entropy intervals.
This argument covers every posterior law, not only two-point laws or
laws whose support happens to lie on the grid.

### 5.4 The outer parameter and the Gaussian region

Use k=1/2,51/100,...,412/100,128/31, a total of 364 grid points. For
each point the right side of (6) is enclosed by the preceding bounds.
The largest such upper endpoint is

`-23553127864847/800000000000000`, at k=321/100.

The exact outer function is p-Lipschitz on [1/2,t/p]. To see this without
differentiating a possibly switching optimizer, for k2>k1 the concave-hull
increment lies in `[0,p(k2-k1)]`: every feasible posterior law has
`0<=E z^2 s^2<=E z=p`. The logarithmic derivative is in
`[0,1/(4k)]`, so adding the term -pk puts every secant slope between
-p and 1/2. Here p=31/32>=1/2. Every k is within half the maximal grid
gap of a grid point. Adding `p/(2*100)=31/6400` gives (7) for k>=1/2.

For 0<k<=1/2, g_k(z)=h(z)+z log2 exactly. The unrestricted Gaussian
optimization has maximizer

`r=(sqrt(1+16t^2)-1)/(4t)`.

Using the rational square-root and logarithm intervals, its contribution
is at most `-86375902199/2000000000000`, which is smaller than (7).
This bounds even the unrestricted Gaussian branch and therefore also
the portion with k<=1/2.

The two regions exhaust every admissible k. Equation (7) certifies that
NO finite latent channel from this lower-obstruction theorem falsifies
the p=31/32,t=4 candidate. It does NOT bound the Bellman value from above,
does NOT rule out another lower policy, and does NOT prove an upper
construction for the original signing problem.
