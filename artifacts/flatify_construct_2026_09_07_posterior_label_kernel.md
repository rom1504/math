# Exact posterior-label kernel: a positive special sector and its limits

Set t=2sqrt15, lambda=15/2, R=H2/sqrt2 and

 M(nu,mu)=E exp[-lambda(||a||²+||b||²)+2t a^T R b],

where a has independent scalar law nu and b independent scalar law mu.
This is the root's critical quadratic-dual kernel, with
M(gamma,mu)=1/16 for every scalar probability law mu.

## 1. Exact four-label formula

For four amplitudes a_i in[0,1], let coordinate i independently have law
N(+a_i,1-a_i²)/2+N(-a_i,1-a_i²)/2. Put

 T=[[0,R],[R,0]], A=lambda I-tT, V=diag(1-a_i²), P=I+2VA.

The kernel is exactly

 det(P)^(-1/2) (1/16) sum_{epsilon in{+/-1}^4}
     exp[-(epsilon*a)^T A P^(-1)(epsilon*a)].

This follows by Gaussian integration, including zero residual variances
by continuity. The matrix symmetric to P is positive definite throughout
the parameter cube: 0<=V<=I and I+2A is positive definite.

When a1=a2=a and a3=a4=b this reduces to

 D^-1 exp(-30a²b²/D) cosh(4sqrt30 ab/D), D=16-15a²b².

CAUTION: averaging this two-amplitude formula does not establish the
kernel bound for scalar Gaussian mixtures. The two coordinates of an
endpoint have independent latent amplitudes, not a shared one.

In fact the four-amplitude value at(0,1,0,1) is

 b0=exp(-15/17) cosh(4sqrt30/17)/(2sqrt34)
    =0.0692623945036... >1/15.

Thus the pointwise four-label extension of a <1/15 bound is FALSE.
The diagnostic script `computations/flatify_construct_2026_09_07_four_label_kernel.py`
also numerically maximizes this formula; its optimizer is not a proof of
a global maximum.

## 2. Label averaging nevertheless proves a strict sector bound

Let nu_p=(1-p)gamma+p nu1, 0<=p<=1, where nu1 is the fair sign law.
This is precisely Gaussianization of the channel that independently
reveals a Boolean input with probability p and otherwise reveals nothing.
For ANY p,q, including heterogeneous row parameters,

 M(nu_p,nu_q)<=3/64+b0/4=0.064190598626... <1/15.        (1)

Proof: expand the four independent labels in{0,1}. The corner values are:

 * Zero or one revealed coordinate: a0=1/16.
 * Two revealed coordinates at the same endpoint: a0.
 * One revealed coordinate at each endpoint: b0.
 * Three revealed coordinates:
   c0=[exp(-15/2)+exp(-45/2)cosh(4sqrt30)]/8 <a0.
 * Four revealed coordinates:
   d0=exp(-30)cosh(4sqrt30) <a0.

All expressions follow by integrating the remaining one or two Gaussian
coordinates. There are four exceptional two-reveal corners, of total
probability 4p(1-p)q(1-q)<=1/4. This proves(1).

For an elementary strict certification without trusting decimals, use
sqrt30<11/2, sqrt30>19/4, sqrt34>23/4, exp(1/2)<5/3,
and exp(2)>7 to get b0<38/483<19/240. Consequently the right side of(1)
is strictly below3/64+19/960=1/15. The displayed formulas also directly
give c0,d0<1/16, for example using 4sqrt30<22.

In the typed certificate this supplies the same strict cap gap associated
with the bound M<1/15, for every heterogeneous collection of these scalar
laws. Unlike a pointwise label bound, it uses the probabilities of the
competing label patterns. This is a proved special-sector step, not a
uniform theorem about all recursively generated Boolean profiles.

## 3. Why the same quadratic kernel does not close the recursion

Independent Boolean inputs under a single H2 step have scalar law
nu0=(1/2)delta0+(1/4)(delta_sqrt2+delta_-sqrt2). Its all-zero kernel atom
and its eight positive-maximal one-hot/balanced atoms give

 M(nu0,nu0)>= [1+exp(-45+8sqrt30)]/16
             =0.081662804934... >1/15.

This is an exact obstruction to propagating the critical quadratic-kernel
bound through even one independent Boolean recursion step. It does not
contradict the exact-flat-sector cap theorem: the optimized typed dual can
use nonquadratic potentials, and the quadratic choice discards the fixed
histogram constraint. A general conditional-mean-and-variance envelope
would have to recover that constraint while retaining label correlations.

The bare four-marginal transport Psi also fails a direct Bellman inequality
at t=4, as proved in recursive_rank_two_frames.md Section6. Neither failure
is an impossibility theorem for the actual rank-two construction.

## 4. Bounded attempt to prove a weaker uniform four-label bound

Root observed that the weaker M<=1/12 would still give a useful strict
sector cap gap, and would correctly pass to all independent-coordinate
symmetric posterior Gaussian mixtures. We attempted a rigorous interval
enclosure using the exact polynomial determinant/adjugate in
`computations/flatify_construct_2026_09_07_four_label_interval.py`.
The determinant in u_i=a_i² is

 D=256-120(u0+u1)(u2+u3)+225u0u1u2u3,

and decreases in each u_i on the unit cube. Coordinate permutations and
signs of H2 allow the ordered domain a0>=a1,a2>=a3,a0>=a2. Interval
arithmetic proves the bound on each accepted box; no floating-point
optimizer is used to accept boxes.

The run was deliberately capped at180seconds. It examined52531 boxes,
accepted24440, and had28 unexamined DFS branches remaining. Its status is
INCOMPLETE. It proves no global four-label bound. The natural polynomial
interval evaluation has substantial dependency overestimation, and a
longer unstructured scan was not treated as a research conclusion.

Another tempting shortcut also fails: the kernel is not maximized at
an endpoint of every amplitude coordinate. With the other amplitudes
fixed at(.075,.98,.8), the exact Gaussian-integration formula evaluates
numerically to .0625000969 at first amplitude0, .0620305294 at1, but
.0643854326 at7/8. The gap is much larger than numerical roundoff, though
these decimal evaluations alone are recorded as a diagnostic. Thus a
coordinatewise convexity/endpoint-reduction proof cannot be assumed.
