# Width-pressure recurrence: audit, temperature iteration, and exact scope

Date: 2026-09-07. The recurrence in
`decisive_bridge_width_pressure_cross_order_2026_09_07.md` passes independent
reconstruction. This note derives its positive iterated consequences and
gives an explicit non-realizability-free countermodel to deducing convergence
from this recurrence and temperature regularity alone. The countermodel is
NOT claimed to be the pressure of actual signings.

## 1. Recurrence audit and legitimate temperature facts

Put f_N(beta)=Psi_N(beta)/N-log2, where Psi is the minimum over actual
signings of one half the sum of the two branch log partition functions.
The exact finite bounds are

    0<=f_N(beta)<=min{(1-1/N)beta^2/4, beta sqrt(log2)}.

Every minimizing parent has Q(beta A/sqrt(N))<=C_beta N, with

    C_beta=min{2log2+beta^2/2, 2log2+2beta sqrt(log2)}.

The first bound is annealing followed by Jensen, not choosing A after a
Gaussian sample. For the second, independent sign rounding and the standard
maximum Gaussian-free subgaussian bound give a signing with
Q(A)<=sqrt(log2)N^(3/2). Either bound therefore applies to the deterministic
minimizer selected at this beta.

For a prescribed cut, switching one whole block represents every partially
deleted parent as a convex combination of two switch-equivalent parents.
Thus its actual cap never exceeds the original cap. The audited quadratic
fluctuation bound holds along this ENTIRE fixed-parent path; its derivative
at zero vanishes by block switching. Integrating it before taking any new
minimum is legitimate. Each fixed energy branch factors at the zero-bridge
endpoint, and the average of its two logarithms is exactly additive.
Consequently, if r=m/N and s=n/N, m+n=N,

    f_N(beta)>=r f_m(beta sqrt(r))+s f_n(beta sqrt(s))
                  +(beta^2/2)rs c(C_beta,1/min(r,s)),       (1)

with the explicit c in the canonical theorem. Factors beta^2 mn/N and
one half from integrating the second derivative both check.

Each fixed-signing f is even, convex, nondecreasing on beta>=0, and has
f(0)=0. After minimizing over signings, nondecreasingness of f_N and of
f_N(beta)/beta survives. Convexity of the MINIMUM itself does not follow
from convexity of the branches and is not assumed below. Its ground slope
is the normalized optimized HALF-WIDTH, not the original absolute cap.

## 2. Genuine iteration with rescaled temperatures

For N=2^k ell, repeated balanced use of (1) gives the exact finite bound

    f_N(beta)>=f_ell(beta 2^(-k/2))
       +(beta^2/8) sum_{d=0}^{k-1} 2^(-d)
                              c(C_(beta 2^(-d/2)),2).      (2)

The same finite number of balanced levels with rounded block sizes applies
to all sufficiently large N. Taking N to infinity first at fixed k,
discarding the nonnegative leaf terms, then increasing k proves

    liminf_N f_N(beta)>=S(beta)
      :=(1/8) sum_{d>=0} u_d^2 c(C_(u_d),2),
      u_d=beta 2^(-d/2).                                  (3)

This is a positive all-order pressure floor at every beta>0. It is not a
ground-value recurrence. Indeed c(C_u,2)<=exp(-16C_u), since K_G>=1.
As log2>=1/2, we have C_u>=1 for u<=1, and C_u>=1+u/2 for u>=1.
Hence c<=e^(-16) in the first region and c<=e^(-16-8u) in the second.
Using e^(-8u)<=u^(-4) for u>=1 and summing the two geometric tails gives

    S(beta)<=e^(-16)/2 for every beta>0.                    (4)

Thus the explicitly accumulated credit itself is bounded as beta grows;
its contribution divided by beta vanishes. Equation (4) is an upper bound
on this particular proven floor, not an upper bound on the actual pressure.

## 3. Homogeneous temperature integrals

For 1<a<2 the Mellin functional

    I_N(a)=integral_0^infinity f_N(beta) beta^(-a-1) d beta

is finite, uniformly in N. The endpoint range is dictated by the quadratic
small-temperature bound and linear large-temperature growth. Integrating
(1), with the exact substitution u=beta sqrt(r), gives

    I_N(a)>=r^(1+a/2) I_m(a)+s^(1+a/2) I_n(a)+K_a(r),     (5)

    K_a(r)=(rs/2) integral_0^infinity
                  c(C_beta,1/min(r,s)) beta^(1-a) d beta>0.

The credit is finite by the bounds used in Section 2. Thus
N^(1+a/2) I_N(a) is superadditive, with a positive extra credit. Its
natural scale, however, has exponent 1+a/2>3/2, not the linear scale
required by ordinary Fekete. Balanced iteration discounts the inherited
I_n by 2^(-a/2) each time; it yields a universal lower floor, not memory
of a liminf optimizer.

An Abelian ground-slope transform makes the same issue explicit:
(log R)^(-1) integral_1^R f_N(beta) beta^(-2) d beta tends to the ground
slope as R grows. The child dilation contributes sqrt(r), so (1) then
retains weights r^(3/2), not r. The integrated positive credit is divided
by log R and vanishes. No temperature dilation is erased for free.

## 4. A rigorous oscillating countermodel to a stronger inference

The following profiles satisfy MORE temperature regularity than is known
for the optimized pressure, together with (1) for its stated c. Nevertheless
their ground slopes do not converge.

For integers n>=1 put

    a_n=9/20+(1/100) sin(log log(n+e)),
    v_n=1-1/n,
    F_a(beta)=a(sqrt(4a^2+beta^2)-2a),
    psi_n(beta)=log2+v_n F_(a_n)(beta).

Here a_min=11/25=.44 and a_max=23/50=.46. The functions are even and
convex, nondecreasing for beta>=0, and

    F_a(beta)=beta^2/4+O(beta^4),
    0<=F_a(beta)<=min(beta^2/4,a beta),
    beta v_n a_n<=psi_n(beta)<=beta v_n a_n+log2.

The lower entropy bound uses 2a_max^2<log2. In particular the EXACT
universal small-beta coefficient (1-1/n)/4 is respected, as are both
displayed annealed/cap upper bounds. But psi_n(beta)/beta tends to
v_n a_n as beta grows, and these slopes have liminf .44 and limsup .46.

Here is the all-order recurrence check. Differentiating in log n gives
|a_N-a_m|<=(1/100)log(N/m). Also F_a increases in a and
partial_a F_a<=F_a/a. Convexity in beta gives F_a(beta sqrt(r))<=sqrt(r)F_a(beta).
With a=a_N and v_m,v_n<=v_N, the weighted child contribution is at most

    v_N F_a(beta)[r^(3/2)+s^(3/2)
          +(1/44){r^(3/2)log(1/r)+s^(3/2)log(1/s)}].

Let Delta=1-r^(3/2)-s^(3/2). The elementary inequalities e^x-1>=x
and sqrt(r)<=(1+r)/2 imply respectively

    r^(3/2)log(1/r)+s^(3/2)log(1/s)<=2 Delta,
    Delta>=rs.

The recurrence gap without its credit is therefore at least
v_N F_a(beta)(21/22)rs. Since v_N>=1/2 for N>=2 and

    F_a(beta)=a beta^2/[sqrt(4a^2+beta^2)+2a]
             >=(2/5) beta^2/(beta+2),

this gap is at least

    (9/50) beta^2 rs/(beta+2).                             (6)

On the other hand eta=1/min(r,s)>=2 and

    c(C_beta,eta)<=e^(-16)             for beta<=1,
    c(C_beta,eta)<=e^(-16-8beta)       for beta>=1.

Both bounds are less than (9/25)/(beta+2): for beta<=1 use e>2,
and for beta>=1 note that (beta+2)e^(-8beta) decreases from 3e^(-8).
Thus (6) dominates the positive credit in (1), for EVERY m,n,beta.

The separate hot-child inequality from
`decisive_bridge_vertex_swap_overlap_convolution_2026_09_07.md`, Section 4,
does not rule out this countermodel either. That exact argument chooses
one balanced gauge after averaging BOTH branch log overlap likelihoods;
Jensen bounds each average by -log binom(N,N/2). It gives

    f_N(beta)>=.5 f_(N/2)(sqrt(2) beta)
                   -[Nlog2-log binom(N,N/2)]/(2N).        (7)

The countermodel satisfies (7) even after OMITTING its nonnegative error:

    .5 F_a(sqrt(2) beta)=F_(a/sqrt(2))(beta),
    a_N>=.44>.46/sqrt(2)>=a_(N/2)/sqrt(2),
    v_N>=v_(N/2).

Thus both cold and hot temperature directions, with their respective
credits, remain compatible with nonconvergent ground slopes.

This is not an example of actual signing pressures or a disproof of
convergence. It proves that (1), even combined with full temperature
convexity, the exact quadratic coefficient, entropy/cap bounds and
monotonicity, cannot alone imply a summable recurrence forcing convergence.
Any successful transform must use an additional actual-signing property,
such as a sharper relation between cross-direction curvature and the
internal temperature response. Merely integrating the existing bounds is
insufficient.

## 5. Finite diagnostic

`computations/decisive_audit_width_pressure_countermodel_2026_09_07.py`
checks the analytic inequalities on all splits through order128 and a
logarithmic temperature grid. The proof of the countermodel is Section 4;
the finite computation checks implementation and normalization only.
It passed all 658368 split/temperature cases. The bridge agent independently
reconstructed the amplitude derivative, order-drift absorption, entropy
band, and comparison against the explicit credit, and reported PASS.
