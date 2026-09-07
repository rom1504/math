# Exact PSD-plus-flip-slack interpolation for optimized width pressure

This is an exact finite-dimensional derivative formula. Its remaining
weighted slack is not proved to have the sign or summable error needed
for convergence.

## 1. Weighted width minimizers and their exact edge slack

Let F(A,lambda)=.5(log Z_A^+(lambda)+log Z_A^-(lambda)), where the
positive edge magnitudes lambda_e may vary. At an actual globally
minimizing signing A define

    C^+_ij=<xi xj>_+, C^-_ij=<xi xj>_-,
    P=C^+ circ C^-,
    d_e=.5 A_e(C^+_e-C^-_e)=partial_(lambda_e)F,
    s_e=tanh(2lambda_e)(1-P_e)-2d_e,
    a_e=s_e/tanh(2lambda_e).

Both C matrices, and hence their Schur product P, are PSD correlation
matrices. The paired edge-flip calculation proves s_e>=0 and a_e>=0
for lambda_e>0.

More explicitly, delete e and let c^+_0,c^-_0 be its two cavity spin
correlations. With t=tanh(lambda_e), edge insertion multiplies the
paired partition by

    cosh(lambda_e)^2
       [1+A_e t(c^+_0-c^-_0)-t^2 c^+_0 c^-_0].

Actual minimizing sign optimality is EXACTLY

    A_e(c^+_0-c^-_0)=-b_e, b_e=|c^+_0-c^-_0|.

Writing q_e=1-t b_e-t^2 c^+_0 c^-_0>0, direct rational algebra gives

    a_e= b_e(1-t^2)^2/(2t q_e).                     (1)

Equivalently, if Delta_e=F(A^e,lambda)-F(A,lambda)>=0 is the actual
width log-pressure cost of flipping e, then

    a_e=[exp(2Delta_e)-1]/sinh(2lambda_e)^2.         (2)

Thus the nonnegative slack is an exact normalized edge-flip cost, not
an unverified remainder in a conditional-variance approximation.

## 2. A log-cosh path with an exactly signed PSD term

Fix a split N=m+n, r=m/N, and put

    z_i=sqrt(n/m) on the first block,
    z_i=-sqrt(m/n) on the second block,
    D=z z^T.

D has zero row sums. Its off-diagonal entries are n/m and m/n
within the respective blocks, and -1 across them. Set

    tau_0=log cosh(2beta/sqrt(N)),
    log cosh(2lambda_e(u))=tau_0(1+uD_e), 0<=u<=1.

For 0<u<1, lambda'_e=tau_0 D_e/[2tanh(2lambda_e)].
At almost every u, use any active signing of the finite optimized
envelope F_*(u)=min_A F(A,lambda(u)). Then

    F_*'(u)
      =-tau_0 z^T P z/8 -tau_0 sum_e D_e a_e/4.       (3)

Indeed d_e=.5 tanh(2lambda_e)(1-P_e-a_e), and

    sum_(i<j) D_ij(1-P_ij)=-.5 z^T Pz.

The diagonal terms vanish because P_ii=1, and the all-one term
vanishes because D has zero row sums. The first term of (3) is
nonpositive by the Schur product theorem. There is NO Taylor remainder.

The second term still has both signs: D_e is positive within blocks
and negative across them. Nonnegativity of a_e by itself does not
control their difference. Formula (1) does not remove it; it expresses
it as the absolute difference of the two edge-deleted branch responses.

## 3. Endpoint normalization and the alternative squared-weight path

The initial point is exactly homogeneous beta/sqrt(N). At u=1 the
cross edges vanish, and the within-block magnitudes satisfy

    log cosh(2lambda_i)=tau_0/r_i.

For comparable block sizes, their normalized temperatures are
sqrt(Nr_i) lambda_i=beta+O_(beta,r)(1/N). Hence the optimized endpoint
is Psi_m(beta)+Psi_n(beta)+O_(beta,r)(1). The O(1) total-pressure
error follows, for example, from the exact homogeneous edge-optimal
bound 0<=partial_beta Psi_d(beta)<=beta(d-1)/2 on compact positive
temperature intervals. The normalized error is O(1/N).

The envelope is continuous at the zero-cross endpoint. Formula (3)
is interpreted a.e. before that endpoint; it can be integrated as
the derivative of this finite optimized pressure. No division at a
zero edge is used to assert an optimizer condition.

For comparison, along the simpler squared-weight path
lambda_e=beta sqrt(1+uD_e)/sqrt(N), the same computation yields

    F_*'(u)=-beta^2 z^T Pz/(4N)
                -.5 sum_e lambda'_e s_e+O_(beta,r)(1). (4)

The error is uniformly integrable through the zero-cross endpoint:
|tanh(2lambda)-2lambda|<=8lambda^3/3, and
|lambda'_e|lambda_e^3=beta^4 |D_e|(1+uD_e)/(2N^2).

## 4. Precisely what is still missing

To obtain a same-temperature almost-superadditive WIDTH comparison
from (3), a sufficient integrated estimate is

    integral_0^1 sum_e D_e a_e(u) du >= -O(N^(1+alpha)),
                         for some alpha<1,

because tau_0=O(1/N). A stronger estimate retaining the favorable
PSD term would also suffice. An unspecified o(N^2) slack budget
alone gives only an o(N) pressure error; a summable modulus is needed
for the usual almost-additive convergence argument.

At u=0 an average over balanced partitions controls the slack at
lower order. At u>0, however, the optimizer already belongs to the
fixed anisotropic profile. Relabeling the whole optimized system
preserves the envelope and does not freely choose a better partition
for its current slack matrix. The earlier adaptive-partition
obstruction therefore still applies.

Finally, even convergence of this width objective is not automatically
convergence of the original absolute-cap minimum. Its common-polarity
endpoint must be handled separately in any original-limit theorem.

## 5. Endpoint reduction to an actual flat-bridge variance problem

For equal children of order m=N/2, fix an active endpoint child pair
and their branch covariance matrices C_1^s,C_2^s. For a full-sign
bridge B define

    V(B)=sum_(s=+,-) Tr[B^T C_1^s B C_2^s].

The child spin means vanish exactly by spin reversal. Thus at cross
amplitude epsilon the width bridge gain is

    (epsilon^2/4)V(B)+O(epsilon^4)

for this fixed finite system. Along the log-cosh path,
epsilon^2=tau_0(1-u)/2+O((1-u)^2). The left endpoint derivative for
these children and B is therefore

    tau_0/(2tanh(2lambda_child)) sum_internal d_e
                                      -tau_0 V(B)/8. (5)

At a generic temperature where active child pressure derivatives agree,
the optimized envelope selects a bridge minimizing V among its active
child pairs, including both possible relative child polarities. With
derivative ties not resolved, it maximizes the full displayed derivative
over active endpoint choices. This follows from taking the minimum at
u=1-delta, not from assuming a single optimizer stays active.

For fixed children, a sufficient nonpositive endpoint derivative is

    min_(B full sign) V(B)
       >= 4 sum_internal d_e/tanh(2lambda_child).     (6)

This is a concrete energy-versus-optimized-flat-bridge variance lemma.
The root's Hubbard theorem proves a positive constant times m^2 lower
bound on V, but not the coefficient required by (6).

The paired edge inequality bounds the right side of (6) from ABOVE by

    2m^2-Tr(C_1^+ C_1^-)-Tr(C_2^+ C_2^-).

Replacing it by this upper bound does not yield a valid generic PSD
lemma. For example take all four matrices equal to
C_rho=(1-rho)I+rho 11^T and choose a full-sign B whose rows and columns
sum to zero. Then V(B)=2(1-rho)^2 m^2, whereas the displayed upper bound
is 2m(m-1)(1-rho^2), which is larger for many rho and m.
This example only refutes the strengthened generic Gram inequality:
when C^+=C^- the true signed energy d is zero, so it does NOT refute
the actual energy comparison (6). Any successful use of (6) must keep
the actual response/energy together with the common covariance.
