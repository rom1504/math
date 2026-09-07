# Cross-order low-cap counting: valid extension and its scale obstruction

Status: elementary proved comparison, not a convergence result. The
same-order positive-entropy thickening is already archived in
`fresh_entropy_support_audit_2026_09_05.md`; the comparison below explicitly
changes matrix order. No novelty is claimed for the random-padding method.

Write d_k=k(k-1)/2, q_k(A)=Q(A)/k^(3/2), and
N_k(C)=#{A: q_k(A)<=C}, counting labeled hollow signings.

## 1. A finite counting inequality

For 2<=n<N, h=N-n and r=n/N,

    N_N(r^(3/2) C + 2 sqrt(h/N))
       >= 2^(d_N-d_n-1) N_n(C).                         (1)

Indeed, fix each good n by n principal block and independently fill all
remaining H=d_N-d_n edges by fair signs. Denote the zero-padded random
remainder by E. For each projective signed spin test, Hoeffding gives

    P(s q_E(x)>u) <= exp[-u^2/(2H)].

There are 2^N such tests. Since H<=hN, taking u=2N sqrt(h) gives

    P(Q(E)>u) <= exp[-(2-log 2)N] < 1/2.

The triangle inequality and the disjointness of extension cylinders give
(1). In particular, where the counts are nonzero, normalized logarithms
s_k(C)=log N_k(C)/d_k satisfy

    s_N(r^(3/2) C + 2 sqrt(1-r))
       >= (d_n/d_N) s_n(C)
          +(1-d_n/d_N)log 2-log 2/d_N.                 (2)

Applying the entropy-compatible core/refill theorem at the output order
additionally makes all counted matrices have operator norm at most
(4 K_G C'/epsilon+8)sqrt(N), at cap threshold C'+2sqrt(epsilon)
and an additional count factor 2^(-N-1), where
C'=r^(3/2)C+2sqrt(1-r). Thus operator regularization is compatible with
this actual cross-order comparison, not merely with one witness.

## 2. Why this does not yield a thermodynamic limit

If h/N=eta is small, the threshold change in (1) is

    2 sqrt(eta) - (3C/2)eta + O(eta^2),

whereas the entropy gain in (2) is only of order eta. Subdividing a
fixed multiplicative change of order into k equal small steps increases
the leading cumulative cap loss like sqrt(k), rather than reducing it.
The bound therefore does not compare distant orders at arbitrarily small
fixed cap relaxation. Positive entropy in a noisy Hamming neighborhood
does not repair that scale mismatch by itself.

This is not a proof that the TOTAL-cap threshold in (1) is optimal.
The square-root scale is, however, unavoidable for a bound based on the
norm of the cross-edge remainder alone: for every n by h sign rectangle R,

    ||R||_(infinity->1) >= n E|S_h|,

by averaging its bilinear test over the h independent column spins and
then optimizing all row spins. Since the absolute quadratic cap of the
zero-internal bridge is precisely this rectangle norm, for h tending to
infinity it is at least (sqrt(2/pi)+o(1)) n sqrt(h). A sharper useful
TOTAL-cap estimate must therefore exploit cancellation with the parent
landscape, not improve independent-noise concentration alone.

## 3. Retention/noise tradeoff does not remove the obstruction

One can also flip each old edge independently so that its mean is t times
the seed, with 0<=t<=1, and fill new edges fairly. The total edge variance
is exactly d_N-t^2 d_n. The strong-variance smoothmax estimate yields

    E q_N(output)
      <= t r^(3/2) C
         + sqrt[2 log(2) (d_N-t^2 d_n)/N^2].           (3)

The right side is a concave function of t on [0,1], hence its minimum is
at t=0 or t=1. Thus this elementary combined channel cannot beat both
full retention and complete resampling. Equation (3) is an EXPECTATION
bound, not by itself the same support-count statement as (1).

An actual near-critical cross-order count inequality would require a
new cancellation or optimized-bridge theorem. This file does not assert
one and does not infer convergence from entropy regularization.
