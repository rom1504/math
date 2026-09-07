# Entrywise-small actual surgery and a full covariance balance law

Date: 2026-09-07. Original version independently audited; target-contraction
improvement independently reconstructed. This is a uniform actual
same-order construction, not yet a cross-order flatification theorem.

Let A be any hollow full signing of order n, and let
G be ANY real hollow symmetric matrix with |G_ij|<=1. Put
0<theta<=sqrt(n), a=(n+2)log2, and t=theta/sqrt(n).

Independently flip edge ij with probability

    p_ij=t(1+A_ij G_ij)/(2(1+t)).

Then 0<=p_ij<=t/(1+t)<=1/2 and the exact mean matrix is

    Abar=(A-tG)/(1+t).

The mean is entrywise feasible and Q(Abar)=Q(A-tG)/(1+t). Since t<=1,
|A_ij-tG_ij|>=1-t, so each centered edge variance is at most
1-((1-t)/(1+t))^2=4t/(1+t)^2. Thus total centered variance is at most
2t n(n-1)/(1+t)^2, and each centered coefficient has absolute value
at most2. Bernstein and a union bound over all Boolean spins give an
actual full signing A' satisfying

    Q(A')<=Q(A-tG)/(1+t)
               +sqrt(4t n(n-1)a)/(1+t)+(4/3)a.           (1)

The existence probability for this cap bound is at least1/2. Markov also
allows selection of a realization changing at most2t n(n-1) edges, with
both events holding simultaneously. The edit bound can be vacuous when
t is not small; the actual coefficient construction remains valid.

In normalized units, write

    e_n(theta)=2sqrt(theta (1-1/n)a/n)n^(-1/4)/(1+t)
        +4a/(3n^(3/2)).                                  (2)

Thus for fixed theta, e_n(theta)=O(n^(-1/4)). More generally theta=o(sqrt(n))
still gives e_n(theta)=o(1). No rank, PSD, incoherence, or input cap bound
is used. The independent adversarial audit supplied this target-contraction
improvement. The initial valid proof instead used p=t(1+A G)/2, paid the
attenuation bias t Q(A), and omitted the denominator in the noise term;
that version required an input bound Q(A)<=C n^(3/2).

Suppose A has normalized excess delta above M_n. Global optimality applied
to the actual signing A' proves, uniformly for ALL entrywise-bounded G,

    Q(A-tG)>=Q(A)-(delta+e_n(theta))n^(3/2).               (3)

Write h(x)=H_A(x)/n^(3/2), q=Q(A)/n^(3/2). Minimizing the left side of (3)
over the compact convex entrywise cube for G, and applying finite minimax,
gives a probability law mu on signed spins (sigma,x) with

    E_mu[q-sigma h(x)]
      +(theta/n^2) sum_(i<j)|E_mu[sigma x_i x_j]|
          <=delta+e_n(theta).                            (4)

All terms in (4) are nonnegative. For fixed theta and exact minimizers this
is an O(n^(-1/4)) normalized ground-slack and averaged entrywise covariance
bound. The measure may depend on theta. Taking theta=n^(1/6), for example,
gives mean slack O(n^(-1/6)) and

    n^(-2) sum_(i<j)|E_mu[sigma x_i x_j]|=O(n^(-1/3)).

The omitted diagonal is essential: it means the theorem does not control
E sigma, nor force equal polarity probabilities. In particular a
one-polarity ensemble with almost orthogonal spins is compatible with (4).

If the negative cap is at least gamma n^(3/2) below Q(A), gamma>0, then
mu{sigma=-1}<=epsilon/gamma, where epsilon=delta+e_n(theta). Writing
C_+=E[1_(sigma=+1)xx^T], the off-diagonal triangle inequality gives

    n^(-2)sum_(i<j)|(C_+)_ij|
       <=epsilon/theta+epsilon/(2gamma).

This is a genuine broad-spread constraint on a one-sided near-ground law,
but it supplies no upper bound on max_x(|H_A(x)|+|v dot x|), and no
cross-order cap inequality is claimed.
