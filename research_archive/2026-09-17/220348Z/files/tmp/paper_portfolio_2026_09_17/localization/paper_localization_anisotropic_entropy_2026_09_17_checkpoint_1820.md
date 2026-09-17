# Flat-mode Boolean entropy from anisotropic Gaussian observation

2026-09-17. New deduction of the localization track, independently
reconstructed by the director. The source mechanism is the Gaussian
observation/information decomposition reconstructed from
[El Alaoui--Montanari, arXiv:2109.00709](https://arxiv.org/abs/2109.00709)
in [the primary localization artifact](paper_localization_2026_09_17.md).
The theorem below is proved directly; it is not attributed to that paper.
External novelty is not established.

The key gain is an explicitly paid coarse revelation of the mode-energy
profile, followed by anisotropic Gaussian noise with CONSTANT PHYSICAL
DIAGONAL. No posterior is replaced by independent coordinates.

## 1. Finite theorem

Let P_1,...,P_k be mutually orthogonal projections resolving the identity
on R^n, each of rank p=n/k. Assume the flat-diagonal condition

    (P_a)_(ii)=1/k for every a and every physical coordinate i. (1)

For Boolean x define

    R(x)=[sum_a ||P_a x||]^2/n,      1<=R(x)<=k.

Let F be any NONEMPTY Boolean set satisfying R(x)/k<=theta for every
x in F. For every zeta>0 and delta>0,

    log|F| <= log binom(floor(k/zeta)+k,k)
       +n (sqrt(theta)+sqrt(zeta))^2/(2delta)
       +n h(Phi(-1/sqrt(delta))).                      (2)

Here Phi is the standard normal distribution function and h is binary
entropy in nats. The first term is a genuine revelation cost. It may
not be omitted in a finite-dimensional assertion.

Hadamard-mode projections
`P_a=(h_a h_a^T/k) tensor I_p` satisfy (1). A common coordinate
permutation or coordinate sign-switching preserves (1). The finite
entropy theorem also applies to other flat orthogonal decompositions;
the later full-sign bridge realization specifically uses Hadamard modes.

## 2. Exact channel and proof

Take X uniform on F and let

    e_a(X)=||P_a X||^2/n,       sum_a e_a=1.

Reveal the histogram label J=(J_1,...,J_k), where

    J_a=floor(k e_a/zeta),   u_a=(J_a+1)zeta/k.

There are at most binom(floor(k/zeta)+k,k) possible labels, by the
weak-composition count. Conditional on J, each e_a<=u_a and

    S_J:=sum_a sqrt(u_a)
       <=sum_a sqrt(e_a)+sqrt(k zeta)
       <=sqrt(k)[sqrt(theta)+sqrt(zeta)].              (3)

All u_a are strictly positive. Given J, choose independent Gaussian
noise G_J with covariance

    Gamma_J=sum_a sigma_a^2 P_a,
    sigma_a^2=delta k sqrt(u_a)/S_J,

and observe Y=X+G_J. The noise is independent of X CONDITIONAL ON J;
the profile-dependent choice has already been paid for by H(J).
Since sum_a sigma_a^2=delta k, (1) implies

    (Gamma_J)_(ii)=delta for EVERY i and J.             (4)

Its trace is delta n and it is positive definite. Let Sigma_J denote
Cov(X|J). The ordinary Gaussian-channel information bound gives

    I(X;Y|J)
      <=(1/2)E logdet(I+Gamma_J^(-1/2)Sigma_J Gamma_J^(-1/2))
      <=(1/2)E tr(Gamma_J^(-1)Sigma_J)
      <=(n/2)E sum_a u_a/sigma_a^2
      =(n/(2delta k)) E S_J^2
      <=n(sqrt(theta)+sqrt(zeta))^2/(2delta).           (5)

The first inequality follows by maximizing differential entropy at a
fixed covariance; the next is logdet(I+M)<=tr M. Conditional means
can only reduce the last covariance trace. These are finite exact
inequalities, not an isotropy assumption.

The EXPLICIT decoder sign(Y_i) has error probability exactly
Phi(-1/sqrt(delta)) for every coordinate, even conditional on X and J,
because X_i is a sign and (4) fixes its Gaussian noise variance. It
need not be the posterior MAP decoder. Correlations between the noise
coordinates do not affect this marginal error calculation. Binary
conditional-entropy subadditivity gives

    H(X|Y,J)<=n h(Phi(-1/sqrt(delta))).                 (6)

Finally J is a deterministic function of X, so

    H(X)=H(J)+I(X;Y|J)+H(X|Y,J).

Combining the label count, (5), and (6) proves (2). The posteriors of
X remain supported on the original set F throughout.

## 3. Asymptotic envelope and its sharp leading form

Define the explicit function

    Psi(theta)=inf_(delta>0)
       [theta/(2delta)+h(Phi(-1/sqrt(delta)))]
      =inf_(t>0)[theta t+h(Phi(-sqrt(2t)))].             (7)

If k=k_n=o(n), then for every FIXED theta>0 and every sequence of
sets as in Section 1,

    limsup_n log|F_n|/n <= Psi(theta).                  (8)

To verify the order, fix zeta and delta first. The histogram term is
at most k log(e(1+1/zeta)), hence is o(n). Send n to infinity, then
zeta down to zero, and finally take the infimum over delta. No uniform
estimate for a vanishing zeta_n is presumed.

As theta decreases to zero,

    Psi(theta)=(1+o(1)) theta log(1/theta).              (9)

Here is an elementary normalization check. Put L=log(1/theta) and
choose t=L+2log L for L>=2. The Gaussian Chernoff bound gives
Phi(-sqrt(2t))<=exp(-t)=theta/L^2, whence

    Psi(theta)
      <=theta(L+2log L)
          +(theta/L^2)(L+2log L+1).                    (10)

This proves the upper asymptotic with leading constant one.

For the lower asymptotic, integration over the interval
[a,a+1/a] shows Phi(-a)>=c exp(-a^2/2)/a for a>=1, with an absolute
c>0. Thus h(Phi(-sqrt(2t)))>=c' sqrt(t)exp(-t) for t>=1.
Fix any nu>0. If t>=(1-nu)L, the first term in (7) is at least
(1-nu)theta L. If 1<=t<(1-nu)L, the second term divided by theta L
tends uniformly to infinity. For t<1 it is bounded below by a positive
constant. Taking the infimum, then nu down to zero, proves (9).

Formula (9) describes the envelope (7); it does NOT assert that every
Boolean code attains that entropy. For n-dependent theta_n, the finite
histogram cost and the choice of zeta_n must still be accounted for.

### 3.1 What flat physical diagonals improve

For arbitrary orthogonal equal-rank projections without (1), the same
noise allocation has trace delta n. Squared-error sign decoding then
only gives average bit error at most delta. For 0<delta<=1/2 the same
proof yields the weaker envelope

    inf_delta [theta/(2delta)+h(delta)]
                       <=(1+o(1))sqrt(theta log(1/theta)).

The exact flat-diagonal Gaussian tail in (4)--(6) is therefore essential
to the theta log(1/theta) gain; replacing it by a trace-only estimate
would lose a square root in this comparison.

## 4. All-order leftovers and direct augmentation consequence

For arbitrary n choose k_n=o(n), let n0=k_n floor(n/k_n), and leave
ell=n-n0<k_n coordinates aside. Apply (2) to the n0-dimensional main
word and count the leftovers by at most 2^ell. Its normalized entropy
cost vanishes. The mean-response bound for the balanced full-sign
bridge is still

    E B_x <=(q+k)sqrt(n0 R(x_main)/k)+q sqrt(ell).

Thus its normalized response coefficient on R/k<=theta is at most
sqrt(theta)+o_n(1). Every physical entry is a sign.

Suppose a sequence A_n admits selected flat Hadamard-mode groupings
with k_n=o(n), and for each fixed small eta>0 its FULL near-extreme set
E_n(eta) satisfies

    limsup_n max_(x in E_n(eta)) R(x_main)/k_n <= theta(eta),
    theta(eta)->0,
    K_theta=limsup_(eta down to0)
                 theta(eta)log(1/theta(eta))/eta <infinity.       (11)

Then the actual full-sign extension conclusion in
[the director's entropy-response theorem](paper_director_spectral_entropy_response_2026_09_17.md)
holds for every slope tau>K_theta/2. Indeed (8)--(9) give its entropy
hypothesis with K_s<=K_theta, while the mean-response hypothesis follows
from sqrt(theta(eta))->0.

For each FIXED eta, first take n to infinity with fixed histogram
precision and channel noise, then refine that precision in the entropy
upper bound. Only afterwards send eta to zero. In the finite hierarchy
used to construct a parent, finitely many precisions can all be chosen
before n tends to infinity. Small strict envelope slacks handle a
limsup mode bound; Psi is continuous from the right, including at zero.

In particular, theta(eta)<=[a+o(1)]eta/log(1/eta) permits every
tau>a/2. On a liminf-realizing minimizing sequence this would force
a>=3c_*, by the director's exact-child extension argument. These are
CONDITIONAL landscape consequences. No such full near-level profile
is currently established for actual minimizing matrices; the existing
Hadamard stress test fails the needed concentration hypothesis.

## 5. Independent audit and replay

The director and discrepancy researcher independently reconstructed the
finite channel, information bound, exact marginal decoder error, both
sides of the envelope asymptotic, and the stated limit order; all passed.

The replay is
`computations/paper_localization_2026_09_17_anisotropic_entropy.py`, with
output `tmp/paper_portfolio_2026_09_17/localization/anisotropic_entropy_audit.json`.
It enumerates twelve small finite channel cases, computes their actual
histogram entropy and conditional covariance capacity bounds, and checks
that every physical noise diagonal is delta (maximum numerical discrepancy
1.12e-16). It also numerically minimizes the explicit envelope at seven
theta values from 1e-2 through 1e-12; the ratio to theta log(1/theta)
decreases from about 1.2001 to 1.0549, consistent with the analytic limit.
Membership cutoffs, logarithms, and Gaussian integrals in the replay are
floating-point diagnostics, not certified proof inputs.

