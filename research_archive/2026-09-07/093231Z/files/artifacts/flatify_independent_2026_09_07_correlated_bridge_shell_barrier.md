# Actual-child shell theorem and an unconditional certificate obstruction

Date: 2026-09-07. This concerns the correlated balanced regular-support
bridge in `flatify_independent_2026_09_07_regular_support_weave.md`.
It does not prove an impossibility for the actual bridge ensemble or for
original convergence. The obstruction is to its uniform fixed-pair
Chernoff estimate followed by a union bound, even with exact child shells
and a separately optimized temperature for each pair of shells.

## 1. Precise sufficient original-value inequality

Let A,D be actual hollow sign children of equal order n, and let B be the
correlated square bridge with n=d floor(pd). Write H_A=x^T A x/2.
The parent has the exact cap identity

    Q(parent)=max_{x,y} (|H_A(x)+H_D(y)|+|x^T B y|).

Indeed replacing y by -y leaves both child energies unchanged and
changes the bridge sign. In this identity no child energies are discarded.

The direct E fixed-pair Chernoff rate available from the weave proof is

    I_p(b)=sup_{t>0} [2t(b sqrt(p)-1)-2 E_t(nu_p)]/p.       (1)

Here the bridge event is |x^T B y|>b n^(3/2); the factor 2 for its two
signs is subexponential. The ensemble is independent of t. Thus a finite
list of temperatures, one for each shell pair, can use a common finite
recursive depth before taking d to infinity.

A precise finite-grid sufficient condition is as follows. Partition the
normalized child energy intervals into finitely many bins J_a,J_d. Let
N_A(a),N_D(d) be their exact spin counts, and let

    b_ad=L-sup{|u+v|: u in J_a, v in J_d}.

If every nonempty pair has b_ad>0 and there are fixed temperatures t_ad
and epsilon>0 such that, for all sufficiently large compatible orders,

    [log N_A(a)+log N_D(d)]/n
      < [2t_ad(b_ad sqrt(p)-1)-2E_{t_ad}(nu_p)]/p-epsilon,  (2)

then a bridge exists with Q(parent)<=L n^(3/2). First approximate each
of the finitely many E values by sufficient finite recursive depth, then
apply the fixed-pair estimate, and sum over bins and their actual spins.
The sum tends to zero exponentially. This is an actual sign-parent
construction criterion, not a width or Gaussian-disorder surrogate.

## 2. Universal upper bound on this available rate

The uninformative scalar channel gives E_t(nu_p)>=g_t(1), where

    g_t(1)=-t(1-rho)+(1/4)log(1-rho^2),
    t=rho/[2(1-rho^2)],  0<rho<1.

Taking the Legendre transform explicitly, for 0<=b<1 and 0<p<=1,

    I_p(b)<=J_p(b):=-log(1-p b^2)/(2p)
                  <=J_1(b)=-log(1-b^2)/2.                (3)

At the maximizer rho=b sqrt(p), the t-dependent linear terms cancel.
The last inequality follows either by differentiating in p or expanding
-log(1-p b^2)/p as its nonnegative power series. Equation (3) upper-bounds
the rate supplied by this certificate; it is NOT a lower bound on actual
tail probabilities. Better ensemble-specific estimates could evade it.

## 3. Large actual shells from noisy extrema

Let x* maximize sigma H_A for sigma in {+1,-1}, with value q n^(3/2).
Switch and orient A so x*=1 and H_A(1)=q n^(3/2). Let independent eta_i
equal -1 with probability delta in (0,1/2), and put alpha=1-2delta.
Writing r_i=sum_j A_ij after this switching, local optimality gives
r_i>=0, while sum_i r_i=2q n^(3/2) and r_i<=n-1. Direct expansion gives

    E H_A(eta)=alpha^2 q n^(3/2),
    Var H_A(eta)=alpha^2(1-alpha^2)sum_i r_i^2
                   +(1-alpha^2)^2 binom(n,2)
                <=2alpha^2(1-alpha^2)(n-1)q n^(3/2)
                   +(1-alpha^2)^2 binom(n,2).            (4)

For uniformly bounded normalized caps, normalized variance is O(n^-1/2).
Consequently there is a sequence of energy windows of width o(n^(3/2))
around sigma alpha^2 q n^(3/2), each containing at least

    exp[n h(delta)-o(n)]                                 (5)

actual spins. To verify counting, intersect the probability-1-o(1)
energy event with |number of flips-delta n|<=n^(3/4). Every string in the
intersection has noise probability exp[-n h(delta)+o(n)], uniformly.
Its cardinality is therefore at least (5). No spectral norm assumption,
ground-state uniqueness, or probabilistic hypothesis on A is needed.

## 4. Width removes every child-polarity hypothesis

Put q_A^+=max H_A/n^(3/2), q_A^-=max -H_A/n^(3/2), and similarly D.
The universal half-range theorem in
`decisive_audit_certified_minimum_width_lower_2026_09_07.md` gives

    (q_A^++q_A^-)/2 >= ell-o(1),
    ell=.4333221116640807,

for every signing sequence, not only width minimizers. For either
relative child polarity eta, choose sigma with

    q_A^sigma+q_{eta D}^sigma
      >= [(q_A^++q_A^-)+(q_D^++q_D^-)]/2
      >=2ell-o(1).                                      (6)

Both noisy extremal clouds then have the same energy sign. Their product
has entropy at least 2h(delta)-o(1), and combined absolute child energy
at least 2alpha^2 ell-o(1). This avoids assuming either child has balanced
extrema. The choice of polarity may depend on n and on the children.

For actual optimizing children, the all-order upper bound gives
Q(A),Q(D)<(U+o(1))n^(3/2), where U=.493608094. A constant-preserving
equal-child closure seeks L<=2sqrt(2)U+o(1). Take delta=1/10. In the
cloud pair from (6), the admissible bridge slack is at most

    b_*=2sqrt(2)U-(32/25)ell = .8414822191338438...,

whereas (3) and (5) give

    available rate <=J_1(b_*)=.6156588543111433...,
    actual pair entropy >=2h(1/10)=.6501659467828965....    (7)

The gap exceeds .03 rigorously. Directed 50-digit interval arithmetic,
with rational inputs U=493608094/10^9 and ell=4333221116640807/10^16,
is in `computations/flatify_independent_2026_09_07_shell_barrier_certificate.py`
and its matching results JSON. Its exact rational interval endpoints,
not displayed decimals, certify the gap.

Thus (2) fails for actual optimal children, regardless of their relative
polarity, for every balanced retention parameter p in (0,1]. It fails
even if every shell count is known exactly and each shell pair receives
its own best temperature in the direct E rate. Shrinking energy windows
or refining finitely many bins cannot repair the fixed positive gap.

## 5. Scope and operative consequence

The sufficient theorem (2) is valid but does not deliver constant-
preserving closure. The obstruction identifies an actual, unavoidable
intermediate-energy family, rather than a hypothetical entropy profile.
It does not forbid a successful correlated bridge: its pairwise tail
events are strongly dependent, and a shell-conditioned bound or a
construction correlating the bridge with the child landscape could be
much stronger than the uniform fixed-pair estimate used here.

In particular this is not a lower bound on the cap of all regular-support
weaves, not a counterexample to the desired parent inequality, and not a
proof of convergence or nonconvergence of M_n/n^(3/2).
