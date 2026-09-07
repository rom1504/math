# Dense fair bridges have a leading floor, even for actual optimal children

Date: 2026-09-07. Status: **proved; root and independent audits passed**; see
`flatify_independent_2026_09_07_iid_bridge_control_audit.md`.
This is a scalable obstruction to typical independent fair bridge filling,
uniform over the internal children. It does not rule out an exceptional
deterministic bridge, exponentially rare sampling, or correlated rounding.

Subsequent stronger explicit control: the independently reconstructed
`flatify_independent_2026_09_07_heat_martingale_sk_bound.md` improves the
bridge floor to Gamma(3/4)/[sqrt(2pi)Gamma(5/4)] = .539352601188...,
with a rational certificate already proving a floor above .53. The
stopped-exit proof below remains valid and records the original derivation.

## 1. Statement and relevance to selected minimizing children

Let B be an m-by-m matrix of independent fair signs. Define
beta(B)=max_(x,y Boolean)|x^TBy|. Then

    liminf_m E beta(B)/(2m)^(3/2)
       >=2240/[57 sqrt(2) pi^(7/2)] > .505.              (1)

The weaker rational shortcut .71/sqrt(2)>.502 also suffices below.

Moreover beta(B)>.5(2m)^(3/2) with probability 1-exp(-Omega(m)).
For ANY internal hollow signings A,D, including choices made after B,

    Q([[A,B],[B^T,D]])>=beta(B).                          (2)

Indeed reversing all spins in the first block preserves both internal
energies and reverses the bridge response; the larger absolute response
is at least the bridge response alone.

Actual equal-order minimizing children have favorable target

    2sqrt((2m-1)/(m-1)) M_m
       <=(.493608094+o(1))(2m)^(3/2).

Thus fair iid bridge filling misses that target by a fixed leading gap
with exponentially high probability, uniformly over those ACTUAL
children and their polarities. The conclusion is not extrapolated from
Hadamard or other suboptimal child families.

## 2. A certified analytic SK lower bound, rather than an approximate constant

Use the standard SK normalization with Gaussian covariance
N xi(q), xi(q)=q^2/2. Write P_SK for its limiting maximum divided by N.
The only imported spin-glass results are the ground-state Parisi formula
and its stochastic control representation, respectively Theorem 1 and
Corollary 2, equation (12), of
[Auffinger–Chen, Parisi formula for the ground state energy in the mixed
p-spin model](https://arxiv.org/pdf/1606.05335). Their definitions identify
SK as xi(q)=q^2/2, so xi''=1. The relevant statements and normalization
were read directly. In particular

    P_SK=inf_gamma [Psi_gamma(0,0)-(1/2)integral_0^1 s gamma(s) ds],
    Psi_gamma(0,0)=sup_(adapted |u|<=1)
       E|W_1+integral_0^1 gamma(s)u_s ds|
       -(1/2)integral_0^1 gamma(s) E u_s^2 ds.

The formulas hold for nonnegative nondecreasing integrable gamma; the
paper's L1 extension covers the full class. No numerical Parisi value
is used below.

Choose a bounded martingale u_s, u_0=0, with terminal sigma=u_1 in
{+1,-1} and E u_s^2=s. Since u_s=E[sigma|F_s], the inequality |z|>=sigma z
gives, for EVERY gamma,

    Psi_gamma(0,0)-(1/2)integral s gamma(s) ds
       >=E[sigma W_1]
          +integral gamma(s)[E(sigma u_s)-(Eu_s^2+s)/2]ds
       =E[sigma W_1].                                  (3)

The same martingale works for all gamma. Thus it supplies a uniform lower
bound before taking the infimum, with no reversal of a variational bound.

### Explicit variance-parametrized stopped Brownian martingale

Let tau be the exit time of standard Brownian motion from (-1,1), starting
at zero, and p(t)=P(tau>t). Its mean is E tau=1. Put
s(t)=E[min(tau,t)]=integral_0^t p(v)dv, which increases from zero to one,
and let t(s) be its inverse. Run Brownian motion at deterministic clock
t(s) and stop upon reaching +/-1. Equivalently, on the control Brownian
filtration the resulting martingale solves

    du_s=sqrt(t'(s)) 1_(not yet stopped) dW_s.

It is bounded, tends to a terminal sign, and E u_s^2=s by optional
isometry. Although t'(s) diverges near one, its stopped stochastic
integral has total expected quadratic variation one, so the construction
and the covariance calculation are legitimate. Ito covariance gives

    E[sigma W_1]=integral_0^1 sqrt(t'(s))p(t(s)) ds
                 =integral_0^infinity p(t)^(3/2)dt.      (4)

### Elementary explicit integral lower bound

The Dirichlet heat equation on (-1,1), with initial survival function one,
gives at the origin

    p(t)=(4/pi)sum_(k>=0) (-1)^k/(2k+1)
                     exp(-(2k+1)^2 pi^2 t/8).

The alternating terms decrease in magnitude. With c=pi^2/8,

    p(t)>=(4/pi)e^(-ct)[1-(1/3)e^(-8ct)].

Using (1-z)^(3/2)>=1-(3/2)z for 0<=z<=1 yields

    integral p(t)^(3/2)dt
      >=(4/pi)^(3/2)[1/((3/2)c)-(1/2)/((19/2)c)]
      =2240/[57 pi^(7/2)] > .71.                       (5)

The last inequality is certified by pi<22/7 and squaring positive sides:

    2240^2 * 7^7 * 100^2 - 71^2 * 57^2 * 22^7
       =468978757537408 > 0.

Thus (3)--(5) prove P_SK>.71. The finer numerical integral near .727 is
not needed for any claim in this note.

The rational-only replay
`computations/flatify_adversary_2026_09_07_iid_bridge_constant.py` and its
computations/results JSON certify every strict constant comparison using
integer square roots and exact fractions, without a numerical Parisi solve.

## 3. Gaussian bridge comparison

For independent standard Gaussian g_ij define
X_(x,y)=m^(-1/2)sum g_ij x_i y_j. Its covariance is m q_x q_y.
On the same 2m-spin index set take an SK Gaussian process with covariance

    m[(q_x+q_y)/2]^2.

Both variances equal m, and the SK covariance exceeds the bridge
covariance by m(q_x-q_y)^2/4. The Gaussian comparison inequality therefore
gives E max X >= E max SK_(2m). The exact covariance convention can be
implemented by adding a spin-independent centered Gaussian to the usual
hollow SK Hamiltonian; this does not change its expected maximum.
Consequently

    liminf E beta(G_m)/m^(3/2)>=2P_SK>1.42.             (6)

The Gaussian comparison itself follows from differentiating the expected
log-sum-exp along the Gaussian interpolation, using the equal variances
and the displayed covariance difference, then letting its temperature
increase. No formula for bipartite SK is assumed.

## 4. Transfer to actual fair signs and concentration

For an m-by-m array J use the smooth maximum

    f_beta(J)=(beta m)^(-1) log sum_(x,y)
                           exp(beta x^TJy/sqrt(m)).

Its difference from beta(J)/m^(3/2) is between zero and 2log2/beta.
Each third coordinate derivative has absolute value at most
8 beta^2/m^(5/2), because it is a third centered moment of a +/-1 Gibbs
observable times the displayed scaling. Replacing the m^2 Gaussian
entries one at a time by independent Rademachers, Taylor expansion with
matching first two moments gives expectation difference O(beta^2/sqrt(m)).
Taking beta=m^(1/6) balances this with the smoothing error, giving the
explicit comparison rate O(m^(-1/6)).
Thus (6) transfers to the actual sign bridge without assuming a bipartite
universality theorem.

Changing a single bridge sign changes beta(B) by at most 2. Bounded
differences therefore gives

    P(beta(B)<=E beta(B)-t)<=exp(-t^2/(2m^2)).

Since .71/sqrt(2)>.5, (1) implies the high-probability statement for all
sufficiently large m. Combining it with (2) proves the actual-child scope.

## 5. What any successful correlated bridge law must pay

Let E_m be the bridge event beta(B)<=c(2m)^(3/2), for any fixed c<.5.
The fair bridge law U_m has U_m(E_m)<=exp(-Omega(m)). Hence EVERY law R_m
supported on E_m obeys

    D(R_m || U_m)>=log(1/U_m(E_m))=Omega(m).

This follows by conditioning/data processing. A successful parent law with
cap below that threshold has a bridge marginal supported on E_m, so this
necessary entropy cost applies even if its children are chosen adaptively.
The lower bound is only linear in m, NOT quadratic; it does not exclude
correlated rounding or an exponentially rare bridge selection mechanism.
