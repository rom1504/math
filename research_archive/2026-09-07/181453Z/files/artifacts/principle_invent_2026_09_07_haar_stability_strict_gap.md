# A strict Haar cap improvement from one-spin stability

2026-09-07. **PASS: independently reconstructed by the director**, in
`principle_director_haar_stability_audit_2026_09_07.md`. The theorem
is for a Haar balanced symmetric involution, NOT an original full-sign
construction. Its purpose is to test the proposed exact Gaussian-control
duality at the shared Gaussian/annealed boundary.

Let n be even and

    U=O diag(I_(n/2),-I_(n/2)) O^T,

with O Haar orthogonal. Define

    c_stab = (1/2) sqrt(1-exp(1/25000)/16).

Then, in probability and in expectation,

    limsup max_(x in {+-1}^n) |x^T U x|/(2n)
       <= c_stab < sqrt(15)/8.                       (1)

The exact expression is about .48412227276536227685, whereas
sqrt(15)/8 is about .48412291827592711065. Only the exact expression
and strict inequality are used in the proof; these decimals are displays.

## 1. Exact one-vector geometry

By sign-switching invariance, fix the Boolean vector x=1. Put

    r=x^T U x/n,      e=x/sqrt(n).

The variable (1+r)/2 has Beta(n/4,n/4) law. Consequently, for each fixed
0<r_0<1,

    Pr(r>=r_0) <= exp[-n I(r_0)+O(log n)],
    I(r)=-log(1-r^2)/4.                              (2)

For example this follows directly from its density
C_n (1-r^2)^(n/4-1), with C_n=exp(O(log n)).

Conditioned on r, invariance under every orthogonal transformation fixing
e shows that

    Ue=r e+sqrt(1-r^2) v,

where v is uniform on the unit sphere in e-perpendicular. Its conditional
law does not depend on r. Equivalently, for iid standard Gaussians g_i,

    sqrt(n) v_i=(g_i-gbar)/sqrt(q),
    gbar=n^(-1)sum_i g_i,
    q=n^(-1)sum_i(g_i-gbar)^2.                        (3)

No independence of the entries of U is assumed. Equation (3) is the
exact conditional spherical law, not a coordinatewise Gaussian ansatz.

## 2. A uniform exponentially small stability probability

Fix

    eta=1/1000,       r_1=31/32,
    a=1/100,         b=101/100,
    delta=1/50000.

If a Boolean x is a local maximum of the hollow quadratic associated
with C=U-diag(U), then flipping its i-th coordinate gives

    x_i(Ux)_i >= U_ii,              every i.         (4)

On the global event max_i|U_ii|<=eta, (4) implies x_i(Ux)_i>=-eta.
For x=1 and 0<=r<=r_1, equations (3)--(4) imply

    (g_i-gbar)/sqrt(q) >= -M,       every i,
    M=(r_1+eta)/sqrt(1-r_1^2).

On |gbar|<=a and q<=b^2 this forces g_i>=-a-bM for every i. The
chosen rational constants satisfy a+bM<4. Indeed

    M=(31+32/1000)/sqrt(63)
      < (311/10)/(79/10)=311/79<399/101=(4-a)/b.

Therefore, uniformly for all r in [0,r_1],

    Pr(stability inequalities (4), diagonal event | r)
      <= Phi(4)^n + Pr(|gbar|>a) + Pr(q>b^2)
      <= 4 exp(-delta n)                              (5)

for all sufficiently large n. Importantly, the diagonal event is used
only to imply a condition on Ux; no conditional independence between
diag(U) and Ux is asserted.

Here are elementary explicit rate checks. Mills' lower bound gives

    1-Phi(4) >= 4 exp(-8)/(17 sqrt(2pi)) > 1/50000.

For the last strict inequality one may use sqrt(2pi)<3 and exp(8)<3000;
then its left side exceeds 4/153000>1/50000. The latter exponential
bound follows, for example, from e<68/25 and (68/25)^8<3000.
The Gaussian mean tail is at most 2 exp(-n/20000).
Finally nq has chi-square law with n-1 degrees of freedom. Writing
x_0=b^2-1=201/10000, its Chernoff rate is at least

    (n-1)(x_0-log(1+x_0))/2
       >= (n-1) x_0^2/[4(1+x_0)].

The coefficient x_0^2/[4(1+x_0)] exceeds 1/11000, so for n>=2 the
displayed rate exceeds n/22000>delta n. These bounds prove (5).

The diagonal event itself has probability tending to one. Every U_ii
has the same Beta-derived law as r, and (2), applied to +/-eta followed
by a union bound over n coordinates, proves this directly.

## 3. Count stable high-energy candidates, not all high-energy vectors

Set

    r_0=sqrt(1-exp(2delta)/16).

Then 0<r_0<sqrt(15)/4<r_1 and, exactly,

    I(r_0)=log 2-delta/2.

For any fixed Boolean x, conditional spherical symmetry and (5) give

    Pr(r in [r_0,r_1], x stable, diagonal event)
      <=4 exp(-delta n) Pr(r>=r_0).

Summing over 2^n vectors and using (2), the probability of any such
candidate is at most exp[-delta n/2+O(log n)], which tends to zero.
For r>r_1, the original tail union bound already gives a vanishing
probability, since I(31/32)>log 2.

A global maximum is always one-spin stable. Thus, on the diagonal
event, no maximum can have normalized Rayleigh value above r_0 with
probability tending to one. The same proof applies to -U, which has
the same balanced-involution law. A further factor two handles both
energy orientations and proves (1) in probability.

Since ||U||op=1, every normalized half-energy lies in [-1/2,1/2].
Boundedness upgrades the assertion to the limsup of the expected cap.

## 4. Hollowing and the fixed-algorithm consequence

Because Tr U=0 exactly, Boolean quadratic energies of C=U-diag(U)
and U are EXACTLY equal: x^T diag(U)x=Tr U=0. There is no diagonal
error in (1). For arbitrary cube-valued outputs, max|diag U|=o(1)
gives the usual o(n) hollowing error.

Consequently, wherever the already audited fixed-GFOM comparison to
this Haar model applies, its ceiling sqrt(15)/8 can be replaced by
c_stab. That is a ceiling on the specified fixed algorithms, not an
upper bound on the Boolean cap of arbitrary Hadamard signings.

## 5. Exact Gaussian-boundary diagnosis

For the Gaussian self-transport terminal potential g_t(1),

    t+g_t(1)
      =sup_(|r|<1) {t r+(1/4)log(1-r^2)}.

This is precisely the limiting log-moment-generating function of the
fixed-vector Haar Rayleigh variable r. Counting all 2^n spin vectors
and optimizing its Chernoff bound gives

    inf_(t>0) [t+log 2+g_t(1)]/(2t)=sqrt(15)/8.

The optimum is t=2sqrt(15), with r=sqrt(15)/4. Thus the shared
Gaussian Bellman floor and old Haar fixed-algorithm ceiling have the
same number because they use the SAME all-vector annealed tail.
The strict improvement (1) shows this coincidence is not an exact
ground-state primal/dual identity, even on the concrete Gaussian/Haar
boundary model.

The missing constraint is explicit: a maximizing spin vector must
jointly satisfy x_i(Cx)_i>=0 at every coordinate. It is not encoded
by the scalar Rayleigh rate or by a one-source mean-variance channel.
An upper control that retains the joint local-field stability event
is therefore a genuine additional operation, unlike reoptimizing the
same terminal Gaussian functional.

Classical ROM metastable-state literature studies this type of constraint,
including [Degli Esposti--Giardina--Graffi](https://arxiv.org/abs/cond-mat/0207681)
and [Cherrier--Dean--Lefevre](https://arxiv.org/abs/cond-mat/0212246).
No replica saddle, metastability enumeration formula, or numerical
ground-state estimate from those papers is used in this elementary proof.
