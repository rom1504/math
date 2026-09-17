# Low actual cap forces a uniform cheap physical response law

2026-09-17. **Proved**, conditional only on the explicitly linked prior
universal half-range and bounded-spectrum Gaussian-sign theorems.
The new combination was independently reconstructed by all three campaign
researchers. External novelty is not claimed. This is an original-signing
structural theorem, NOT a convergence theorem or an improved scalar bound.

## 1. Statement, including every uniformity requirement

For a symmetric hollow full sign matrix A of order n, write

    H_A(x)=x^T A x/2,   Q(A)=max_x |H_A(x)|,
    W(A)=(max_x H_A(x)-min_x H_A(x))/2,
    kappa=sqrt(2/pi),   K=pi/[2 log(1+sqrt(2))].

Suppose the universal all-order statement

    W(B)>=(ell-o(1)) k^(3/2)                               (1)

holds for EVERY hollow full signing B of order k. Fix constants C,c>0
with

    C < ell+c/(2K).                                       (2)

There are delta>0 and n0, depending only on C,c and the convergence
modulus in (1), such that for EVERY n>=n0 and EVERY A with
Q(A)<=C n^(3/2), there is ONE probability law nu on physical columns
h in {+-1}^n satisfying

    E_nu h=0,   E_nu hh^T=I,
    E_nu exp(v.h)<=exp(3||v||_2^2/4)  for every real v,
    sup_{x: |H_A(x)|>=c n^(3/2)} E_nu |h.x|/sqrt(n)
        <= kappa-delta.                                   (3)

If the indicated high-energy set is empty the conclusion is vacuous.
Otherwise nu is a mixture of paired Gaussian-sign laws with correlation
spectra in [1/2,3/2]. No bound on ||A||op, on a dual query covariance,
or on the cardinality of the high-energy set is assumed. The internal
edges of A are not changed. Exact minimizers are included but are not
required. The law may depend on all of A and need not have a presently
known polynomial-time construction or compressed description.

The inherited [half-range theorem](decisive_audit_certified_minimum_width_lower_2026_09_07.md)
gives ell=.4333221116640807. For the safe choices K<=2, ell=.43,
C=.5, c=.30, the strict margin in (2) is .005. An independently
checked rational calculation in the companion audit gives the fully
explicit choice delta=2^(-67) for all sufficiently large n.
The threshold n0 is not made numerically effective here.

## 2. The signed-energy alternative retained from the paper tracks

We reconstruct the input needed from
[the signed-energy concentration alternative](paper_discrepancy_signed_energy_concentration_alternative_2026_09_17.md).
Let mu be ANY law on the high-energy set and define

    sigma(x)=sign H_A(x),   Sigma=E_mu xx^T,
    Kmu=E_mu sigma(x)xx^T,  s=E_mu sigma(x),
    K0=Kmu-sI,  Qh=1_(Sigma>2),  P=I-Qh,
    a=tr(Sigma Qh)/n,   b=||P K0 P||_F^2/n,
    J=E_mu sigma(x)(Qh x)^T A(Qh x).

If every paired law has mu-average normalized response at least
kappa-o(1), uniformly along a putative sequence of counterexamples,
then

    a=o(1),   b=o(1),
    J=2 E_mu |H_A(x)|+o(n^(3/2)).                           (4)

Here is why this implication contains actual information. For a hollow
symmetric D with ||tD||op<=1/2, take the equal mixture of signs of
N(0,I+tD) and N(0,I-tD). It is exactly isotropic by the arcsine
identity. Uniform bounded-spectrum scalar replacement gives

    E_nu |h.x|/sqrt(n)=kappa f(v_x)+e_(n,x),
    f(v)=(sqrt(1+v)+sqrt(1-v))/2 <= 1-v^2/8,
    v_x=(kappa^2/n)x^T asin(tD)x,
    sup_x |e_(n,x)|<=O(n^(-1/6)sqrt(log(en))).               (5)

Arcsine is entrywise. The reconstruction and hypotheses of this
replacement are in
[the actual-energy response proof](paper_bernoulli_signed_covariance_response_2026_09_17.md),
Section 4, and its linked full proof; matching covariances alone would
NOT justify it.

Taking D=Qh-diag(Qh), t=1/2, the PSD arcsine series and rank(Qh)<=an/2
give E_mu v_x >= kappa^2(1/2-pi/12)a. Thus absence of a fixed response
discount forces a->0. Since -Sigma<=Kmu<=Sigma, B=P K0 P has operator
norm at most three. Taking D=B-diag(B), one has

    <Kmu,D>=||B||_F^2,
    ||asin(tD)-tD||op<=216t^3,  0<t<=1/12.

Therefore a fixed positive b likewise forces a response discount by
(5), with t chosen sufficiently small depending on b. Finally

    |E_mu sigma(Px)^T A(Qh x)|<=sqrt(2a)n^(3/2),
    |s<A,Qh>|<=sqrt(a/2)n^(3/2),
    |<A,B>|<=sqrt(b)n^(3/2).

These estimates use ||A||F^2=n(n-1), not ||A||op. Exact block expansion
proves (4), and quantitatively

    J/n^(3/2) >= 2c-sqrt(b)-(5/sqrt(2))sqrt(a).              (6)

## 3. The additional mechanism: covariance mass localizes to coordinates

Put Y=Qh x and R=E_mu YY^T=Qh Sigma Qh. The spectral projection
commutes with Sigma, so 0<=R<=Sigma. In particular R_ii<=1, even
though individual coordinates Y_i need not lie in [-1,1]. Fix epsilon>0
and set

    S={i:R_ii>epsilon},   T=[n]\S.

Then |S|<=an/epsilon. In the Hilbert space L2(mu), the vectors
u_i=sigma Y_i and v_j=Y_j have norms sqrt(R_ii), sqrt(R_jj).
The real bilinear Grothendieck inequality gives, writing beta for
the Boolean bilinear norm,

    |E sigma Y_T^T A_TT Y_T|<=K epsilon beta(A_TT),
    |E sigma Y_S^T A_ST Y_T|<=K sqrt(epsilon) beta(A_ST),
    |E sigma Y_S^T A_SS Y_S|<=K beta(A_SS).                  (7)

No same-map rounding assertion is used: u and v are different families,
and the Grothendieck constant is explicitly paid. The physical paired
law is constructed separately in (5).

For every principal restriction, independent rounding puts the energy
of any cube point in its Boolean energy interval. Hollow polarization
therefore proves

    beta(A_SS)<=4W(A_SS),   beta(A_TT)<=4Q(A).               (8)

Moreover beta(A_ST)<=Q(A): for fixed x,y, global reversal of y leaves
both internal energies unchanged and reverses the bridge. Consequently
max_{epsilon=+-1}|H_S(x)+H_T(y)+epsilon x^T A_ST y|
=|H_S(x)+H_T(y)|+|x^T A_ST y|. The two cross terms in J must BOTH
be included. Combining (7)--(8) yields

    J <=4K W(A_SS)+K Q(A)[2sqrt(epsilon)+4epsilon].          (9)

This converts a thin abstract covariance subspace into a small
COORDINATE range budget. Neither small rank nor small trace alone
would justify discarding J.

## 4. The actual full support pays for the complementary range

The exact half-range inequality is

    W(A)>=W(A_SS)+W(A_TT),  hence
    Q(A)>=W(A_SS)+W(A_TT).                                 (10)

For each pair of endpoint witnesses, bridge reversal makes the parent
maximum at least the sum of internal maxima, and the parent minimum
at most the sum of internal minima. Subtract and divide by two. No
assumption on either interval's midpoint is needed. This elementary
range-budget fact already appeared in the archive; (7)--(9) supply
the new use of it.

Under (4), for each fixed epsilon>0 we have |S|=o(n). By (1),
W(A_TT)>=(ell-o(1))n^(3/2). Equations (9)--(10) imply

    2c <= 4K(C-ell)+KC[2sqrt(epsilon)+4epsilon]+o(1).

First send n to infinity with epsilon fixed, THEN epsilon to zero.
This contradicts (2). Thus no sequence of dual query laws can approach
the independent response kappa under the stated fixed margin.

Equivalently, for some fixed d>0, every such mu admits a paired law
with average response at most kappa-d. The compact convex hull of
paired laws sits in the finite-dimensional probability simplex on
physical columns. Finite minimax exchanges the mu and nu optimizations
and gives ONE mixed law satisfying (3) simultaneously at every query.
Centering, exact isotropy, and the common subGaussian MGF bound are
preserved under this mixture. The last bound follows from Gaussian
product Holder for correlation R<=3I/2 and cosh(z)<=exp(z^2/2).

## 5. Quantitative version and fixed limit order

Before taking limits, when a<epsilon, (6), (9), and (10) give

    C >= (ell-o(1))(1-a/epsilon)^(3/2)
       +[2c-sqrt(b)-(5/sqrt(2))sqrt(a)]/(4K)
       -C[sqrt(epsilon)/2+epsilon].                        (11)

All constants epsilon,a,b can be selected FIRST using the strict
margin in (2). The uniform replacement error in (5) is then removed
by n->infinity. This supplies an explicit positive discount, not a
diagonal extraction using uncontrolled covariance bounds. A sharper
tracking with epsilon=a^(2/3) gives a margin-to-discount lower bound
of order margin^6, with explicit fixed constants in the independent
localization audit. No useful numerical parent slope is asserted.

## 6. Primary-source reconstruction and verification

The real inequality used in (7) follows from Theorem 4.1 of
[Friedland--Lim--Zhang, arXiv:1711.10595](https://arxiv.org/abs/1711.10595).
The decisive mechanism was reconstructed: put z=asinh(1), lift each
unit vector into the direct sum of its odd tensor powers with weights
sqrt(z^(2j+1)/(2j+1)!), alternating the signs on one of the two
families. Both lifts have squared norm sinh(z)=1, and their cross
inner product is sin(z<u,v>). Gaussian signs turn this into
(2z/pi)<u,v>, since z<pi/2. Averaging the resulting Boolean bilinear
forms proves K=pi/(2z). Vectors of norm at most one are padded in
mutually orthogonal extra directions; finite families reduce to a
finite Gram space. This imports no unverified numerical improvement
to Grothendieck's constant.

The half-range dependency was reread, including its full lower-chain
reconstruction, and the frozen rational lower certificate was replayed
in this campaign. The new finite normalization checks exhaust all
full signings through order five, plus additional samples: 49,736
principal profiles/partitions verify (8) and (10) exactly. See
`computations/paper_bernoulli_2026_09_17_range_localization_audit.py`.
These finite tests are checks of algebra, not proofs of the asymptotics.

## 7. Exact gain and exact remaining gap

Removed obligation: bounded operator norm of the actual signing, or
bounded covariance of every maximizing query law, is NOT needed to
obtain a fixed cheap physical response on the complete high-energy
set of a sufficiently low-cap signing. Large covariance mass is
directly exploitable; thin covariance energy would spend more range
on a negligible coordinate set than the low parent cap allows.

Not removed: the useful near-order slope would need response below
approximately 3Q(A)/(2n^(3/2)), which is at most .741 in the current
regime. A discount of 2^(-67) below kappa=.797884... does not meet
that slope. Nor does a mean-response inequality alone control the
maximum over all old spins after several new rows interact. The
current theorem therefore does NOT improve the reported interval,
establish favorable all-order realization, or settle convergence.

The existing actual-sign high-energy counterexample has cap at least
.5+kappa/4-o(1), outside (2), and is consistent with this theorem.
This distinction between high energy in arbitrary signings and high
energy in low-cap full signings is essential.
