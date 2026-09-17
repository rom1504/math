# Independent reconstruction: low-cap full signings admit a uniform physical response discount

2026-09-17, Bernoulli track. Independent mathematical reconstruction of
the director's [canonical low-cap response theorem](paper_director_low_cap_uniform_response_2026_09_17.md).
This audit records the argument and its imported lower-bound dependency;
the elementary range lemma itself is already archived. It does not
claim a new extremal constant, convergence, or a sufficient parent slope.

## 1. Statement and normalization

Write H_A(x)=x^T A x/2, Q(A)=max|H_A(x)|, and
W(A)=(max H_A-min H_A)/2 for a hollow full signing. Let w be any
certified universal asymptotic minimum-half-range lower constant, so

    liminf_m min_(full signings B of order m) W(B)/m^(3/2)>=w.

Let K be a valid real Grothendieck constant upper bound. Fix C,c>0
with the strict inequality

    C < w+c/(2K).                                      (1)

Then there are delta>0 and n_0, depending only on this fixed data,
such that every order-n full signing with Q(A)<=C n^(3/2), n>=n_0,
admits one centered exactly isotropic physical sign law nu with
subGaussian covariance proxy (3/2)I and

    sup_(|H_A(x)|>=c n^(3/2)) E_nu |h.x|/sqrt(n)
      <=sqrt(2/pi)-delta.                              (2)

An empty query code is harmless. The law may depend on A and on the
chosen threshold. It is a single law for the ENTIRE code, not a law
chosen separately at each query.

For example w=.4333221116640807 and the safe bound K=2 allow C=.5
and c=.4 with considerable strict slack. No assumption on ||A||_op,
the covariance of the query law, or exact minimizing optimality is used.

## 2. Imported paired-sign response and its no-discount consequence

The exact imported scalar mechanism is reconstructed in
[the signed-covariance response theorem](paper_bernoulli_signed_covariance_response_2026_09_17.md)
and [the no-operator-bound concentration alternative](paper_discrepancy_signed_energy_concentration_alternative_2026_09_17.md).
For any symmetric hollow H with ||tH||_op<=1/2, the equal mixture of
signs of Gaussians with covariances I+tH and I-tH is centered, exactly
isotropic, and (3/2)-subGaussian. Uniformly in every Boolean x,

    E|h.x|/sqrt(n)=kappa f(v_x)+O(n^(-1/6)sqrt(log(en))),
    v_x=(kappa^2/n)x^T asin(tH)x,
    f(v)=(sqrt(1+v)+sqrt(1-v))/2<=1-v^2/8,
    kappa=sqrt(2/pi).

The arcsine is entrywise. The error has an absolute constant because
all latent spectra stay in [1/2,3/2]. This is a physical sign theorem,
not replacement of the column by a Gaussian in the final conclusion.

Suppose, toward a contradiction, that (2) fails for all fixed delta
along a sequence n->infinity. Compact finite-dimensional minimax over
the convex hull of these paired laws supplies a query law mu_n on
the complete high-energy code whose average response against EVERY
paired law is at least kappa-o(1). Define

    Sigma=E_mu xx^T, sigma=sign H_A(x),
    K0=E_mu sigma xx^T-(E_mu sigma)I,
    Q=1_(Sigma>2), P=I-Q, M=tr(Sigma Q), B=PK0P.

The audited concentration alternative gives

    M=o(n), rank(Q)=o(n), ||B||_F^2=o(n),
    J:=E_mu sigma(Qx)^T A(Qx)
        =2E_mu|H_A(x)|+o(n^(3/2)).                  (3)

For completeness, the two key quantitative reasons are as follows.
The pair H=Q-diag Q, t=1/2 has mean v at least
kappa^2(1/2-pi/12)M/n, using PSD Schur powers in asin(Q/2).
Thus M/n must vanish. Also ||B||_op<=3; for H=B-diag B,

    <K0,H>=||B||_F^2,
    ||asin(tH)-tH||_op<=216t^3,  t<=1/12.

Any fixed positive subsequential lower bound on ||B||_F^2/n would
therefore give a fixed response discount, contrary to the chosen mu.
Finally the full-unit Frobenius identity ||A||_F^2=n(n-1) bounds

    |E_mu sigma(Px)^T A(Qx)|<=sqrt(2n(n-1)M),
    |(E_mu sigma)<A,Q>|<=n sqrt(rank Q),
    |<A,B>|<=n||B||_F,

which proves the last identity in (3). No A-operator bound enters.

## 3. Small covariance mass localizes the marked energy to few coordinates

Put Y=Qx and R=E_mu YY^T=Q Sigma Q. Because Q is a spectral projector
of Sigma, 0<=R<=Sigma. In particular R_ii<=1 and tr R=M=o(n).
Fix epsilon>0 BEFORE sending n to infinity and let

    S={i:R_ii>epsilon}, T=S^c.

Then |S|<=M/epsilon=o(n), while ||Y_i||_(L2(mu))<=sqrt(epsilon)
for i in T and <=1 for all i.

Apply real Grothendieck to the Hilbert vectors sigma Y_i and Y_j.
The bilinear norm beta(A)=max_(u,v Boolean)|u^T A v| obeys
beta(A)<=4W(A)<=4Q(A). Rectangular and principal restrictions cannot
increase beta, by averaging the omitted Boolean coordinates. Hence

    |E sigma Y_T^T A_TT Y_T|<=4K Q(A)epsilon,
    |E sigma Y_S^T A_ST Y_T|<=4K Q(A)sqrt(epsilon).   (4)

The displayed cross bound is deliberately valid but coarse. Bridge
reversal gives the sharper rectangular norm beta(A_ST)<=Q(A), hence
the cross bound improves to K Q(A)sqrt(epsilon). The canonical proof
uses this stronger version. There are TWO cross terms in the symmetric quadratic expansion.
Combining (3)--(4),

    E sigma Y_S^T A_SS Y_S
      >=[2c-4KC(epsilon+2sqrt(epsilon))-o(1)]n^(3/2).

The same vector inequality on S gives

    E sigma Y_S^T A_SS Y_S<=K beta(A_SS)<=4K W(A_SS).

Consequently

    W(A_SS)/n^(3/2)
      >=c/(2K)-C(epsilon+2sqrt(epsilon))-o(1).        (5)

No assertion is made that Qx is Boolean, that Q is a coordinate
projector, or that its large energy is small merely because M=o(n).
The coordinate set S is selected only AFTER the covariance estimate.

## 4. Half-ranges add, even when the favorable sectors differ

For any principal partition S,T, averaging the two global reversals
of the T spin leaves the two internal energies unchanged and cancels
the bridge. Choosing the internal maxima and, separately, the internal
minima proves

    max H_A>=max H_(A_SS)+max H_(A_TT),
    min H_A<=min H_(A_SS)+min H_(A_TT),
    W(A)>=W(A_SS)+W(A_TT).                          (6)

This is why a one-sided bulk lower bound is unnecessary. The two
range endpoints need not use the same energy polarity.

The other factor check is exact: for Boolean u,v in S, let
a=(u+v)/2 and b=(u-v)/2. They lie in the coordinate cube, and

    u^T A_SS v=2[H_(A_SS)(a)-H_(A_SS)(b)].

Independent coordinate rounding keeps both energies inside the full
Boolean interval, so beta(A_SS)<=4W(A_SS), as used in (5).

Since |T|=n-o(n), the universal width lower bound gives
W(A_TT)>=(w-o(1))n^(3/2). Now (5)--(6) and Q(A)>=W(A) imply

    C>=w+c/(2K)-C(epsilon+2sqrt(epsilon)).

First finish the n-limit at each FIXED epsilon; only then send
epsilon down to zero. This contradicts (1) and proves (2).

## 5. Exact scope of the imported .4333221116640807 width constant

This audit read the complete
[minimum-width proof](decisive_audit_certified_minimum_width_lower_2026_09_07.md)
and the complete
[fresh lower-chain reconstruction](decisive_audit_fresh_full_lower_chain_2026_09_07.md).
Their relevant content is universal, not optimizer-specific:

1. A simultaneous Grothendieck diagonal majorant controlled by width
   permits deletion of a fixed small coordinate fraction, leaving an
   actual principal FULL signing with fixed normalized operator cap.
2. Principal width is monotone, by averaging omitted fair spins.
3. The two literal feasible means F+H sign(BF) and -F+H sign(BF)
   have energy DIFFERENCE 2sqrt(n-1) sum_i H_i |(BF)_i|. This is
   bounded by the complete energy range 2W, without opposite-polarity
   or midpoint assumptions.
4. The same frozen marked-response policy and rational certificate
   therefore give the width constant .4333221116640807. The matrix
   and approximation limits are completed at fixed deletion fraction
   before that fraction tends to zero.

The present audit does not relabel reading that complete reconstruction
as a new independent implementation of its interval engine or as a
fresh proof of every tree-chaos source cited therein. The original
certificate and dependency chain remain the stated imports.

Range superadditivity and the consequent sublinear-block range budget
already appear in
[the September 6 optimizer-law audit](resumed_bound_audit_original_lower_sampling_and_optimizer_laws_2026_09_06.md),
Section 2. The new composition is with the paired-response obstruction
and its Hilbert-vector localization, not a newly named range lemma.

## 6. Finite normalization replay and remaining value gap

After the canonical director file was frozen, this track read all its
Sections 1--7 and reconstructed the sharper rectangular bound as well:
PASS. The discrepancy track's explicit rational corollary was also
independently read and replayed: C=1/2,c=3/10, eventual width lower
.43, K<=2 imply eventual discount delta=2^(-67). Its fixed cap-budget
contradiction margin is exactly 457/409600. The uniform asymptotic
comparison gives existence of the finite threshold, not a practical
numerical n_0.

`computations/paper_bernoulli_2026_09_17_range_localization_audit.py`
passes 49,736 exact principal-profile and partition checks on 1,138
full signings, including EVERY signing through order five. It verifies
beta<=4W and (6) with integer arithmetic, including both one-sided
endpoint inequalities. Python compilation passes. Output:
`tmp/paper_portfolio_2026_09_17/bernoulli/range_localization_audit.json`.
This tests finite normalization; it is not a numerical proof of
Grothendieck's inequality or the asymptotic width lower bound.

The result establishes a fixed response discount below kappa, not a
bound at the stronger slope 3Q(A)/(2n), not a sharp subGaussian proxy,
and not a full-parent upper estimate. Sampling or realizing a favorable
bridge still has to retain all new spins and pay selection/escape.
No change in the certified original extremal interval follows here.
