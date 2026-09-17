# Full-parent cost of the operator-free cheap-column law

2026-09-17. Localization-track composition of the low-cap uniform-response
theorem with cloned all-energy profiles. Every child and bridge below is
an ACTUAL full signing. The result is a complete quantitative upper
certificate, followed by two precise limitations; it is not a favorable
cross-order recurrence.

## 1. The response law has a global energy envelope

Let W be a full signing of order n with Q(W)<=C n^(3/2), and write
C_W=Q(W)/n^(3/2). Assume 0<c<C and the
[low-cap response theorem](paper_localization_low_cap_response_quantitative_2026_09_17.md)
has produced one centered isotropic physical column law nu such that

    E_nu exp(theta dot h)<=exp(3||theta||^2/4),
    m(x):=E_nu|h dot x|<= (kappa-delta)sqrt(n)
                      when |H_W(x)|>=c n^(3/2).     (1)

The law lies in the convex hull of equal bounded-spectrum Gaussian-sign
pairs. Their uniform scalar comparison also gives, on EVERY Boolean x,

    m(x)<= (kappa+e_n)sqrt(n),
    e_n=O(n^(-1/6)sqrt(log(en))).

Put alpha=delta/(C-c), and e(x)=|H_W(x)|/n^(3/2). Then the director's
elementary envelope holds over the ENTIRE cube:

    m(x)/sqrt(n)<=kappa+e_n-alpha[e(x)-c].           (2)

For e(x)<=c the right side is at least the global upper bound. For
c<=e(x)<=C it is at least kappa+e_n-delta, hence at least the
high-energy upper bound. This is a deliberately weak affine envelope,
not an interpolation theorem for the exact response function.

If d_x=Q(W)-|H_W(x)| and q columns will be added, (2) implies exactly

    |H_W(x)|+q m(x)
       <=Q(W)+q sqrt(n)mu_W-(1-alpha q/n)d_x,
    mu_W=kappa+e_n-alpha(C_W-c).                    (3)

This removes the need for a separate low-energy MEAN argument. It does
not remove the centered random supremum, which is paid next.

## 2. A completely unconditional full-spin iid-parent certificate

Take q iid columns h_1,...,h_q from nu, and ANY actual full signing D
of order q. They form the bridge of an actual parent P of order n+q.
For every old x and new y,

    |H_P(x,y)|<=|H_W(x)|+sum_j|h_j dot x|+Q(D).

Set Z_x=sum_j[|h_j dot x|-m(x)]. For a single column and independent
copy V'=h' dot x, Jensen and reverse triangle give

    E exp(lambda(|V|-E|V|))
       <=E cosh(lambda(|V|-|V'|))
       <=E cosh(lambda(V-V'))
       <=exp((3/2)lambda^2 n).

Thus Z_x has variance proxy3qn. If alpha q/n<=1, (3) and a union
over the full2^n old cube prove, with probability at least1-exp(-u),

    Q(P)<=Q(D)+Q(W)+mu_W q sqrt(n)
                       +sqrt[6qn(n log2+u)].        (4)

Every new-spin maximum and every old word has been paid. In particular
there is no hidden assumption that a maximizing old word remains in
the protected energy code. The last term has normalized size
O(sqrt(q/n)), not o(q/n), at a macroscopic added density.

## 3. Common-column absolute increments really are subGaussian

The law of h can have dependent coordinates. This does not obstruct the
following increment estimate. For w=x-y, put
D(h)=|h dot x|-|h dot y| and let D' use an independent column. Then

    E exp(lambda(D-ED))<=E cosh(lambda(D-D'))
       <=E cosh(2lambda h dot w)
       <=exp(3lambda^2||w||^2).                    (5)

The middle inequality follows from
|D-D'|<=|h dot w|+|h' dot w| and
cosh(a+b)<=[cosh(2a)+cosh(2b)]/2. It does not require D itself to
be symmetric. Equation (5) gives the process Z_x subGaussian increments
with variance proxy6q||x-y||^2, even though all its queries use the
same q physical columns. There is no coordinate-independence claim.

Dirksen's primary Theorem3.2 and Gaussian majorizing measures therefore
give, on any antipodal code F containing a chosen anchor x_0,

    sup_(x in F)|Z_x-Z_(x_0)|
       <=C sqrt(q)w_G(F)+C sqrt(qn u)               (6)

with probability at least1-exp(-u), up to universal constant changes.
The fixed anchor separately has variance proxy3qn. The primary chaining
proof and Gaussian-comparison input were reconstructed in the campaign;
no unsupported Gaussian replacement of h is used in (5)--(6).

## 4. Clone preparation supplies a sharper, still fully paid parent bound

Start with any actual signing A of order n with Q(A)<=C_in n^(3/2),
where C_in<C is fixed. Choose eta=o(1), eta>=n^(-2/3), and prepare
W using cloned blocks, with

    p=floor(eta^(-1/2)),  r=floor(sqrt(n)eta^(3/4)),
    B=n/sqrt(p),  ell=log(e/eta).

For all sufficiently large n, the same actual W obeys

    Q(W)<=Q(A)+C_0 sqrt(eta)n^(3/2)<=C n^(3/2),
    w_G(E_W(D))<=C_0 sqrt(ell)(D/r+B)
                                      for ALL D>=0. (7)

Apply the operator-free response theorem to W itself; the original A
need not satisfy an operator bound. Sample iid bridge columns from its
law (1). Suppose

    1<=q<=c_0 r^2/ell,
    alpha q/n<=1/2,
    T_q=C_1 sqrt(q ell)B.                           (8)

There are universal choices of sufficiently small c_0 and large C_1
such that, with probability at least

    1-C exp[-c ell B^2/n],

the resulting ACTUAL full parent satisfies

    Q(P)<=Q(D)+Q(W)+mu_W q sqrt(n)+C T_q.            (9)

The constants in this paragraph do not denote the response constants
C_1,C_2 in the separate quantitative-law artifact. Here they are
universal chaining/profile constants. If the integer interval in (8)
is empty, the assertion has no q to which it applies.

To prove (9), (6)--(7) give an anchored centered level supremum bounded
by C sqrt(q ell)(D/r+B), with upper deviations C sqrt(qn u).
Condition (8) makes its D coefficient at most1/32. Choose the large
constant in T_q to absorb the intercept. A dyadic shell argument over
EVERY d_x>=T_q then ensures

    |Z_x-Z_(x_0)|<=d_x/4 for all d_x>=T_q,

except with probability C exp[-cT_q^2/(qn)]. On the inner level the
centered supremum and scalar anchor are both O(T_q), with the same
failure order. Equation (3) retains at least d_x/2 as deterministic
loss, so all outside levels are paid as well. This proves (9) without
a separate high/low mean partition or unexamined wrong-polarity branch.

At the largest allowed q, one has

    q/n=O(eta^(3/2)/ell),
    T_q=O(eta n^(3/2)),
    failure exponent=Omega(n sqrt(eta)ell).

The bound includes the unaltered cap Q(D) of all fresh internal edges.
No upper control of D has been inferred from its being a full signing.

## 5. Why this does not yet yield a favorable extension

First, write epsilon=q/n. The compatibility condition in (8) requires
eta^(3/2)>=C epsilon ell. The certified preparation cost in (7) is
of order sqrt(eta)n^(3/2), hence no better than the scale
(epsilon ell)^(1/3)n^(3/2) in this combined certificate. This is much
larger than epsilon n^(3/2) as epsilon->0. This is a limitation of
the paid bounds, not a lower bound on actual preparation cost.

Second, the particular latent band itself has a sharp mean-response
floor, independently identified by the discrepancy researcher. If
R_+-I=T, then ||T||op<=1/2 and the entrywise arcsine power series plus
Schur tensor compression imply

    ||kappa^2 asin(T)||op<=kappa^2 asin(1/2)=1/3.

Thus every query has |v_x|<=1/3 in (4) of the response artifact, and
every paired law, and every convex mixture of them, has

    E|h dot x|/sqrt(n)>=kappa f(1/3)-e_n
          =0.786393873897...-e_n.                   (10)

This exceeds3/4, and hence exceeds the first-order allowance1.5c_0
for every normalized child cap c_0<=1/2. The floor is asymptotically
attained by taking T equal to one half a perfect-matching permutation
matrix and querying the all-ones word: each pole is a sum of independent
two-coordinate sign pairs, so the elementary scalar CLT gives (10).

Consequently even perfect elimination of the stochastic selection term
would not make this Gaussian-pair class alone certify the desired
linear parent slope. This is stronger than saying that the explicit
discount happened to be numerically small. It is restricted to the
latent spectral interval[1/2,3/2], not all physical isotropic laws.

The operator-free response theorem remains substantive: it gives one
uniformly subGaussian, exactly isotropic law with a fixed discount on
a full macroscopic energy code of every low-cap actual signing. Equations
(4) and (9) explain precisely what its direct iid deployment can certify
and what extra geometry or a different response class must still supply.

Independent audit: the Bernoulli researcher read Sections 1--5 in full
and returned PASS, including the common-column increment normalization,
global affine mean envelope, every-energy shell absorption, failure
exponent, full fresh-child cost, preparation limitation, and class floor.

## 6. An actual wrong-slope theorem, even for an adaptively chosen bounded-cap child

The class floor yields more than a limitation of the upper certificate.
Let beta_*=kappa f(1/3), and independently sample the q bridge columns
from ANY members of the half-band paired class; the laws may differ.
Fix one old signed ground word x_0 with sign sigma_0. For EVERY
realized bridge, EVERY actual new child D, and regardless of whether
D was chosen before or after observing that bridge,

    Q(P)>=Q(W)+sum_j|h_j dot x_0|-Q(D).              (11)

Indeed set y_j=sigma_0 sign(h_j dot x_0), completing zero ties
arbitrarily. Multiplying the parent energy at(x_0,y) by sigma_0 makes
the old energy Q(W) and every bridge term positive; the internal fresh
energy is at least-Q(D). This evaluates one literal Boolean parent
word and is a lower bound, not a mean-field approximation.

Each column mean at x_0 is at least(beta_*-e_n)sqrt(n) by (10).
Its centered absolute response has variance proxy3n, as proved in
Section 2. Therefore, with probability at least1-exp(-q t^2/6),

    Q(P)>=Q(W)+(beta_*-e_n-t)q sqrt(n)-Lq^(3/2)     (12)

SIMULTANEOUSLY for every possible adaptive child D satisfying
Q(D)<=Lq^(3/2). No independence of D from the bridge is assumed.
The restriction on its normalized cap is essential: an arbitrary
adaptive child with cap of order q^2 is not covered by this bound.

For an explicit wrong-slope conclusion, suppose Q(W)<=C n^(3/2)
with gamma=beta_*-3C/2>0. Put t=gamma/4 and take n large enough that
e_n<=gamma/4. If epsilon=q/n<=1 obeys

    L sqrt(epsilon)+(3C/8)epsilon<=gamma/4,

then (12) and(1+epsilon)^(3/2)<=1+3epsilon/2+3epsilon^2/8 imply

    Q(P)>= [Q(W)/n^(3/2)](n+q)^(3/2)
                           +(gamma/4)q sqrt(n)      (13)

with probability at least1-exp(-q gamma^2/96), uniformly over those
adaptive children. Thus the normalized parent cap strictly increases
by at least(gamma/4)epsilon/(1+epsilon)^(3/2). This applies throughout
C<=1/2 and every fixed L, at all sufficiently small fixed added density.

The discrepancy track's
[class-floor artifact](paper_discrepancy_paired_spectral_band_response_floor_2026_09_17.md)
has an even cleaner independent-child statement: the aligned fresh
spins are then iid fair signs, so their child energy is negligible
on the n^(3/2) scale without any cap assumption on that independent
child. Section 6 supplies the complementary adaptive-child version
by explicitly paying its cap. Both leave rare globally selected
bridge realizations and different physical-law classes outside scope.
