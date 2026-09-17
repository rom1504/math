# Exact non-Gaussian isotropic laws: quartic feature response versus energy variance

2026-09-17. Localization-track physical-law construction. Pure Walsh
degree four enforces exact second moments, while the absolute-overlap
kernel gives an exact response formula. The construction is elementary;
no Gaussian replacement or approximate covariance repair is needed for
its response. External novelty is not established.

## 1. Two identities used throughout

Let h be uniform on{+-1}^n, n>=5, and let mu_n=E|sum_i h_i|. Set
m=n-1 when n is even and m=n when n is odd. If F is ANY homogeneous
Walsh polynomial of degree four, then

    E F=0,  E F h_i=0,  E F h_i h_j=0,
    E |h dot x| F(h)=-mu_n F(x)/[m(m-2)]             (1)

for EVERY Boolean x. The first line is Walsh orthogonality. The
degree-four Fourier coefficient of the symmetric kernel|sum h_i|
is-mu_n/[m(m-2)], obtained by conditioning on four coordinates and
taking the fourth finite difference; the central-binomial calculation
is fully reconstructed in the
[fixed-power tilt artifact](paper_bernoulli_fixed_power_tilts_2026_09_17.md).
Convolution proves the second line of (1).

Therefore ANY nonnegative density1+lambda F relative to the uniform
cube is exactly normalized, centered, and isotropic. Its response is

    E_nu|h dot x|=mu_n[1-lambda F(x)/(m(m-2))].       (2)

There is also a useful uniform subGaussian bound. If a nonnegative
even density w has Ew=1 and ||w||_2<=L, then its physical law has
linear subGaussian proxy2 max(1,L). Indeed Cauchy--Schwarz and the
coefficientwise Rademacher/Gaussian even-moment inequality give

    E_w(theta dot h)^(2k)
       <=L sqrt((4k-1)!!) ||theta||^(2k)
       <=L 2^k(2k-1)!! ||theta||^(2k).

The second inequality follows by induction: the squared ratio of
the two double factorials grows by
(4k+1)(4k+3)/(2k+1)^2<4. For k>=1 the last expression is at most
[2 max(1,L)]^k(2k-1)!!||theta||^(2k). Evenness removes odd moments;
summing proves the claimed MGF. Exact isotropy alone would not give
this subGaussian conclusion.

## 2. Positive construction for a fixed-rank diffuse feature space

Let U be an n-by-r real matrix with orthonormal columns, P=UU^T,
r>=1. Write u_i for its row vectors and define

    p_i=||u_i||^2=P_ii,  d=max_i||u_i||,
    S(h)=U^T h,  R(h)=||S(h)||^2,
    V(h)=sum_i p_i u_i h_i,  s_2=sum_i p_i^2.

The following polynomial is PURE Walsh degree four:

    F_U(h)=R(h)^2-2(r+2)R(h)+r(r+2)
                              +8 S(h) dot V(h)-6s_2. (3)

For a direct check, take A=P-diag(P) and H_A=(R-r)/2. For any
hollow weighted matrix, the pure degree-four part of H_A^2 is
H_A^2-||Ah||^2+sum_(i<j)A_ij^2. Multiplying by four, using P^2=P,
gives (3). Its coefficient on a four-set{i,j,k,l} is

    8(P_ij P_kl+P_ik P_jl+P_il P_jk).

### Positivity, including an explicit finite leverage bound

Since ||V(h)||<=rd and s_2<=rd^2, put

    E_r(d)=8r sqrt(r+2)d+3(2rd)^(4/3)+6rd^2.

Then for EVERY physical h,

    F_U(h)>=-2(r+2)-E_r(d).                          (4)

To prove this, write R^2-2(r+2)R+r(r+2)
=(R-r-2)^2-2(r+2), and use
sqrt(R)<=sqrt(r+2)+sqrt(|R-r-2|). With v=sqrt(|R-r-2|),
v^4-8rd v>=-3(2rd)^(4/3), proving (4).

Thus for any

    0<=lambda<=[2(r+2)+E_r(d)]^(-1),
    d nu_U/dUniform(h)=1+lambda F_U(h)              (5)

is an ACTUAL centered exactly isotropic sign law. If r is fixed and
d->0, lambda can approach1/[2(r+2)]. A fixed multiplicative slack
in (5) gives a strictly positive density at every sign word.

For CONSTANT-LEVERAGE frames p_i=r/n there is a sharper exact statement:

    F_U(h)=R^2-[2(r+2)-8r/n]R+r(r+2)-6r^2/n,
    F_U(h)>=-2(r+2) for every real R>=0.             (6)

The unrestricted quadratic minimum exceeds-2(r+2) by
2r(r+8)/n-16r^2/n^2>=0, since n>=r. Hence lambda=1/[2(r+2)]
is already valid at finite n in (5), without leverage slack or an
asymptotic correction. Rank-one Boolean directions x_0/sqrt(n) have
constant leverage at EVERY order.

### A dimension-free MGF, not merely covariance control

One has

    E F_U(h)^2<=8r(r+2).                             (7)

An elementary way to see this is to consider a standard Gaussian g
and the fourth Wiener-chaos polynomial

    ||U^Tg||^4-2(r+2)||U^Tg||^2+r(r+2).

Its multivariate Hermite expansion has the same distinct-coordinate
degree-four coefficients as F_U; dropping repeated-coordinate Hermite
terms only decreases its squared L2 norm. In r independent Gaussian
coordinates it equals sum_i H_4(g_i)+2sum_(i<j)H_2(g_i)H_2(g_j).
Orthogonality gives24r+16 binom(r,2)=8r(r+2), proving (7).

For lambda<=1/[2(r+2)], the density in (5) therefore has squared
L2 norm at most1+2r/(r+2)<3. Section 1 shows that ALL these actual
laws have linear subGaussian proxy at most2sqrt(3), independent of
n,r, and the chosen feature space. This is a genuinely non-Gaussian
physical law outside the earlier paired-Gaussian spectral class.
Moreover its information cost relative to the uniform cube obeys the
finite bound

    D(nu_U || Uniform)<=log E_Uniform[(dnu_U/dUniform)^2]<log 3.

Thus the response discount comes from a bounded-information global
rewrite of a physical column, not from independent edge rounding after
a low-information channel. Those are different operations. At the
constant-leverage endpoint, the limiting Gaussian feature density is
the explicit nonnegative function
(||g||^2-r-2)^2/[2(r+2)], with Gaussian mass one and covariance I.

### Exact scalar response and the protected queries

Equations (1)--(3) give the EXACT finite formula

    E_(nu_U)|h dot x|
       =mu_n[1-lambda F_U(x)/(m(m-2))].              (8)

If r is fixed and a=U^T x/sqrt(n), then uniformly over all Boolean x,

    F_U(x)/(m(m-2))=||a||^4+O_r(1/n),
    E_(nu_U)|h dot x|/sqrt(n)
       =kappa[1-lambda||a||^4]+O_r(1/n).             (9)

For clarity, R(x)=n||a||^2<=n, ||V(x)||<=rd<=r, and the lower-degree
terms in (3), divided by n^2, are O_r(1/n). Also m(m-2)=n^2+O(n)
and mu_n/sqrt(n)=kappa+O(1/n). This proof is uniform and uses NO CLT
for the feature coordinates, no nonsingular joint covariance, and no
assumption that query words themselves are random.

In particular all words with

    ||(I-P)x||^2/n<=eta

have response at most kappa[1-lambda(1-eta)^2]+O_r(1/n). For fixed
rank and diffuse rows, the limiting in-subspace response can approach

    kappa[1-1/(2(r+2))].                            (10)

For r=1,2,3,4 this is approximately .664904, .698149, .718096,
.731394, respectively. These values cross the .7404 response target
near the reported upper cap for ranks up to four. This is a real
physical response primitive; it does NOT assert that a complete
nearcode of actual minimizing signings lies near such a fixed-rank
diffuse feature space, or that independent deployment pays parent escape.

## 3. Negative discriminator: actual marked-energy variance can move without cheap response

Now let A itself be an actual full signing, d_A=binom(n,2), and put

    P_4(h)=H_A(h)^2-||Ah||^2+d_A.                   (11)

This too is pure Walsh degree four. If ||A||op<=L sqrt(n), L>=1,
choose0<lambda<=1/(2L^2) and define

    d nu_A/dUniform=1+lambda P_4/n^2.              (12)

It is at least1/2, hence a genuine centered exactly isotropic law with
full physical support. Its degree-four coefficient on a four-set is
2 times the sum of the three matching sign products. That sum is odd,
so its magnitude is either one or three. Thus

    4 binom(n,4)<=E P_4^2<=36 binom(n,4).            (13)

The density L2 norm is at most sqrt(1+(3/2)lambda^2), so Section 1
also gives a dimension-free subGaussian proxy. This construction does
NOT assume such an operator bound for unknown optimizing signings;
bounded-spectrum actual families, including the known-half examples,
are a sufficient concrete domain for the experiment.

The marked-energy variance REALLY changes by a fixed amount:

    E_(nu_A)[H_A(h)/n]^2
       =d_A/n^2+lambda E P_4^2/n^4
       >=d_A/n^2+(lambda/6-o(1)).                  (14)

This follows from orthogonality of P_4 to the constant and degree-two
parts of H_A^2. Yet the exact physical response is

    E_(nu_A)|h dot x|
       =mu_n[1-lambda P_4(x)/(n^2 m(m-2))].          (15)

If additionally Q(A)<=C n^(3/2), then uniformly over EVERY Boolean x,

    E_(nu_A)|h dot x|/sqrt(n)=kappa+O_(C,L)(1/n).    (16)

Indeed |P_4(x)|<=C^2 n^3+(L^2+1/2)n^2. Thus (12) is a precise
non-Gaussian exact-isotropy variance-allocation mechanism that survives
physical sign sampling but fails to lower the leading scalar response,
even at actual ground words. It is not enough merely to alter central
marked-energy fluctuations; the feature construction succeeds because
its quartic polynomial is of order n^2, not merely order n, on the
specific normalized directions one wants to protect.

## 4. Verification and exact scope

The director independently reconstructed the feature polynomial,
level-two cancellation, positivity mechanism, and response formula.
The Bernoulli researcher independently reconstructed all of the
marked-energy construction, including (13)--(16), and returned PASS.
The same researcher subsequently read all of Sections1--2 and4 and
independently reconstructed the feature positivity, Gaussian-chaos L2
bound, dimension-free MGF and uniform exact response: PASS.
The standard pure-degree-four kernel identity is inherited from the
previous fully checked exact Walsh calculation, not newly attributed.

This artifact separates original ingredients (Walsh orthogonality,
the exact radial kernel coefficient, Gaussian moments) from the new
physical constructions and their proved response/cost consequences.
No feature-rank property of actual minimizing nearcodes is assumed.

Finite replay:
`computations/paper_localization_2026_09_17_quartic_isotropic_laws.py`.
It checks12 feature-frame cases (including nonconstant leverage),4
actual full-sign energy cases, and380 exact physical response identities;
all PASS. It also records the finite-n approach to the four feature
response constants. Results are preserved in
`tmp/paper_portfolio_2026_09_17/localization/quartic_isotropic_laws.json`.

## 5. A sharp scope boundary from an isotropic query dual

The Bernoulli researcher supplied the following additional deduction;
the localization track independently reconstructed it. Suppose mu is
an exactly isotropic law on BOOLEAN QUERY words x. Then, for every
frame U of every rank r,

    E_mu R=r,   E_mu S.V=s_2,
    E_mu F_U=E_mu R^2-r(r+2)+2s_2
              <=r(n-r),                           (17)

using0<=R<=n and s_2<=r. Consequently EVERY feature law above with
0<=lambda<=1/[2(r+2)] has dual average response at least

    E_mu E_nuU |h.x|
       >=mu_n[1-r(n-r)/(2(r+2)m(m-2))]
       >=mu_n[1-n/(2m(m-2))].                      (18)

The bound is uniform over ALL ranks, frames, and convex mixtures of
the admissible feature laws. Its normalized leading value is kappa,
not a fixed discount. In particular this entire mixture class cannot
protect an actual ground code carrying an exactly isotropic dual law.
The fixed-rank positive coverage theorem is not contradicted: a query
law with E R=r cannot place all its mass at R approximately n when
r is fixed. This identifies precisely the additional geometric
condition the positive feature primitive needs.
