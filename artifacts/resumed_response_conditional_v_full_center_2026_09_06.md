# A finite-observable full-center certificate: conditional projection onto V

Date: 2026-09-06. Status: exact derivation and rational replay complete,
submitted for independent director audit. This sharpens the same full
response using conditional information; it does not identify the full
distribution of its inverse creation polynomial or change the theorem.

The exact 64-bin lower endpoint is

    .431460392823700540631522256404469656378004801252582084430582.

Subject to director integration, the safely rounded original-problem
statement is liminf M_n/n^(3/2)>=.4314603928237005. The original
convergence question remains open.

## 1. Conditioning retains useful inverse-response information

Use the same refined 21-anchor finite-resolvent mask as the preceding
certificate: alpha=3623/5000, resolvent a=3479/1000, degree D=200, and
the unchanged 21 rational rho coefficients. Let

    V=Ug,  H=1_{|V|<=alpha},  W=UH,
    EV^2=1,  EW^2=p,  EVW=c,
    d=sqrt(p-c^2),  W=cV+dN_1,
    F=sign(W)(1-H).

The full first projection is a_V V+b_N N_1, with exact coefficients
given in `resumed_response_full_center_gain_and_purification_2026_09_06.md`.
Thus

    K=U*F=lambda g+gamma H,
    lambda=a_V-b_N c/d,  gamma=b_N/d,
    tau^2=1-p-a_V^2-b_N^2.

The full theorem gives E H Gamma_tau(K), where
Gamma_tau(k)=E|k+tau N|. Instead of averaging K over the entire mask,
first condition on V:

    E H Gamma_tau(K)
      >= E 1_{|V|<=alpha} Gamma_tau(lambda gbar(V)+gamma),  (1)
    gbar(v)=E[g|V=v].

The convexity is in the shift K; no conditional Gaussian law for g or
K is assumed. The canonical inverse g is a degree-200 polynomial in
its finite Gaussian feedback frame, even though V has countably many
canonical tree coefficients. Hence gbar is an EXACT univariate
degree-at-most-200 even polynomial.

## 2. Exact normalized Hermite coefficients of gbar

Write V=rho.G+sZ, where the 21 anchors and Z are independent standard
Gaussians, R^2=sum rho_T^2, and s^2=1-R^2. Let

    beta_0=p,
    beta_d=-2phi(alpha) He_{d-1}(alpha)/sqrt(d!),  d>=2 even.

For each anchor child monomial h_T, let d_T be its total local degree
and put

    v_T=E[h_T h_{d_T}(V)]
       =sqrt(d_T!/product_j m_j!) product_j rho_j^(m_j),

where m_j are its child multiplicities. Define

    A_d=sum_{T:d_T=d} rho_T v_T,
    B_d=sum_{T:d_T=d} v_T^2,
    L_d=sum_{ell=0}^d binom(d,ell) R^(2(d-ell)) s^(2ell)/(a+ell),
    S_d=sum_{ell=0}^d binom(d,ell) R^(2(d-ell)) s^(2ell)/(a+ell)^2.

The anchor-deleted unnormalized resolvent has exact squared norm

    N=sum_{0<=d<=D even} beta_d^2 [S_d-B_d/a^2].

The d=0 term vanishes because the edge constant is removed. Applying
the innovation resolvent to each multivariate Hermite coefficient and
then taking its inner product with h_d(V) gives

    gbar(v)=sum_{0<=d<=D even} g_d h_d(v),
    g_d=A_d+(s/sqrt(N)) beta_d [L_d-B_d/a].             (2)

This is the multinomial identity for the innovation degree, not an
approximation to a multivariate integral. In particular g_0=rho_edge.
The script independently verifies the normalization-sensitive identity

    sum_{d<=D even} g_d beta_d=E[gH]=c

to an outward interval of width less than 3e-58.

## 3. Exact bin masses and signed first moments

Partition [0,alpha] into 64 equal rational intervals [l_i,r_i], and
combine each with its negative reflection. Let p_i be its Gaussian
mass and k_i the integral of lambda gbar(V)+gamma over that symmetric
bin. They are exactly

    p_i=2[Phi(r_i)-Phi(l_i)],
    k_i=gamma p_i+lambda integral_bin gbar(v)phi(v)dv.

For d>=1 the normalized Hermite identity gives

    integral_l^r h_d(v)phi(v)dv
      =[phi(l)h_{d-1}(l)-phi(r)h_{d-1}(r)]/sqrt(d).

The constant term integrates as g_0 times the mass. Thus every k_i
is evaluated using finite Hermite endpoint recurrences, not quadrature.
Jensen on each bin turns (1) into the rigorous lower bound

    sum_i p_i Gamma_tau(k_i/p_i).                    (3)

The exact replay checks p_i>0 and k_i>0 in all bins. Summing the bin
masses reproduces p; summing the signed first moments reproduces the
old value J=ca_V+db_N. These checks are independent of the Gaussian
absolute-value evaluation in (3).

There is no numerical integration remainder to estimate. The loss
from forgetting within-bin variation is exactly the permitted Jensen
loss, making (3) a lower bound for any fixed partition.

## 4. Exact arithmetic and result

The executable is
`computations/resumed_response_conditional_v_certificate_2026_09_06.py`,
with output
`computations/results/resumed_response_conditional_v_certificate_2026_09_06.json`.
The JSON contains all 201 conditional Hermite coefficient intervals
and all 64 symmetric bin certificates. It uses the inherited exact
Fraction outward-grid arithmetic and finite Gaussian-series enclosures.
No optimizer, floating-point value, or numerical quadrature enters
the proof.

One implementation pitfall was caught by the normalization check:
forming the reciprocal of sqrt(200!) in a fixed absolute 10^-60 grid
would produce a useless interval. The final code instead constructs

    sign(He_{d-1}(alpha))*sqrt(He_{d-1}(alpha)^2/d!)

before interval conversion. The ratio inside the square root is an
exact rational of moderate magnitude, so no tiny factorial reciprocal
is rounded away. This is a numerical-enclosure repair, not a change
to the mathematical coefficient formula.

The exact output is

    conditional-V 64-bin lower bound in
    [.431460392823700540631522256404469656378004801252582084430582,
     .431460392823700540632427984625288259518371795940743872444826].

Its lower endpoint exceeds 8629/20000=.43145, checked exactly. The
creation derivative energy is still less than
.999904932505320245853669704564681545270483582860128726229320<1.

For comparison only, a noncertifying floating conditional integral
gives approximately .431460824655. The exact finite-bin certificate
deliberately stops short of that value. No claim that the full
functional depends only on V, or that conditioning on V is optimal,
is made; (1) is a rigorously useful information reduction.
