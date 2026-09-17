# Growing-rank cold feature tilts: covariance and uniform physical response

2026-09-17. Independent reconstruction of the director's Fourier-feature
proposal. The rank is genuinely allowed to grow. All trig-ratio estimates
are confined to a good event; its complement is estimated without ratios.

## 1. The physical law and quantitative theorem

Fix 0<v<1 and L<infinity, and put beta=1/v-1. Let U be n-by-r with
orthonormal columns, P=UU^T, and let u_i denote its rows. Assume

    p_i=||u_i||^2,   p_*=max_i p_i<=Lr/n,
    S_2=sum_i p_i^2,   r=o(n^(1/3)).

Define the ACTUAL sign law

    dnu_v/dUniform(h)=exp[-beta||U^T h||^2/2]/Z,
    Z=E_Uniform exp[-beta||U^T h||^2/2].             (1)

It is centered by global sign symmetry. There are constants depending
only on v such that, with

    tau=2n exp(C_v r-c_v/p_*),

the following conclusions hold for all sufficiently large n:

    |Z/v^(r/2)-1|<=C_v(S_2+tau),                   (2)

    E_nu hh^T=I+(v-1)(P-diag P)+E,
    E_ii=0,
    |E_ij|<=C_v sqrt(p_i p_j)(S_2+p_i+p_j+tau),
    ||E||_F<=C_v r(S_2+p_*+tau)=O_(v,L)(r^3/n).   (3)

In particular the Frobenius error tends to zero. Since r=o(n^(1/3)),
tau decays faster than any inverse power of n. Define
a_x=U^T x/sqrt(n). Uniformly over EVERY Boolean x,

    E_nu |h.x|/sqrt(n)
      =kappa sqrt[1+(v-1)||a_x||^2]
                              +O_(v,L)((r^2/n)^(1/4)). (4)

The asymptotic response conclusion also follows without a claimed rate
from the characteristic functions and uniform second moments below.
The stated rate uses their elementary absolute-value Fourier integral.
No random-query assumption or nondegenerate joint feature Gram is used.

This theorem does NOT claim that (1) is exactly isotropic, nor infer a
dimension-free MGF for arbitrary real vectors from (3). It provides the
cold branch that can be combined with a separately proved hot branch
and an explicitly paid physical covariance correction. A subsequent
independent argument DOES establish cold subGaussian proxy1+o(1) in
this same rank range: see the
[separate cold concentration proof](paper_discrepancy_cold_feature_subgaussian_2026_09_17.md),
independently reconstructed by this track. That argument uses external-
field Fourier determinants and evenness, not covariance alone.

## 2. Partition function, including bad Fourier regions

Let g be standard Gaussian in R^r and z_i=sqrt(beta) u_i.g. Fourier
inversion of a Gaussian, followed by independence of the physical signs,
gives the exact identity

    Z=E_g prod_i cos(z_i).                          (5)

Fix a small numerical c_0>0, for example1/8, and put
G={max_i|z_i|<=c_0}. On G all cosines are positive and

    prod_i cos(z_i)=exp[-beta||g||^2/2-R(g)],
    0<=R(g)<=C sum_i z_i^4.

Thus the difference from exp[-beta||g||^2/2] on G is at most that
Gaussian weight times C sum_i z_i^4. The reference integral is
Z_0=v^(r/2), and dividing by Z_0 changes g's law to N(0,vI). Hence
the integrated good-event error relative to Z_0 is at most C_v S_2.

On G^c use |prod cos|<=1, not a ratio or a logarithm. A Gaussian
union bound gives

    P(G^c)<=2n exp[-c_0^2/(2 beta p_*)].

The reference integrand is also at most1. Dividing these two tail
errors by v^(r/2) and enlarging constants gives at most C_v tau.
This proves (2), including possible negative cosine products outside G.
The denominator is therefore at least Z_0/2 eventually, despite Z_0
itself being exponentially small in the growing rank.

## 3. Exact off-diagonal numerator and Frobenius error

For i!=j the exact Fourier numerator is

    E_Uniform h_i h_j exp[-beta||U^Th||^2/2]
      =-E_g sin(z_i)sin(z_j)prod_(k!=i,j)cos(z_k).  (6)

On G, and ONLY there, rewrite this as
-tan(z_i)tan(z_j)prod_k cos(z_k). Since

    |tan(z_i)tan(z_j)-z_i z_j|
        <=C|z_i z_j|(z_i^2+z_j^2),

the difference from -z_i z_j exp[-beta||g||^2/2] is bounded by

    C exp[-beta||g||^2/2]|z_i z_j|
                          [z_i^2+z_j^2+sum_k z_k^4].

Gaussian moments under N(0,vI) bound its integral divided by Z_0 by

    C_v sqrt(p_i p_j)(p_i+p_j+S_2).                (7)

For example E|z_i z_j|z_k^4 is at most
sqrt(E z_i^2 z_j^2) sqrt(E z_k^8)
<=C_v sqrt(p_i p_j)p_k^2. Correlations among these coordinates are
retained in this Cauchy--Schwarz estimate.

On G^c, use |sin(z_i)sin(z_j)prod_others cos|<=|z_i z_j|. The
same bound applies to the absolute reference numerator. Another
Cauchy--Schwarz inequality bounds the normalized tail error by

    C_v sqrt(p_i p_j) P(G^c)^(1/2)/Z_0
       <=C_v sqrt(p_i p_j) tau,

after reducing c_v in tau. This also handles p_i=0 exactly; no division
by p_i occurs. The reference numerator divided by Z_0 is
-beta v P_ij=(v-1)P_ij. Normalizing by Z rather than Z_0 adds at most
C_v sqrt(p_i p_j)(S_2+tau), proving the entry estimate in (3).
Diagonal correlations are exactly1 under every sign law, which explains
the subtraction of diag P in (3).

For completeness, the Frobenius bounds are

    ||(sqrt(p_i p_j))_ij||_F=r,
    ||(sqrt(p_i p_j)(p_i+p_j))_ij||_F<=2r p_*.

Thus ||E||_F<=C_v r(S_2+p_*+tau). Since S_2<=rp_*<=Lr^2/n,
the advertised O(r^3/n) error follows. This deliberately conservative
Frobenius estimate is enough for rank o(n^(1/3)); no unproved cancellation
between its entries is used.

## 4. All-query characteristic functions and absolute response

For a Boolean x and real t, the exact characteristic numerator is

    E_Uniform exp[it h.x/sqrt(n)]
                    exp[-beta||U^Th||^2/2]
       =E_g prod_i cos(z_i+t x_i/sqrt(n)).          (8)

For |t|<=c_0 sqrt(n), the event G places every shifted argument inside
[-2c_0,2c_0]. The Gaussian reference integral is exactly

    E_g exp[-sum_i(z_i+t x_i/sqrt(n))^2/2]
       =Z_0 exp[-t^2 sigma_x^2/2],
    sigma_x^2=1+(v-1)||a_x||^2 in[v,1].            (9)

Relative to that integral, g has Gaussian law
N(-sqrt(beta)v t a_x,vI). For the shifted arguments, their variances
are beta v p_i and their means are
t[x_i/sqrt(n)-beta v u_i.a_x]. The fourth-moment sum is bounded by

    C_v[S_2+t^4(p_*+1/n)]
       <=C_v(1+t^4)(S_2+p_*+1/n),                 (10)

because sum_i(u_i.a_x)^4<=p_*||a_x||^2<=p_*. On the bad event G^c,
both the cosine product and the Gaussian reference integrand have
absolute value at most1; their normalized tails are again bounded by
C_v tau. Thus, also paying the partition denominator error,

    sup_x |phi_(nu,x)(t)-exp[-t^2 sigma_x^2/2]|
       <=C_v(1+t^4) delta_n,
    delta_n=S_2+p_*+1/n+tau.                       (11)

This estimate is uniform over x, even if ||a_x||=1. In addition (3)
implies a uniform second-moment bound for h.x/sqrt(n): its difference
from sigma_x^2 is at most (1-v)r/n+||E||_op, which is bounded and tends
to zero under the stated growing-rank hypothesis.

Finally use the exact symmetric absolute-value identity

    E|X|=(2/pi)int_0^infinity [1-E cos(tX)]/t^2 dt.

The difference of the integrals below epsilon is O(epsilon) by the
uniform second moments. Above R it is O(1/R). On [epsilon,R], (11)
gives O(delta_n/epsilon+delta_n R^3). Choosing
epsilon=delta_n^(1/4) and R=delta_n^(-1/4) gives O(delta_n^(1/4)).
Here R<=n^(1/4), since delta_n>=1/n, so (11)'s required t-range is
valid eventually. Since delta_n=O_(v,L)(r^2/n), this proves (4).

## 5. Provenance and limitation

The director proposed the growing-rank Fourier route and its target
Frobenius scale. The localization track independently reconstructed
the complete argument above, in particular the tails in (6), the
normalization by exponentially small Z_0, and uniform shifted-query
control. The only imported identities are elementary Gaussian Fourier
integration and the absolute-value characteristic-function identity.
The resulting law is physical at every n; no approximate-sign output
or weighted-edge relaxation is used. External novelty is not asserted.

The discrepancy researcher independently read Sections1--5 in full and
reconstructed the tail, covariance, shifted fourth-moment, and quantitative
absolute-response steps: PASS. Finite replay is
`computations/paper_localization_2026_09_17_growing_rank_cold_tilt.py`;
18 physical-cube cases and156 Gaussian-quadrature Fourier identities all
PASS. It explicitly avoids trig ratios in the covariance numerator.
Results: `tmp/paper_portfolio_2026_09_17/localization/growing_rank_cold_tilt.json`.

## 6. Information cost at growing rank

There is a quantitative entropy consequence not requiring differentiation
of a trig approximation. The partition estimate is uniform as beta ranges
over a fixed compact neighborhood of its positive value. Put

    f_n(beta)=log Z(beta),
    f_0(beta)=-(r/2)log(1+beta).

By (2), |f_n-f_0|<=C_(v,L)r^2/n throughout that neighborhood; the
exponential tail is negligible. Both functions are smooth and f_n is
convex, since it is a log moment-generating function. Convex difference
quotients on either side, and |f_0''|<=C_v r, give

    |f_n'(beta)-f_0'(beta)|
       <=C_(v,L)[r h+(r^2/n)/h].

Taking h=sqrt(r/n) proves

    |E_nu ||U^Th||^2-vr|<=C_(v,L)r^(3/2)/sqrt(n)=o(1).

Since the physical density is known exactly,

    D(nu_v||Uniform)=beta f_n'(beta)-f_n(beta)
       =(r/2)(v-1-log v)+O_(v,L)(r^(3/2)/sqrt(n)).  (12)

Thus the Gaussian variance-change information cost is recovered with
an absolute o(1) error throughout the claimed growing-rank range.
This statement is for the cold law before any covariance correction.
A later correction of mass o(1) does not by itself preserve absolute
o(1) entropy accuracy when the main cost grows like r; that needs a
separate error estimate.

The same partition estimate also gives, for fixed t sufficiently small,

    log E_nu exp[t||U^Th||^2]
       =-(r/2)log(1-2vt)+O_(v,L)(r^2/n),

by the exact ratio Z(beta-2t)/Z(beta). In particular if r grows, the
feature energy divided by r concentrates around v by ordinary Chernoff
bounds. This concerns the feature energy, not an unproved arbitrary-
linear subGaussian bound for the physical sign vector.
