# Growing-rank cold feature tilts: covariance and uniform physical response

2026-09-17. Independent reconstruction of the director's Fourier-feature
proposal. The rank is genuinely allowed to grow. All trig-ratio estimates
are confined to a good event; its complement is estimated without ratios.

**Subsequent stronger range:** Section7 upgrades covariance error to
O_op(p_*) and O_F(sqrt(n)p_*), and extends physical realization and the
uniform response to r=o(sqrt(n)). Sections1--6 preserve the original
fully elementary r=o(n^(1/3)) reconstruction and its quantitative bounds.

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

The discrepancy researcher observed that the SAME argument is uniform
over every deterministic scalar shift s. Both physical and Gaussian
queries are symmetric, and the exact difference formula is

    E|X-s|-E|sigma_x G-s|
      =(2/pi)int_0^infinity cos(ts)
             [exp(-t^2 sigma_x^2/2)-phi_X(t)]/t^2 dt.

The extra factor has absolute value at most one. Thus the same bound
O_(v,L)((r^2/n)^(1/4)) holds uniformly over ALL Boolean x and ALL real
s. Section7 extends this shifted statement to r=o(sqrt(n)) as well.
This is a scalar all-offset estimate; it does not identify a joint
coupling over queries or justify replacing a maximum of correlated
query fields by independent Gaussian surrogates.

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

## 7. Convex Fourier extension: the stronger range r=o(sqrt(n))

The director proposed the following operator-norm upgrade; this track
reconstructed it independently, as did both the Bernoulli and discrepancy
tracks. Let v,L be fixed as before, but now assume only

    r=o(sqrt(n)),       p_*<=Lr/n.

Then the actual cold law (1) satisfies

    Cov_nu(h)=I+(v-1)(P-diag P)+E,
    ||E||_op<=C_v p_*+negligible,
    ||E||_F<=C_v sqrt(n)p_*+negligible=O_(v,L)(r/sqrt(n)). (13)

The errors are uniform over frames satisfying these assumptions. Its
Boolean-query response still obeys (4), now throughout this larger range.
Its information cost has the stronger absolute approximation

    D(nu_v||Uniform)=(r/2)(v-1-log v)+O_(v,L)(r^2/n)+negligible. (14)

Here every negligible term decays faster than an inverse power of n.
The actual cold law additionally has proxy1+O_(v,L)(r/sqrt(n)) by the
independently audited strong-convex interpolation in
[the cold concentration theorem](paper_discrepancy_cold_feature_subgaussian_2026_09_17.md),
Section8. That MGF theorem is a separate argument, not inferred from (13).

### 7.1 A globally positive auxiliary density

Fix c=1/8 and set q(z)=-log cos(z)-z^2/2 for |z|<=c. Extend q evenly
outside this interval by its quadratic Taylor polynomial at c:

    q(z)=q(c)+q'(c)(|z|-c)+(q''(c)/2)(|z|-c)^2, |z|>=c.

Then q is even, C^2, convex, and, with an absolute constant,

    0<=q''(z)<=C min(z^2,c^2),
    |q'(z)|<=C|z|^3,     |q'(z)|<=C|z|,
    0<=q(z)<=C z^4.

For z_i=sqrt(beta)u_i.g put R(g)=sum_i q(z_i), and define rho by the
globally positive density proportional to

    exp[-||g||^2/(2v)-R(g)].

Its potential V has Hessian at least v^(-1)I and it is even. The ordinary
Brascamp--Lieb variance inequality applies. A primary proof was read in
Carlen--Cordero-Erausquin--Lieb,
[Asymmetric Covariance Estimates of Brascamp--Lieb Type](https://arxiv.org/pdf/1106.0709),
Theorem1.1 at p=2 and the full Section2. The C^2 uniformly convex
potentials here meet its assumptions; smooth approximation preserves
the same bounds. All integrations below have Gaussian decay at infinity.

The rho normalizer relative to N(0,vI) is between1-C_v S_2 and1, since
E_(N(0,vI)) R<=C_v S_2. On the good event G={max_i|z_i|<=c}, its
unnormalized Fourier weight equals the original cosine product exactly.
Outside G, the true partition integrand is bounded by1; the positive
extension integrand is bounded by exp(-beta||g||^2/2)<=1. After the
same normalization as in Section2, their discrepancy is at most

    epsilon_n=C_v n exp(C_v r-c_v/p_*).

It is negligible for r=o(sqrt(n)). The same assertion for covariance
numerators is proved explicitly in 7.3 rather than using trig ratios
outside G.

### 7.2 Latent covariance and cubic residual

Brascamp--Lieb gives Cov_rho(g)<=vI. Write
H=E_rho Hess R. Its defining sum, and E z_i^2<=beta v p_i, imply

    ||H||_op<=C_v p_*,       tr H<=C_v S_2.         (15)

There is also a useful elementary covariance lower bound. If s=grad V,
integration by parts gives E[g s^T]=I and E[s s^T]=E Hess V=v^(-1)I+H.
The block covariance matrix of (g,s) is positive semidefinite, so its
Schur complement yields

    Cov_rho(g)>=(v^(-1)I+H)^(-1).

Consequently

    0<=vI-Cov_rho(g)<=v^2 H,
    ||Cov_rho(g)-vI||_op<=C_v p_*,
    |E_rho||g||^2-vr|<=C_v S_2.                    (16)

Let d_i=q'(z_i) and tau_i=z_i+d_i. On G, tau_i=tan z_i. Both vectors
are centered under the even law rho. For any a in R^n,

    Var_rho(sum_i a_i d_i)
      <=v beta E||U^T diag(q''(z_i))a||^2
      <=C_v sum_i a_i^2 E z_i^4
      <=C_v p_*^2 ||a||^2.                         (17)

The last moment estimate follows directly from Brascamp--Lieb applied
to z_i^2: E z_i^4<=5(beta v p_i)^2. Applying it to z_i^3 similarly
gives E z_i^6<=45(beta v p_i)^3. These estimates use no independence
of the z_i. Equation (17), the covariance Cauchy--Schwarz inequality,
and (16) give, with M=E_rho tau tau^T,

    ||M-beta v P||_op<=C_v p_*,
    M_ii<=C_v p_i.                                (18)

### 7.3 Returning to physical signs, with exact diagonal restoration

The off-diagonal physical Fourier numerator is (6). On G it is exactly
-tau_i tau_j times the positive auxiliary Fourier weight. Outside G,
the true numerator is bounded in absolute value by |z_i z_j|; the
auxiliary numerator is at most C|z_i z_j| because |tau_i|<=C|z_i|.
Cauchy--Schwarz and the Gaussian union tail therefore give, including
the denominator correction, an entrywise discrepancy bounded by

    C_v sqrt(p_i p_j) epsilon_n.

The corresponding matrix has operator and Frobenius norms at most
C_v r epsilon_n. Since physical diagonal moments equal1 exactly,

    Cov_nu(h)=I-M+diag M+T,
    ||T||_op+||T||_F<=C_v r epsilon_n,  T_ii=0.

Combining with (18) proves (13); adding or removing the explicit
(v-1)diag P costs only O_op(p_*) and O_F(sqrt(n)p_*). The identity
Cov_nu(h)_ii=1 is retained, not approximated by a latent diagonal.

### 7.4 Response and entropy in the enlarged range

The partition and characteristic-function proof in Sections2 and4 only
requires S_2=o(1) and negligible Fourier tails. Those now hold because
r=o(sqrt(n)). Equation (13) supplies the uniform second moments needed
for the absolute-value integral. Thus (4), with error O((r^2/n)^(1/4)),
extends unchanged to this larger range.

For the sharper entropy, use the exact physical covariance representation:

    E_nu||U^Th||^2=r-tr(PM)+sum_i p_i M_ii+negligible.

Since U^T tau=sqrt(beta)g+U^Td, (16) gives

    tr(PM)=beta E||g||^2+2sum_i E z_i d_i+E||U^T d||^2
           =beta vr+O_v(S_2).

Indeed |z_i d_i|<=C z_i^4 and
E||U^T d||^2<=sum_i E d_i^2<=C_v sum_i p_i^3<=C_v p_* S_2.
Also sum_i p_i M_ii=O_v(S_2) by (18). Therefore

    E_nu||U^Th||^2=vr+O_v(S_2)+negligible.

Finally log Z=(r/2)log v+O_v(S_2)+negligible, and the exact density (1)
has entropy -beta E||U^Th||^2/2-log Z. This proves (14), including the
absolute o(1) error in the larger rank regime.

The positive auxiliary density and continuous integration by parts are
only proof devices; the final law remains the original physical sign
tilt (1). No boundary reflection process or undefined conditional density
on a trig-sign region is introduced.

### 7.5 Sharper Frobenius repair cost and independent audits

The Bernoulli and discrepancy researchers each independently read all of
Section7 and reconstructed the operator, trace-entropy and tail steps:
PASS. Both also supplied the following sharper Frobenius deduction,
which this track independently verified:

    ||E||_F<=C_v sqrt(r)p_*+negligible
            =O_(v,L)(r^(3/2)/n).                   (19)

Indeed the latent term beta U(Cov g-vI)U^T has rank at most r and
operator norm O(p_*). Equation (17), summed coordinatewise, gives

    tr Cov(d)<=C_v sum_i p_i^3<=C_v r p_*^2,
    ||Cov(d)||_op<=C_v p_*^2.

Therefore ||Cov(d)||_F<=C_v sqrt(r)p_*^2. For the cross term use the
positive block covariance factorization: if A=Cov(z), B=Cov(d), then
Cov(z,d)=A^(1/2)T B^(1/2) for a contraction T, and hence

    ||Cov(z,d)||_F^2<=||A||_op tr B<=C_v r p_*^2.

Summing proves ||M-beta vP||_F<=C_v sqrt(r)p_*. The exact physical
error relative to I+(v-1)(P-diag P) is minus the off-diagonal part of
this matrix, plus the negligible Fourier-tail error. Removing a
diagonal cannot increase Frobenius norm. This proves (19).

The improvement changes the covariance-repair cost, not the present
range of the Fourier argument: the good-region replacement still
requires r=o(sqrt(n)). No extension beyond that range is claimed.
