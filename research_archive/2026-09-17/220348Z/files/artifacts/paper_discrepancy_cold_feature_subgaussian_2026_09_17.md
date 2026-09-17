# Cold growing-rank feature laws are asymptotically 1-subGaussian

2026-09-17. Proved directly for all real test directions. Independent
director, Bernoulli, and localization reconstructions PASS. The proof
uses an absolute Fourier integral, a Gaussian determinant comparison,
and a symmetry lemma that removes an additive MGF error. It does NOT
infer concentration from covariance convergence or import a rank-one
Ising result as a growing-rank theorem.

**Stronger rank range:** Section 8 upgrades the result to r=o(sqrt(n)),
with proxy 1+O_(v,L)(r/sqrt(n)), by strongly convex Gaussian interpolation.
Sections 1--7 retain the simpler initial r=o(n^(1/3)) proof and its audit.

## 1. Statement and precise hypotheses

Fix 0<v<1 and beta=1/v-1. Let U be n-by-r with U^T U=I_r, let
P=UU^T, and let p_i=P_ii, p_*=max_i p_i. Consider the actual physical law

    dnu/dUniform(h)=exp[-beta||U^T h||^2/2]/Z,
    Z=E_Uniform exp[-beta||U^T h||^2/2].                    (1)

Suppose p_*<=Lr/n for a fixed L, and r=o(n^(1/3)). Then nu is centered,
has full support, and is K_n-subGaussian in EVERY real direction, where

    E_nu exp(t.h)<=exp(K_n||t||^2/2) for all t in R^n,
    K_n=1+O_(v,L)(sqrt(r/n)(r+log n))=1+o(1).             (2)

This includes growing rank and arbitrary diffuse real features. It does
not require sign-valued columns, block structure, or a finite query set.
The rank-zero case is the independent cube and is immediate.

The only previously proved input is the relative partition lower bound
from `paper_localization_growing_rank_cold_tilt_2026_09_17.md`, Section 2:

    Z >= v^(r/2)(1-delta_n),
    0<=delta_n<=C_v[S2+tau_n]=o(1),
    S2=sum_i p_i^2,
    tau_n=2n exp(C_v r-c_v/p_*).                           (3)

That bound is a complete elementary cosine-product reconstruction,
including the tail after the exponentially small normalization. It is
uniform under the hypotheses above. Taking a nonnegative upper error
delta_n in (3) avoids any assertion about the sign of the actual error.

## 2. An arbitrary external test vector and a Fourier modulus

Fix ANY t in R^n, with no bound on its norm or coordinates. Let
m_i=tanh(t_i), D_i=1-m_i^2=sech^2(t_i), and D=diag(D_i). Tilting the
independent cube by exp(t.h) gives independent signs with means m_i.
The exact Gaussian Fourier identity therefore gives

    E_nu exp(t.h)
      =[product_i cosh(t_i)] N_t/Z,
    N_t=E_g product_i[cos(z_i)+i m_i sin(z_i)],
    z_i=sqrt(beta) u_i.g,   g~N(0,I_r).                    (4)

Although the integral in (4) is written in complex form, N_t is the
positive real expectation of exp[-beta||U^T h||^2/2] under that biased
product law. Its absolute value is bounded using

    |cos(z_i)+i m_i sin(z_i)|
       =sqrt(1-D_i sin^2(z_i))<=1.                         (5)

For |z|<=1, sin^2(z)>=z^2-z^4/3. Since log(1-u)<=-u,

    log sqrt(1-D_i sin^2 z_i)
       <=-D_i z_i^2/2+D_i z_i^4/6.                       (6)

No division by a cosine and no logarithm of a signed cosine is used.
The estimate remains valid at zeros by continuity.

## 3. A uniform Gaussian ball, with the full normalization paid

For n>=2 put

    R_n^2=[2 log(2/v)+4]r+16 log n,
    B_n=beta^2 p_* R_n^4/6.                               (7)

On the event ||g||^2<=R_n^2, every z_i obeys
z_i^2<=beta p_* R_n^2=o(1), so (6) applies for all sufficiently large n.
Also, using the projection identity sum_i(u_i.g)^2=||g||^2,

    sum_i z_i^4<=beta^2 p_* ||g||^4<=6 B_n.

Consequently the modulus of the product in (4), on this event, is at
most exp(B_n) exp[-beta g^T U^T D U g/2]. The reference Gaussian integral
over the whole space is

    A_t=det(I_r+beta U^T D U)^(-1/2)>=v^(r/2),             (8)

because 0<=U^T D U<=I_r. On the complementary event use (5). Gaussian
Chernoff at parameter 1/4 gives the explicit estimate

    P(||g||^2>R_n^2)
       <=2^(r/2) exp(-R_n^2/4)
       =v^(r/2) exp(-r)n^(-4)
       <=v^(r/2)n^(-4).                                 (9)

Combining (8)--(9) yields the MULTIPLICATIVE upper estimate

    N_t <=[exp(B_n)+n^(-4)] A_t,                         (10)

uniformly over t. The factor v^(r/2) in (9) is essential; an unnormalized
small Gaussian tail alone would not suffice in growing rank.

## 4. Determinant comparison gives a defective MGF bound

Let H=U^T D U. Its eigenvalues lie in [0,1], so

    log[A_t/v^(r/2)]
      =(1/2){r log(1+beta)-log det(I_r+beta H)}
      <=(beta/2) Tr(I_r-H)
      =(beta/2) sum_i p_i tanh^2(t_i)
      <=(beta p_*/2)||t||^2.                             (11)

The scalar inequality used here is
log(1+beta)-log(1+beta s)<=beta(1-s), 0<=s<=1.
Using product cosh(t_i)<=exp(||t||^2/2), (3), and (10), for delta_n<=1/2
we obtain the exact convenient upper bound

    log E_nu exp(t.h)<=a_n||t||^2+b_n,
    a_n=(1+beta p_*)/2,
    b_n=B_n+2delta_n+n^(-4)>0.                            (12)

Under the stated leverage and rank hypotheses,

    b_n=O_(v,L)(r(r+log n)^2/n)=o(1).                     (13)

The harmless n^(-4) can be included in this bound for r>=1. The
relative partition and tail errors in (3) are smaller on this scale.
Equation (12) is NOT yet the desired subGaussian inequality, because
its additive error cannot simply be discarded at small t.

## 5. Symmetry removes the additive defect

Here is a general elementary lemma. Suppose a globally symmetric random
vector X has

    log E exp(t.X)<=a||t||^2+b for all t, with b>0.

Then for all t,

    log E exp(t.X)<=c(a,b)||t||^2,
    c(a,b)=[exp(a sqrt(b)+b)-1]/sqrt(b).                   (14)

For ||t||>=b^(1/4), the original bound gives coefficient a+sqrt(b),
which is at most c(a,b). For 0<||t||<b^(1/4), take u in the same direction
with ||u||=b^(1/4), and lambda=||t||/||u||. The even power series gives

    cosh(lambda z)-1<=lambda^2[cosh(z)-1],  0<=lambda<=1.

Symmetry, the assumed estimate at u, and log(1+s)<=s therefore give

    log E exp(t.X)
      <=[||t||^2/sqrt(b)] [exp(a sqrt(b)+b)-1],

which proves (14). The case t=0 is exact.

Apply (14) to (12). The physical cold law is globally symmetric, and
a_n stays bounded. Its subGaussian proxy is therefore

    K_n=2c(a_n,b_n)
       =1+beta p_*+O_v(sqrt(b_n))=1+o(1),                 (15)

with the rate in (2). This proof is uniform in every real t, including
large, sparse, dense, or n-dependent directions. No external-field
covariance theorem has been assumed.

## 6. Consequences and exact scope

Combining this with the independently proved hot-component proxy v for
v>1, a fixed finite cold/hot mixture has a common dimension-free proxy
max{1+o(1),max_j v_j}. The exact Gaussian-sign covariance repair used
in `paper_director_growing_rank_variance_realization_2026_09_17.md` has
proxy 3/2. Consequently its exactly isotropic physical output has the
common finite proxy max{3/2,max_j v_j,1+o(1)}. All-query response and
information bounds are unchanged.

The same proof covers the adaptive covariance-capture rank budget:
actual rank q<=r and p_*<=r/(eta n) give the same o(1) bounds, with fixed
eta. Its compact minimax class can therefore impose this common MGF
bound as well as exact covariance and its information budget. Mixing
and a final independent-cube component preserve the bound.

This closes a tail-control obligation for the growing-rank realization;
it does not prove that an actual minimizing nearcode has low-rank
covariance capture. Nor does it make the common proxy uniform when the
desired response epsilon tends to zero: the hot variance and the sharp
scalar subGaussian price must still be paid.

## 7. Primary-literature boundary

Mikulincer--Sohn, [Fast mixing in Ising models with a negative spectral
outlier via Gaussian approximation](https://arxiv.org/abs/2512.22803),
v2 (2026), Theorem 9, proves a stronger all-external-field covariance
statement for one diffuse negative mode. That statement does not by
itself cover this growing-rank negative projection, and it was not used
in the proof above. The direct MGF argument avoids that extension.
External novelty of this elementary combination remains unestablished.

Finite replay is in
`computations/paper_discrepancy_2026_09_17_cold_feature_sg.py`.
It checks 1,800 determinant inequalities, finite-cube covariance/MGF
calculations on three frames, and two rank-one binomial sequences through
n=512. These are floating diagnostics, not the proof of the uniform
asymptotic bound. In particular small non-diffuse frames can have
covariance norm substantially above one; (2) is not an exact finite-n
1-subGaussian assertion.

## 8. Strongly convex interpolation extends the range to r=o(sqrt(n))

The director proposed replacing the crude latent-ball quartic bound by
a globally controlled interpolation. The following reconstruction closes
that step. Under the SAME fixed-v leverage assumption p_*<=Lr/n, it gives

    r=o(sqrt(n))  =>  K_n=1+O_(v,L)(r/sqrt(n))=1+o(1).    (16)

This is a tail-control statement for the actual cold law; any improvement
to the covariance-repair theorem must be justified separately.

### 8.1 A convex truncated quartic and its exact derivative bounds

Fix 0<c<=min{1/2,1/(2sqrt(beta))}. Define the even function

    f_c(z)=z^4,                                      |z|<=c,
    f_c(z)=c^4+4c^3(|z|-c)+6c^2(|z|-c)^2,            |z|>=c.

The values and first two derivatives match at |z|=c. Directly,

    f_c is C^2, even and convex,
    0<=f_c(z)<=z^4,       0<=f_c''(z)<=12c^2.             (17)

The inequality f_c<=z^4 outside the central interval follows from the
exact remainder 4c(|z|-c)^3+(|z|-c)^4. Put

    F(g)=(1/6)sum_i f_c(sqrt(beta)u_i.g).

The orthonormal-column identity gives

    0<=Hess F(g)<=2beta c^2 I_r<=I_r/2.                  (18)

Also F(0)=0 and grad F(0)=0, so F(g)<=||g||^2/4 globally. This supplies
integrable domination in every interpolation and differentiation below.

### 8.2 Uniform exponential-quartic bound under every reference covariance

Let Q be ANY symmetric matrix with Q>=I_r, and write gamma_Q for the
centered Gaussian law of precision Q. For s in[0,1] let

    drho_s/dgamma_Q=exp(sF)/E_(gamma_Q) exp(sF).

Its potential has Hessian Q-s Hess F>=I_r/2. It is even, hence centered.
The Brascamp--Lieb variance inequality applies after every linear tilt,
whose Hessian is unchanged. Therefore every tilted covariance is at most
2I_r. Integrating the log-Laplace Hessian from zero proves that rho_s
is 2-subGaussian, uniformly in s, Q, and r.

For z_i=sqrt(beta)u_i.g this implies
P_(rho_s)(|z_i|>t)<=2 exp[-t^2/(4beta p_i)]. Layer-cake integration gives

    E_(rho_s) z_i^4<=64 beta^2 p_i^2,
    E_(rho_s) F<=(32/3)beta^2 S2.

The zero-row case is exact and needs no division by p_i. Integrating
d/ds log E_(gamma_Q)exp(sF)=E_(rho_s)F yields

    log E_(gamma_Q) exp(F)<=(32/3)beta^2 S2.              (19)

The variance principle used here is the classical Brascamp--Lieb bound
Var(f)<=E[grad(f)^T(Hess V)^(-1)grad(f)]. A primary proof, including the
uniformly convex approximation and integration-by-parts argument, is
Carlen--Cordero-Erausquin--Lieb,
[Asymmetric covariance estimates of Brascamp--Lieb type](https://arxiv.org/pdf/1106.0709),
Theorem 1.1 at p=2 and Section 2. That proof was read for this application.
Its regularity conditions hold for the C^2 uniformly convex potentials
above; alternatively smooth approximation preserves the uniform Hessian
bounds. This use is on the auxiliary continuous Gaussian interpolation,
not an unproved discrete-cube Brascamp--Lieb assertion.

### 8.3 The bad Fourier region and partition lower bound

Return to the arbitrary external direction t from Section 2 and its
reference Q=I_r+beta U^TDU>=I_r. On

    G_c={max_i |sqrt(beta)u_i.g|<=c},

(6) and f_c(z)=z^4 give

    product_i |cos(z_i)+i m_i sin(z_i)|
      <=exp[-beta g^TU^TDUg/2+F(g)].

Integrating and using (19) bounds the good-region contribution by
A_t exp(C_v S2), where A_t is (8). On the complement, (5) and a union
bound give

    P(G_c^c)<=2n exp[-c^2/(2beta p_*)].

After dividing by A_t>=v^(r/2), this tail is at most

    tau'_n=2n exp[(r/2)log(1/v)-c^2/(2beta p_*)].         (20)

It vanishes faster than any power of n under p_*<=Lr/n and
r=o(sqrt(n)). Thus, uniformly in ALL t,

    N_t<=[exp(C_v S2)+tau'_n] A_t.                       (21)

For clarity, the partition lower bound (3) also holds throughout this
larger range independently of covariance approximation. On G_c the
unbiased cosine product is positive and its difference from
exp[-beta||g||^2/2] is bounded by that Gaussian weight times
C_v sum_i z_i^4. The normalized integral error is O_v(S2). Outside G_c
bound both original integrands by one, and divide by v^(r/2), paying
(20). Consequently

    Z>=v^(r/2)[1-C_v(S2+tau'_n)],                        (22)

and the bracket tends to one since S2<=Lr^2/n=o(1). No r=o(n^(1/3))
covariance estimate has been reused to justify this step.

### 8.4 Completing the improved MGF estimate

Insert (21)--(22) into the exact Fourier identity (4). The determinant
comparison (11) is unchanged. We obtain, for all t,

    log E_nu exp(t.h)<=a_n||t||^2+b'_n,
    a_n=(1+beta p_*)/2,
    0<b'_n<=C_v(S2+tau'_n)+n^(-4).                       (23)

The positive n^(-4) is harmless and lets the same symmetry lemma apply
without a zero-defect case distinction. Equations (14) and (23) give

    K_n=1+beta p_*+O_v(sqrt(S2+tau'_n+n^(-4)))
       =1+O_(v,L)(r/sqrt(n)),                            (24)

which proves (16). The bound remains uniform with an actual rank q<=r
and leverage p_*<=r/(eta n), for fixed eta: replace S2 by its upper bound
r^2/(eta n), and the tail's dimension q by the rank budget r.

This completes the stronger cold-MGF ingredient. The enlarged-rank
exact-isotropy realization additionally needs the separately audited
covariance operator/Frobenius improvement; (24) alone does not imply it.
