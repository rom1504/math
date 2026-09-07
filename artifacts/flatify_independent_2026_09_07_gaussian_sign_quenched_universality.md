# Quenched pressure universality for bounded-spectrum Gaussian signs

2026-09-07. Director/independent-agent synthesis, using conditional
Bernoulli replacement and a block-local Gaussian Stein kernel. This
transfers the ACTUAL quenched spin pressure to its covariance-matched
Gaussian process. It does not by itself compare that Gaussian process
to the required original parent target.

## 1. General theorem

Fix epsilon>0 and K<infinity. Let G be a centered Gaussian vector in
R^d with correlation matrix epsilon I<=R<=KI. Put C=sign G and let
Y be a centered Gaussian vector with covariance

    S=Cov(C)=(2/pi)arcsin[R].

Let Omega be ANY finite configuration set, let u_i(omega) be arbitrary
real features with |u_i(omega)|<=1, and let b(omega) be arbitrary real
deterministic offsets. Define

    Phi(z)=log sum_(omega in Omega)
                exp[b(omega)+lambda sum_i z_i u_i(omega)].

Then

 |E Phi(C)-E Phi(Y)|
    <=C_(epsilon,K) |lambda|^3 d [log(d+1)]^(3/2).     (1)

The bound is uniform in the size and structure of Omega and in all
offsets. In particular it does NOT use a first-moment union over spin
pairs or a small Hessian bound on the pressure itself.

## 2. Conditional independent replacement

Write G=W+sqrt(epsilon)Z, with Cov(W)=R-epsilon I and Z independent
standard normal. Set

    f(w)=2Phi_normal(w/sqrt(epsilon))-1,
    sigma(w)=sqrt(1-f(w)^2).

Conditional on W, the coordinates C_i are independent, with means
f(W_i) and variances sigma(W_i)^2. Replace them one by one by

    F_i=f(W_i)+sigma(W_i)xi_i,                         (2)

where xi_i are independent standard Gaussians independent of W.
The first two conditional moments match. Every single-coordinate
third derivative of Phi has magnitude at most8|lambda|^3: it is the
third Gibbs cumulant of a feature bounded by one. The centered sign
variables have uniformly bounded third absolute moments and the
replacement Gaussians have variances at most one. Taylor's formula
therefore gives, conditionally and then unconditionally,

    |E Phi(C)-E Phi(F)|<=C|lambda|^3 d.                (3)

Also E F_i=0 and Cov(F)=S EXACTLY. For distinct coordinates the
independent xi's leave E f(W_i)f(W_j), the same conditional-sign
covariance; on the diagonal f^2+sigma^2=1.

## 3. Local smoothness, including Gaussian auxiliary tails

Group U_i=(W_i,xi_i) into d Gaussian blocks of dimension two. Their
joint covariance T has blocks

    T_ij=diag((R-epsilon I)_ij, 1_(i=j)),

so ||T||op<=max(K,1). The block variances may be singular in the W
coordinate; the Gaussian identities below hold by a square-root
representation or regularization and passage to the limit.

Let g(w,xi)=f(w)+sigma(w)xi. Its gradient and Hessian satisfy

    ||grad g(w,xi)||+||Hess g(w,xi)||op
                         <=C_epsilon(1+|xi|).         (4)

For clarity, the nontrivial scalar point is boundedness of sigma'
and sigma''. Write z=w/sqrt(epsilon). Then
sigma=2sqrt(Phi_normal(z)Phi_normal(-z)). On compact z intervals
smoothness is immediate. For z>=1, the elementary Mills inequalities
phi(z)z/(1+z^2)<=Phi_normal(-z)<=phi(z)/z show that sigma decays at
least as C exp(-z^2/4)/sqrt(z). Its logarithmic first derivative is
O(1+z), and its logarithmic second derivative is O(1): the inverse
Mills ratio r(z)=phi(z)/Phi_normal(-z) satisfies z<=r(z)<=z+1/z and
r'=r(r-z)<=2 for z>=1. Hence sigma',sigma'' are bounded. Evenness
treats negative z, and rescaling gives constants depending on epsilon.
The derivatives of f are bounded directly from its Gaussian density.

No truncation of xi is needed. Its Gaussian maximum moments will be
paid explicitly after the pointwise contraction below.

## 4. Exact multivariate Stein kernel with coordinate-local factors

Let P_s denote the OU semigroup for the Gaussian block vector U. Put

    a_i(U_i)=integral_0^infinity e^(-s)P_s grad g(U_i) ds,
    b_i(U_i)=grad g(U_i).

The transition marginal of block i depends only on U_i, so a_i is a
function of that block alone. From (4) and the derivative identity
grad P_s=e^(-s)P_s grad,

    ||a_i||+||D a_i||+||b_i||+||D b_i||
                          <=C_epsilon(1+|xi_i|).       (5)

Define

    Gamma_ij(U)=a_i(U_i)^T T_ij b_j(U_j).

The OU covariance identity gives the exact integration by parts

    E[F_i psi(F)]=sum_j E[Gamma_ij partial_j psi(F)],
    E Gamma_ij=S_ij.                                  (6)

One way to check the orientation in (6) is to use the OU resolvent
of F_i, whose gradient is a_i and is supported on block i; the
gradient of psi(F) on block j is b_j partial_j psi. Their covariance
metric is T_ij. Gamma need not be symmetric; Hessians below are
symmetric, and (6) requires no artificial symmetrization.

## 5. Gaussian comparison and the two-star contraction

Take Y independent of U with covariance S, and set

    Z_t=sqrt(t)F(U)+sqrt(1-t)Y,       0<t<1.

Differentiating E Phi(Z_t), using (6) and ordinary Gaussian integration
by parts for Y, gives

    d/dt E Phi(Z_t)
       =(1/2)sum_ij E[(Gamma_ij-S_ij) Phi_ij(Z_t)].     (7)

For fixed Y, apply the Gaussian OU covariance identity to Gamma_ij
and B_ij(U)=Phi_ij(Z_t). Its derivatives satisfy

    partial_(U_k,delta) B_ij
         =sqrt(t) b_(k,delta)(U_k) Phi_ijk(Z_t).

The derivative of Gamma_ij is supported ONLY on blocks i and j.
Writing U_s=e^(-s)U+sqrt(1-e^(-2s))U' with an independent copy U',
the two types of terms in the sum of covariances in (7) are therefore

  integral e^(-s) sum_(i,j,k,alpha,beta,gamma,delta)
      (T_ij)_(alpha,beta)(T_ik)_(gamma,delta)
      E[(partial_gamma a_i,alpha)(U_i)
         b_j,beta(U_j) b_k,delta(U_s,k)
         sqrt(t) Phi_ijk(Z_t(U_s,Y))] ds,              (8)

and the analogous expression with the second T block centered at j,
arising when the derivative hits b_j. All Greek indices range over
two values. This finite block-coordinate expansion avoids hiding a
dimension-dependent operator estimate in the notation.

Exactly Phi_ijk=lambda^3 Cum_Gibbs(u_i,u_j,u_k). Expand this third
cumulant into its five expectation products, with total absolute
coefficient sum six. Each is represented using at most three Gibbs
replicas. For fixed U,U_s,Y and fixed replicas, the corresponding
feature factors have the form v_i w_j z_k with every coordinate
bounded by one.

For each fixed set of Greek indices, the star sum in (8) thus has
the form

    sum_i d_i v_i
         [T^(alpha,beta) diag(e) w]_i
         [T^(gamma,delta) diag(f) z]_i.

Every d-by-d block-coordinate submatrix T^(alpha,beta) has operator
norm at most ||T||op. By Cauchy--Schwarz, the displayed sum has
magnitude at most

    ||d||infinity ||T||op^2 ||e||infinity||f||infinity d.

The same estimate applies to the star centered at j. It is pointwise:
the Gibbs replicas may depend arbitrarily on the Gaussian environments.
No independence between the spins and the local coefficients is used.

By (5), the product of the three coefficient maxima is bounded by
C_epsilon M_s^3, where

    M_s=1+max_i |xi_i|+max_i |xi_(s,i)|.

Both auxiliary arrays have standard Gaussian marginals; their possible
dependence is irrelevant. The elementary Gaussian union tail and its
integral give uniformly in s

    E M_s^3<=C[log(d+1)]^(3/2).                       (9)

Combining the finitely many Greek-coordinate and cumulant terms,
integrating e^(-s), and then integrating sqrt(t) from zero to one
proves

    |E Phi(F)-E Phi(Y)|
       <=C_(epsilon,K)|lambda|^3 d[log(d+1)]^(3/2).   (10)

All differentiations are justified first with smooth cutoffs and then
by Gaussian moment domination from (4)--(5); the final bound (9) is
independent of those cutoffs. Equations (3) and (10) establish (1).

### Explicit gap dependence

For 0<epsilon<=1, (4)--(5) can be sharpened to gradient bounds
C epsilon^(-1/2)(1+|xi|) and Hessian bounds
C epsilon^(-1)(1+|xi|). Each star in (8) uses exactly one Hessian
factor and two gradient factors. Consequently the constant in (1)
can be replaced by the explicit dependence

    C K^2 epsilon^(-2)|lambda|^3 d[log(d+1)]^(3/2),   (1a)

with an absolute C. The conditional replacement (3) is absorbed.
There is no small-source requirement in this pressure comparison;
all Gibbs derivative bounds hold at every lambda.

## 6. Actual two-child pressure and ground-state consequence

For a bridge with n-by-n entries, take d=n^2 and configurations
omega=(sigma,x,y), sigma=+/-1 and x,y Boolean. Set

    u_ij(omega)=sigma x_i y_j,
    b(omega)=beta sigma[H_A(x)+H_D(y)]/sqrt(n),
    lambda=beta/sqrt(n).

The offset is deterministic but otherwise unrestricted; A,D may be
actual optimizing children. The theorem gives the ABSOLUTE-cap pressure
comparison

    |E log Z_C(beta)-E log Z_Y(beta)|
       <=C_(epsilon,K) beta^3 sqrt(n)[log(n+1)]^(3/2). (11)

This pays the configuration entropy through quenched pressure, not
through a separate pairwise union. For either disorder law, dividing
log Z by beta/sqrt(n) differs from its absolute parent maximum by at
most (2n+1)log 2 times sqrt(n)/beta. Consequently their expected
normalized absolute parent caps differ by at most

    O(1/beta)
      +C_(epsilon,K) beta^2 n^(-1/2)[log(n+1)]^(3/2). (12)

For fixed epsilon,K, beta=n^(1/6)/sqrt(log(n+1)) balances the two
terms, giving normalized error O_(epsilon,K)(n^(-1/6)sqrt(log(n+1))).
Thus the expected ground-state problem for this ACTUAL Gaussian-sign
bridge is asymptotically equivalent to the Gaussian bridge with its
EXACT sign covariance, including the unchanged deterministic children.

For opposite-spectral R=I-rho A tensor D/p, with
p=||A||op||D||op and fixed rho<1, epsilon=1-rho,K=1+rho and

    Cov(vec Y)=I-[(2/pi)arcsin(rho/p)] A tensor D.

No flat-child hypothesis is required for this transfer. It does not
show that the Gaussian parent maximum meets 2sqrt(2)M_n. That is a
new, precise Gaussian-process optimization obligation, and no favorable
comparison or recurrence is inferred from the covariance formula alone.

## 7. A paid moving-gap approach to the endpoint Gaussian process

For the p-normalized opposite-spectral law choose

    epsilon_n=1-rho_n=n^(-1/7),
    beta_n=n^(1/14)/sqrt(log(n+1)).

Equation (1a) and the soft-max error show that the normalized expected
parent maxima for the actual rho_n sign law and its covariance-matched
Gaussian law differ by O(n^(-1/14)sqrt(log(n+1))).

Let S_1=I-[(2/pi)arcsin(1/p)]A tensor D denote the endpoint GAUSSIAN
covariance. Since p>=n-1, ||S_(rho_n)-S_1||op=O(epsilon_n). For all
Boolean rank-one features of norm n, covariance entries therefore
change by at most O(epsilon_n n^2). Standard Gaussian max comparison
(add independent Gaussian noise to dominate the increment variance
error, then apply Sudakov--Fernique) bounds their expected maxima
difference by O(n^(3/2)sqrt(epsilon_n)); the index set has logarithm
O(n), including the absolute-cap polarity. Deterministic child-energy
offsets are identical on both sides and do not change this comparison.

Thus actual nonsingular rho_n sign laws realize the endpoint Gaussian
parent maximum up to O(n^(10/7)sqrt(log(n+1))) in expectation, uniformly
over the original full-sign children. This is NOT a direct theorem
about the singular rho=1 Gaussian-SIGN law, nor a proof that the
endpoint Gaussian maximum is sufficiently small.
