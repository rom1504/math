# Correlated Gaussian rounding: polynomial slack gain and endpoint rigidity

Date: 2026-09-06. Status: independently derived by the director and
response agent; all constants and limiting endpoint conventions were
independently reconstructed by the bound-audit agent. This is a finite
theorem for arbitrary bounded-operator hollow signings. It strengthens
the earlier independent-rounding slack bound and uses no Haar law.

## 1. Finite statement

Let A be a symmetric hollow n-by-n sign matrix, n>=2, and set

    B=A/sqrt(n-1), ||B||op<=L, Q_B(x)=x^T Bx/2,
    Lambda(B)=max_(x in {+1,-1}^n)|Q_B(x)|/n.

For EVERY u in [-1,1]^n, define

    e_B(u)=Q_B(u)/n, d=1-||u||^2/n.

Then

    Lambda(B) >= |e_B(u)|+d^2/(4pi L)-err_(n,L),             (1)
    err_(n,L)=1/[pi L(n-1)]+1/[2L^2 sqrt(n-1)].

In original normalization Lambda(B)=Qabs(A)/(n sqrt(n-1)), where
Qabs(A)=max_x |x^T A x|/2. This improves the spin vector for the same
fixed matrix; it does not change the coefficient signing.

More precisely, for either chosen orientation B'=+B or -B there is a
correlated Boolean rounding v with E v_i=u_i EXACTLY and

    E Q_(B')(v)/n
      >= e_(B')(u)+(2/L)[(1/n)sum_i phi(t_i)]^2-err_(n,L), (2)
    t_i=Phi^(-1)((1+u_i)/2).

The conventions t_i=+/-infinity and phi(t_i)=0 apply at u_i=+/-1.
The finite error in (2) can be slightly sharpened as in the proof.

## 2. Covariance-matched Gaussian rounding

The matrix I+B'/L is positive semidefinite and has diagonal one. Take
a centered Gaussian vector G with this covariance and put

    v_i=2 1{G_i<=t_i}-1.

This has mean u_i, including the endpoint conventions. For standard
normal N, the first normalized Hermite coefficient of
2 1{N<=t_i}-1 is -2 phi(t_i), and its centered L2 norm squared is
w_i=1-u_i^2. Hence, by the Gaussian Hermite expansion, for i!=j,

    Cov(v_i,v_j)=4 phi(t_i)phi(t_j) rho_ij+R_ij,
    rho_ij=B'_ij/L,
    |R_ij|<=rho_ij^2 sqrt(w_i w_j).                         (3)

Indeed the terms of degree k>=2 are bounded using |rho|^k<=rho^2 and
Cauchy--Schwarz on their Hermite coefficients. Degenerate Gaussian
pairs at |rho|=1 are covered by the same L2 expansion or continuity.

The normalized half-energy gain is

    (1/2n)sum_(i!=j) B'_ij Cov(v_i,v_j).

Its first-chaos term is nonnegative and equals

    (2/[nL(n-1)])[(sum_i phi(t_i))^2-sum_i phi(t_i)^2]
       >= (2/L)[avg_i phi(t_i)]^2-1/[pi L(n-1)],           (4)

using phi(t)^2<=1/(2pi). The remainder has absolute value at most

    (1/[2nL^2(n-1)^(3/2)]) sum_(i!=j) sqrt(w_i w_j)
       <= d/[2L^2 sqrt(n-1)]
       <= 1/[2L^2 sqrt(n-1)].                             (5)

The first inequality in the second line uses
2 sqrt(w_i w_j)<=w_i+w_j. Equations (4)--(5) prove (2). No independent
rounding or a conditional field CLT is being invoked: the chosen
coordinate correlations themselves produce the positive energy gain.

## 3. An elementary Gaussian-profile bound

Let I(p)=phi(Phi^(-1)(p)) on (0,1), extended by zero at the endpoints.
Differentiation gives I'(p)=-Phi^(-1)(p) and I''(p)=-1/I(p)<0.
It is concave, symmetric about 1/2, and I(1/2)=1/sqrt(2pi).
The chord from zero to 1/2 therefore yields

    I(p)>=sqrt(2/pi) min(p,1-p).

Applying this at p=(1+u_i)/2 gives

    phi(t_i)>=(1-|u_i|)/sqrt(2pi)
              >=(1-u_i^2)/(2sqrt(2pi)).                   (6)

Thus (2/L)(avg phi)^2>=d^2/(4pi L). Choose B'=B when e_B(u)>=0,
and B'=-B otherwise. The Boolean maximum of absolute energy dominates
the expected oriented energy, proving the finite statement (1).

## 4. Quantitative rigidity of nearly optimal fractional endpoints

If |e_B(u)|>=Lambda(B)-epsilon, then (1) implies

    d <= sqrt(4pi L [epsilon+err_(n,L)]).                  (7)

In particular, along any fixed operator cap, asymptotically optimal
fractional mean vectors have vanishing average slack. This is stronger
than the mere existence of a Boolean maximizer of a multi-affine
quadratic form. It is a finite quantitative statement about EVERY
near-optimal fractional mean vector.

## 5. Relation to the current policy and algorithmic scope

The successful certified pure policy has ternary F, binary H=1-|F|,
and terminal C=H sign(BF). With a Boolean tie convention, its two
endpoints +/-F+C are exactly Boolean, so d=0. Thus (1) does not itself
increase the banked .4333221116640807 lower bound.

Deliberately damping a Boolean vector x to u=(1-rho)x creates
d=2rho-rho^2 and loses d times its positive normalized energy e.
The guaranteed polynomial gain is d^2/(4pi L)<=d/(4pi L), which cannot
pay that loss at e=.433322... and L>=1. This limits this certificate,
not all possible correlated-rounding or damping schemes.

The exact construction samples a Gaussian with covariance I+B'/L.
It uses the full matrix covariance, not a new independent seed that
may be substituted into an already exposed response. No finite query
complexity claim is required for (1). A fixed-accuracy matrix-polynomial
approximation to the positive square root can be considered separately;
it must not be silently counted as a free one-query operation.
