# A fixed gap for scalar-affine Gaussian rounding of actual near-minimizers

Status: proved and independently reconstructed, including the absolute-energy
concentration strengthening. No conclusion about convergence.

Let A be a hollow symmetric sign matrix of order n. Write P=max H_A,
R=max(-H_A), W=(P+R)/2, Q=max(P,R). Let G=Id+sA be positive semidefinite
and X=sign(Z), Z a centered Gaussian vector with covariance G. Singular
global G is allowed because every marginal coordinate has variance one.

## 1. Exact mean inequality

For s>=0, the bivariate Gaussian sign identity gives exactly

    E H_A(X)=n(n-1) arcsin(s)/pi.

PSD forces s<=1/(-lambda_min A). The elementary spectral bound
R<=n(-lambda_min A)/2 implies s<=n/(2R). Consequently, whenever R>n/2,

    E H_A(X)<=n(n-1) arcsin(n/(2R))/pi.                 (1)

For s<=0 replace A by -A and R by P; this bounds the absolute mean.
There is no polarity choice made after observing X in this calculation.

Now let A_n be ANY sequence with limsup Q(A_n)/n^(3/2)<=u, where

    c=.4333221116640807,     u=.494515125.

The verified universal width lower bound gives W(A_n)>=(c-o(1))n^(3/2).
Since min(P,R)=2W-Q, both endpoints are at least
(2c-u-o(1))n^(3/2). Thus every admissible scalar s has |s|=O(n^(-1/2)),
uniformly over its choice, and (1) proves

    limsup |E H_A(X)|/n^(3/2)
       <=1/[2pi(2c-u)] < .427688.                       (2)

In particular its mean is separated from W(A_n) by more than
.005634 n^(3/2) asymptotically. The threshold for this deduction is
u<2c-1/(2pi c)=.4993540740..., crossed by this campaign's new upper
bound but not by the previous .499432211 bound.

This excludes the scalar-affine Gaussian separating measure used in the
Steiner-ETF ramp obstruction on actual near-minimizers. It does not prove
that another separating measure is absent or that a small-trace ramp exists.

## 2. Quadratic concentration for this particular Gaussian law

The following calculation upgrades absolute MEAN to expected ABSOLUTE
energy. It is specific to the affine covariance and its uniformly small
off-diagonal correlations, not an arbitrary Gaussian response law.

First, for every nonzero eigenvalue lambda with unit eigenvector v,

    ||v||_infinity <=sqrt(n-1)/|lambda|,
    Q(A)>=|H_A(v/||v||_infinity)|
          >=|lambda|^3/[2(n-1)].

The first inequality is row Cauchy--Schwarz; the second uses multilinearity
to bound the quadratic form on the whole cube [-1,1]^n by Q. Thus

    ||A||op <=[2(n-1)Q(A)]^(1/3).                      (3)

For four distinct indices and |s| small, Gaussian density differentiation
on the compact set of four-dimensional correlation matrices with eigenvalues
in [1/2,3/2] gives the uniform Taylor expansion

    E XiXjXkXl = (2/pi)^2 s^2
       [Aij Akl+Aik Ajl+Ail Ajk]+O(|s|^3).             (4)

All third derivatives are uniformly integrable: they are a bounded-degree
polynomial times a Gaussian density with uniformly bounded precision.
The constant and first-order terms vanish by coordinate parity. The only
second-order terms are the three disjoint matchings, each with coefficient
(E|N(0,1)|)^4=(2/pi)^2. This also proves (4) for singular global G since
the four-coordinate marginal remains positive definite.

Subtracting E XiXj E XkXl cancels the Aij Akl term; its additional arcsin
remainder is O(s^4). Multiplying by the two energy coefficients and summing
the remaining matching terms over distinct indices gives two fourth-cycle
sums. Each is Tr(A^4) with at most O(n^3) repeated-index terms removed.
Shared-index edge pairs contribute at most O(n^3 |s|+n^3 s^2), and identical
edges contribute O(n^2). Thus for a numerical constant C,

    Var H_A(X) <= C[n^2+n^3 |s|+n^4 |s|^3
                         +s^2 Tr(A^4)]               (5)

for all sufficiently small |s|. Factors from ordered/unordered edges are
absorbed only into the absolute C, not into the leading mean in (1).

For the actual near-minimizers of Section 1, (3) gives ||A||op=O(n^(5/6)).
Since Tr A^2=n(n-1), Tr A^4<=||A||op^2 Tr A^2=O(n^(11/3)). With
|s|=O(n^(-1/2)), (5) gives Var H_A(X)=O(n^(8/3)), uniformly in s. Hence

    E|H_A(X)| <= |E H_A(X)|+O(n^(4/3)),
    limsup E|H_A(X)|/n^(3/2)<.427688.                  (6)

Chebyshev also bounds the probability of any fixed normalized excess above
this ceiling by O(n^(-1/3)); the improvement below gives O(n^(-1/2)).
Mixtures chosen before sampling retain the
expectation bound by averaging. Adaptive maximization over samples or over
coupled Gaussian laws is NOT covered.

### 2.1 Classical interpolation sharpens the error to O(n^(5/4))

Let beta_R=max_(x,y in [-1,1]^n)|x^T A y|. The corresponding complex
bilinear norm is at most 2 beta_R: rotate the output to be real, then split
the two inputs into real and imaginary parts. For real unit vectors u,v,
omit their zero coordinates and apply the three-lines lemma to

    F(z)=sum_ij A_ij sign(u_i)sign(v_j)|u_i|^(2z)|v_j|^(2z).

On Re z=0 it is bounded by 2 beta_R; on Re z=1 by max|A_ij|, since
the two absolute input sums are one. At z=1/2 it is u^T A v. Therefore

    ||A||op^2<=2 beta_R max|A_ij|<=8W(A)               (7)

for hollow signs. The last step uses
x^T A y=2[H_A((x+y)/2)-H_A((x-y)/2)] and the interval [-R,P], hence
beta_R<=4W. This is classical complex interpolation, not a new operator
theorem. Its hypotheses agree with
[Riesz--Thorin, Theorem 21 in Tao's notes](https://terrytao.wordpress.com/2009/03/30/245c-notes-1-interpolation-of-lp-spaces/).

Now Tr A^4<=8W(A)n(n-1)=O(n^(7/2)); (5) gives Var H_A(X)=O(n^(5/2))
and the error in expected absolute energy is O(n^(5/4)). The earlier
elementary bound remains a valid independent, weaker proof. The constant
ceiling in (2)/(6) is unchanged. No bounded normalized spectral norm was
assumed: (7) only gives ||A||op=O(n^(3/4)).

## 3. What this does and does not remove

The lower-response and upper-construction theorems jointly force a real
gap for a familiar rounding method on actual near-minimizing signings.
It is not a limitation of the signings themselves: their true cap is at
least c n^(3/2). The successful marked/nonlinear law is outside this
scalar-affine covariance class. Nor does the result exclude anisotropic
Gaussian covariances or prove an all-order recurrence.

### 3.1 Signed-mean robustness under low-complexity covariance changes

There is a restricted extension, NOT a mean-absolute-energy extension.
Let G_0=Id+sA be PSD with |s|<=1/6, let G be any other correlation matrix, and
put D=G-G_0. For |u|<=1/6 and |v|<=1,
`|arcsin(v)-arcsin(u)|<=2|v-u|`: same-sign outward secants are maximized
at (1/6,1), bounded by 3pi/5<2; inward ones have derivative at most
6/sqrt(35)<2; opposite-sign secants average slopes from zero bounded by
pi/2. The exact Gaussian sign identity and off-diagonal Cauchy--Schwarz give

    |E H_A(sign N(0,G))-E H_A(sign N(0,G_0))|
       <=(2/pi) sqrt(n(n-1)) ||D||F.

Thus an o(sqrt(n)) Frobenius change cannot remove the signed-mean gap.
In particular a uniformly bounded-operator, rank-o(n) covariance correction
cannot do so. If the target signed-mean improvement is delta*n^(3/2)
and ||D||op<=K, then rank(D)>=(pi*delta/(2K)-o(1))^2 n is necessary.
The correlations G need not be affine, but the variance proof above does
not automatically extend to them. No claim about their expected ABSOLUTE
energy is made here. This corollary was independently checked, including
the full versus off-diagonal Frobenius factor.

Dependencies: the freshly audited universal width lower bound, this
campaign's strict all-order upper construction, the elementary Gaussian
sign identity, and the finite Taylor argument (4). No random-disorder
universality theorem, spectral flatness, or optimizer stationarity is used.

Independent reconstructions:
[mean and exact constants](decisive_audit_nearmin_scalar_affine_gram_exclusion_2026_09_07.md),
[four-sign remainder and concentration](decisive_independent_four_sign_covariance_audit_2026_09_07.md).
The latter gives the explicit remainder 55296|s|^3 by differentiating the
four-dimensional density three times. Root independently checked its
331776 integrated derivative bound and exact fourth-cycle cancellation.
The exact constant and finite cycle checkers are in computations/.
