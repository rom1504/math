# Full-response variational structure: purification coordinates, stationarity, and a moment dual

Date: 2026-09-06. This is an independent calculation after the full-noise
channel theorem. The purification lemma itself is due to the parallel
response audit; the first variation and the dual below are reconstructed
here. No stationarity condition is asserted to prove a global maximum.

## 1. A useful convex parameter space before purification

Let G=(G_1,...,G_d) be independent standard Gaussian tree coordinates,
and let h_j=U^{-1}G_j be their even orthonormal child features. Include
any ancestors needed to make these features measurable in the old
Gaussian probability space. Put

    Psi(k,t)=E_N |k+sqrt(t) N|,   t>=0.

For odd f and even q with |f|<=q<=1, define

    p=E q,  a_j=E[f G_j],  K=sum_j a_j h_j,
    t=p-|a|^2,
    J(f,q)=E[(1-q) Psi(K,t)].                         (1)

Here f,q are the conditional mean and the conditional support probability
of a ternary response, not necessarily the response itself. Let U be an
independent even Gaussian uniform gate and let epsilon be the product
of the signs of three further independent Gaussian coordinates. Set

    F'=sign(f) 1{U<=|f|}+epsilon 1{|f|<U<=q}.

Then E[F'|old]=f and E[(F')^2|old]=q. Global parity is odd. Each added
Gaussian coordinate has zero first-chaos correlation with F': the gate
is even, and every noise-coordinate correlation leaves two independent
mean-zero signs. Thus the added coordinates create no new first-chaos
coefficients. This realizes (1) exactly in an
enlarged finite Gaussian system. Finite approximation precedes every
matrix limit.

The simpler purification of an already given response takes q=|f| and
uses only an independent even gate. It preserves a and increases residual
variance from E f^2-|a|^2 to E|f|-|a|^2. The response audit establishes
equality of the relevant suprema after finite approximation. The larger
convex domain |f|<=q<=1 is particularly convenient for variations and
line search; it is an actual ternary randomized-gate domain, not a
relaxation with an unproved realizability gap.

If f depends on all old coordinates, a denotes its entire first-chaos
vector. All formulas below hold with convergent Hilbert-space sums.
For finite-dimensional computation retain a finite ancestor-closed bank.

## 2. Derivatives and the exact threshold constants

At t>0, set z=K/sqrt(t) and write phi,Phi for the standard normal density
and distribution function. Then

    Psi_k=2Phi(z)-1,   Psi_t=phi(z)/sqrt(t).

Define the scalar and the coefficient vector

    B=E[(1-q) phi(z)]/sqrt(t),
    A_j=E[(1-q)(2Phi(z)-1) h_j]-2B a_j.             (2)

For any bounded feasible direction (delta f,delta q), differentiation of
p and a gives

    DJ[delta f,delta q]
      =E[(A dot G) delta f +(B-Psi(K,t)) delta q].   (3)

The factor 2 multiplying B a_j comes from differentiating -|a|^2.
There is only B, not 2B, in the support score. The Gaussian integration
in Psi creates this distinction. For finite directions the derivative
is ordinary dominated differentiation near t>0. The Hilbert-space form
follows because |2Phi(z)-1|<=1 and the h_j are orthonormal.

Pointwise maximization of (3) over |f'|<=q'<=1 gives the best linearized
response. Write V=A dot G and D=|V|+B-Psi(K,t). Then

    q'=1, f'=sign V   where D>0;
    q'=0, f'=0        where D<0.                    (4)

At D=0 any mixture of these choices is allowed; when V=0, either sign
or a mixture is allowed if support is selected. All choices preserve
global parity because V is odd and D is even.

Every local maximum of (1) on the convex feasible domain, with t>0,
must maximize its own linearization. Hence it obeys (4), up to these
ties. To see necessity without assuming differentiable boundary paths,
join the proposed maximizer to the measurable pointwise maximizer (4)
by a convex segment. Unless (4) already holds, its initial derivative
is strictly positive. Conversely (4) only says that all feasible first
derivatives are nonpositive. The objective is not generally concave;
this condition is not global sufficiency or even second-order sufficiency.

## 3. A constructive monotone update, with a genuine stopping criterion

Starting from any (f,q) with t>0, compute (2) and the linearized maximizer
(f_best,q_best) of (4). Its Frank--Wolfe gap is

    gap=E[(|V|+B-Psi)_+ -V f -(B-Psi)q]>=0.          (5)

If gap>0, the direction toward this pair has strictly positive derivative,
so a sufficiently small positive step increases J. Choose a step by
backtracking or by maximizing the exact scalar function

    J((1-theta)f+theta f_best,(1-theta)q+theta q_best),
          0<=theta<=1.                             (6)

Every intermediate pair remains feasible and is exactly realizable by
an additional even randomized gate. The current value and the endpoint
value alone are insufficient to justify taking the full step: theta=1
can decrease the nonlinear objective. With an exact positive-gap test
and backtracking, monotone increase is justified. Vanishing gap is a
first-order certificate only, not a certificate of a global upper bound.

## 4. Exact dual for the fixed moment slice

Fix p,a and hence t=p-|a|^2 and K. Optimize (1) over |f|<=q<=1 satisfying
E q=p and E fG=a. The following gives an upper bound for every c in R
and b in R^d:

    J <= E Psi(K,t)-c p-b dot a
           +E[(|b dot G|+c-Psi(K,t))_+].            (7)

Proof: add c(Eq-p)+b dot(EfG-a), which is zero on the slice,
and maximize pointwise over |f|<=q<=1. The pointwise maximum is the
positive part in (7). This is a convex dual minimization in c,b for each
fixed p,a. It supplies useful independently checkable UPPER bounds on
moment slices, unlike the stationarity equations.

The infimum of (7) equals the moment-slice optimum, allowing limiting
multipliers at boundary moments. Indeed |f|<=q<=1 is weak-star compact
in L-infinity squared, and integration against 1,G,Psi is weak-star
continuous since these functions are integrable. Its finite-dimensional
moment/objective image is compact and convex. Supporting-hyperplane
separation of its hypograph, or the concave biconjugate identity for
this compact image, gives exactly (7). Interior slices have the usual
finite-multiplier optimality statement; boundary slices can require
an infimum rather than an attained minimum.

At a global or local stationary pair, (2) supplies the distinguished
dual multipliers c=B,b=A. The pair then attains the pointwise dual
maximum and is globally optimal ON ITS OWN FIXED (p,a) SLICE. This is
a genuine sufficiency statement, but it does not compare different
moments p,a. Global optimization can therefore be separated into a
finite moment search and a convex slice optimization.

## 5. The exact feasible moment body and positive residual variance

For finite isotropic Gaussian G, the feasible moment pairs are exactly

    0<=p<=1,    |a|<=m(p),
    m(p)=2 phi(Phi^{-1}(1-p/2)),                    (8)

with m(0)=0. Necessity follows by testing in the direction a/|a|:
|a|<=E[q |N|], and among 0<=q<=1 of mean p the upper tail of |N|
maximizes this integral. Sufficiency follows by taking q to be that
upper-tail indicator and f=r sign(N) q with r=|a|/m(p), along the
desired direction. These choices are even/odd as required.

For every p>0, m(p)^2<p. This is strict Cauchy--Schwarz because the
Gaussian magnitude is not constant on a positive-probability support.
Thus every nonzero feasible ternary moment slice has t>0: the derivatives
above do not encounter an exact zero-noise boundary except at p=0.
Small t may still require numerical care; no uniform positive lower
bound on t is asserted as p tends to zero.

## 6. A low-dimensional star-bank specialization

For the bank G_0 and independent star coordinates G_1,...,G_d, choose

    h_0=1,   h_r=normalized H_(2r)(G_0).

Then K=a_0+sum_(r>=1) a_r h_r(G_0) is a one-variable even polynomial,
while V=A_0 G_0+sum_(r>=1) A_r G_r. Conditional on G_0=x, V is Gaussian
with mean A_0 x and variance sum_(r>=1) A_r^2. Equation (4) is therefore
a two-dimensional Gaussian integral, irrespective of bank size.

At a nondegenerate stationary point the star components of a and A are
parallel: Gaussian integration by parts in G_r gives a_r=A_r C for a
common nonnegative scalar C determined by the conditional threshold.
This reduces the dimensionality of a stationarity solver, but does not
remove the nonlinear shape of K or justify discarding additional tree
coordinates in the unrestricted variational problem.

The full response functional is at most 1/2 by Cauchy--Schwarz:
J<=sqrt(E(1-q)^2) sqrt(p)<=sqrt((1-p)p)<=1/2.
The separately audited original-class fixed local-algorithm ceiling is
stronger for GFOM-realizable constructions. Neither bound proves that
the stationary equations above have a particular globally optimal value.

## 7. One new Gaussian innovation implements the UNRESTRICTED gradient

There is an exact generative update stronger than a fixed-bank search.
Maintain a finite orthonormal Gaussian frame G_1,...,G_d whose inverse
features h_j=U^{-1}G_j are measurable functions of that same frame.
The frame need not consist of canonical individual tree coordinates;
it may consist of deterministic orthonormal first-chaos combinations.
Ancestor-closed canonical banks supply initial examples. Suppose f,q
are measurable in the frame, and set

    u=(1-q)(2Phi(K/sqrt(t))-1),
    w=u-2B K,    V=U w.

Then the FULL first variation is

    DJ=E[V delta f +(B-Psi) delta q].               (9)

Its old-frame projection coefficients are A_j=E[w h_j], agreeing with
(2). But V need not lie in the old Gaussian frame. Define

    sigma^2=||w||_2^2-sum_j A_j^2.

If sigma>0, append the Gaussian coordinate and its even inverse feature

    G_new=(V-sum_j A_j G_j)/sigma,
    h_new=(w-sum_j A_j h_j)/sigma.                 (10)

The new coordinate is standard Gaussian and independent of the entire
old frame, because U is an isometry into Gaussian first chaos. The new
feature is an explicit old-frame even function, so the measurability
invariant is preserved and no circular definition is present. If sigma
is zero, no new coordinate is needed. The exact unrestricted linearized
maximizer is (4) with

    V=sum_j A_j G_j+sigma G_new.

Every new response and every line-search mixture depends on a finite
Gaussian frame, and its first-chaos projection lies in that frame. Thus
its new K is again a finite linear combination of the known h features.
This gives a closed generative variational algorithm, adding at most one
independent Gaussian per step. Each fixed finite number of steps admits
the earlier canonical-tree L2 approximation and matrix realization;
one must not let the number of steps grow with n without a new theorem.

A zero fixed-bank gap is not an unrestricted stationary certificate
when sigma>0. This missing Gaussian innovation is the precise structural
reason that a shallow finite-bank optimizer can stall well below the
hierarchical construction. At an unrestricted stationary point V=Uw is
already among the coordinates on which the threshold rule depends; the
fixed-moment dual still applies, but global optimality is not automatic.

For a single update, conditional on the old frame, V is Gaussian with
mean A dot G and variance sigma^2. Its threshold support probability,
conditional sign mean, and new-Gaussian first moment are elementary
normal tails and densities. Also h_new depends only on the old frame.
Consequently the exact line-search objective for this update requires
only integration over the OLD frame; the newly introduced Gaussian can
be integrated out analytically. This is useful for low-dimensional
certified starts, before the old-frame dimension becomes large.

## 8. An explicit positive gain from any nonzero innovation

The innovation improvement is quantitative without an unevaluated
high-dimensional integral. Define

    T(z)=phi(z)-z Phi(-z)>0,   z>=0.

This function is decreasing and convex. For every mean mu, variance
sigma^2>0, and real threshold b, direct one-dimensional normal integration
and the cases |mu|>=max(b,0) or |mu|<b give

    E[(|mu+sigma N|-b)_+]-(|mu|-b)_+
      >=sigma T((|mu|+max(b,0))/sigma).             (11)

The old-frame feasible pair has linearized value at most its pointwise
old-frame maximum. Apply (11) with mu=A dot G and b=Psi(K,t)-B, then
Jensen to T. Since E|A dot G|=sqrt(2/pi)|A| and
E max(Psi-B,0)<=E Psi<=sqrt(p), the unrestricted gap satisfies

    gap >=g0:=sigma T((sqrt(2/pi)|A|+sqrt(p))/sigma)>0. (12)

One can turn this derivative bound into an explicit actual increase.
Let (f_new,q_new) be the full linearized maximizer, a_new its first-chaos
vector, and p_new=E q_new. Put

    v=|a_new-a|,
    Lp=|p_new-p|+2 max(|a|,|a_new|)v,
    d0=t/2,
    M=2v+(2phi(0)Lp+4phi(0)v^2)/sqrt(d0)
          +2phi(1)v Lp/d0+phi(0)Lp^2/(2d0^(3/2)).  (13)

For the exact feasible convex line theta in [0,1/2],

    t_theta=(1-theta)t+theta t_new
                 +theta(1-theta)|a_new-a|^2>=t/2.

The derivative formulas

    |Psi_kk|<=2phi(0)/sqrt(d0),
    |Psi_kt|<=phi(1)/d0,
    |Psi_tt|<=phi(0)/(2d0^(3/2))

and |Psi_k|<=1, Psi_t<=phi(0)/sqrt(d0) yield |J''(theta)|<=M.
Only L2 norms of the inverse features enter this calculation. Therefore
the completely explicit step

    theta=min(1/2,g0/M)

improves J by at least theta*g0/2. If M=0, use theta=1/2 and the line
is linear. This bound can be tiny, but it is rigorous and removes a
separate numerical-certificate obligation for strict improvement whenever
the innovation variance is known positive.

## 9. One non-shallow numerical sanity check, not a certificate

The script `computations/resumed_bound_audit_single_hierarchical_update.py`
defines a two-coordinate star threshold response by displayed rounded
constants and recomputes its actual Gaussian moments. Its one full-
gradient update introduces a genuinely non-star feature depending on
the previous center indicator. The new Gaussian is integrated out
analytically; the remaining old-frame rectangle is normal-CDF quadrature.

At quadrature sizes (512,128) and (1024,256), respectively, the values
agree to approximately 2e-11:

    starting value          0.3748605582060,
    innovation variance     0.1368862510019,
    unrestricted gap        0.0553931707667,
    full endpoint theta=1   0.3733215796253,
    best line step theta    0.4385575482,
    best line value         0.3863515547518.

The full step actually decreases the objective; line search is not a
cosmetic safeguard. These numbers are exploratory and are below the
banked hierarchical lower bound, so they are not offered as a new
universal decimal certificate. Equations (11)--(13), rather than agreement
of quadratures, certify the structural strict-improvement conclusion.

## 10. What the moment dual does, and why a uniform finite tail is not free

For a fixed finite bank, (7) is a rigorous computable upper bound on each
moment slice, and strong duality removes the optimization over arbitrary
measurable f,q once p,a are fixed. Taking the supremum over the feasible
moment body (8) is an exact finite-moment characterization of that bank.
It is not automatically an upper bound on the unrestricted full functional.

For a particular response whose omitted first-chaos norm is eta, dropping
those coefficients changes the functional by at most

    (1+sqrt(2/pi)) eta,                             (14)

using the L2 contraction for K and the Lipschitz bound for the residual
standard deviation. Conditional expectation onto an ancestor-closed old
bank preserves the retained moments and support mass; whenever the
retained inverse features are measurable in that bank, it also preserves
the corresponding truncated objective. Thus (7) plus (14) is a valid
a-posteriori bound for a class with a separately proved uniform tail eta.

There is no uniform small canonical tail over all bounded feasible
responses. Given any finite bank and any c in (0,1), choose an unused
Gaussian tree coordinate G_T and set f=c sign(G_T), q=c. This is feasible
and its first-chaos mass c sqrt(2/pi) lies entirely outside that bank.
Even selecting the largest finitely many coefficients does not uniformly
solve the issue: replace G_T by a normalized sum of arbitrarily many
unused coordinates. Boundedness and Gaussian L2 compactness do not yield
the missing tail estimate. A global dual bound obtained by finite
truncations therefore needs an additional coercivity or structural
argument; it cannot be claimed from the currently banked approximation
lemmas alone.

The rigorous obligations removed by this note are narrower but concrete:
functional optimization within a fixed moment slice is exactly dualized;
an unrestricted ascent direction is generated with one Gaussian rather
than an unspecified infinite feature search; and nonzero innovation
gives an explicit positive gain without high-dimensional quadrature.
The companion finite-causal stationarity theorem rules out high-value
termination of this generative process at any finite causal frame.
