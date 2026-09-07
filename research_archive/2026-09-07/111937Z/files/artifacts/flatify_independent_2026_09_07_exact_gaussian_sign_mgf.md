# Exact-covariance small-field MGF for Gaussian signs

2026-09-07. Director's smoothing/Stein argument, independently reconstructed
by the independent, adversary, and construction agents. This is the
canonical proof. It supersedes the Holder-inflated bound for covariance
matrices with a fixed positive spectral gap. It supplies a genuine
energy-sensitive exponential bound for an actual sign bridge, not a
convergence theorem.

## 1. The dimension-uniform theorem

Fix 0<epsilon<=1 and K>=1. There are constants c_(epsilon,K)>0 and
C_(epsilon,K)<infinity such that, for EVERY dimension d, correlation
matrix R satisfying

    epsilon I <= R <= K I,

centered Gaussian G with covariance R, and real h with
||h||infinity<=c_(epsilon,K),

 | log E exp<h,sign G> - (1/2) h^T S h |
       <=C_(epsilon,K) ||h||infinity ||h||^2,           (1)

where S=Cov(sign G)=(2/pi)arcsin[R], entrywise. Neither entrywise
small correlations nor a bound on ||Rh||infinity is assumed.

In particular, for d=n^2 and h_ij=t x_i y_j/sqrt(n), with t bounded
and x,y Boolean, the error is O_(epsilon,K,t)(sqrt(n)), uniformly in
ALL spin pairs. The spectral-gap constants are fixed before n grows.

## 2. A smooth Gaussian small-Hessian lemma

Let Z be standard Gaussian in R^d. Suppose a smooth real function U
has global bounds

    sup||grad U||<=L,       sup||Hess U||op<=H<1/2.

Then an absolute constant C satisfies

 |log E exp U - E U -(1/2)Var U| <=C H L^2.            (2)

Linear growth of U guarantees all exponential integrability needed
below. Smooth approximation extends the statement to the corresponding
C^2 functions; the functions used in Section 3 are already smooth.

### Stein identity and a global derivative bound

Write P_s for the Ornstein--Uhlenbeck semigroup, and define

 tau_U(z)=integral_0^infinity e^(-s)
              <grad U(z),P_s grad U(z)> ds.

Gaussian integration by parts, or the OU resolvent of U-EU, gives

    E[(U-EU) phi(U)]=E[tau_U phi'(U)],
    E tau_U=Var U.                                    (3)

For completeness, -L_OU^(-1)(U-EU)=integral_0^infinity P_s(U-EU)ds;
differentiating uses grad P_s=e^(-s)P_s grad. Gaussian integration by
parts against phi(U) then gives (3). The bounded gradient permits
dominated convergence throughout the resolvent differentiation.

Differentiating the displayed formula for tau_U and using the same
semigroup derivative identity gives

    sup||grad tau_U||
       <=HL integral_0^infinity(e^(-s)+e^(-2s))ds
       =(3/2)HL.                                     (4)

### A proof of the tilted Poincare bound used here

For 0<=t<=1 let mu_t have density proportional to exp(tU(z)) against
standard Gaussian measure. Its potential is
V_t(z)=||z||^2/2-tU(z), with Hess V_t>=(1-H)I. Thus

    Var_(mu_t) f <=(1-H)^(-1) E_(mu_t)||grad f||^2.    (5)

One direct justification, recalled to specify the imported convexity
step, uses the Langevin generator L_t=Delta-grad V_t dot grad. Its
diffusion flow contracts gradients at rate 1-H (differentiate the
flow, whose derivative obeys J'=-Hess V_t J). Consequently
|grad P_s^t f|^2<=exp[-2(1-H)s]P_s^t|grad f|^2.
The variance identity
Var_mu f=2 integral_0^infinity E_mu|grad P_s^t f|^2 ds
then proves (5). Approximation by smooth bounded functions gives the
globally Lipschitz functions U and tau_U used below.

By (4), (5), and covariance Cauchy--Schwarz,

 |Cov_(mu_t)(tau_U,U)|
     <=(1-H)^(-1)(3/2)HL^2 <=3HL^2.                  (6)

### Integrating the log MGF

Let ell(t)=log E exp(tU). Equation (3), applied to exp(tU), yields

    ell'(t)=EU+t E_(mu_t) tau_U.

Also d/dt E_(mu_t)tau_U=Cov_(mu_t)(tau_U,U). Therefore (6) implies

    |E_(mu_t)tau_U-Var U|<=3tHL^2.

Integrate ell'(t) from zero to one to obtain (2), for example with
constant one in this last bound after integrating 3t^2. Only an
absolute constant is needed subsequently.

## 3. Independent Gaussian smoothing of every sign coordinate

Decompose

    G=W+sqrt(epsilon)Z',

where W is centered Gaussian of covariance R-epsilon I and Z' is an
independent standard Gaussian vector. W may be singular. Write
W=(R-epsilon I)^(1/2)Z with a standard Gaussian Z; all derivative
bounds below are taken after this whitening representation.

Set

    f(w)=2Phi(w/sqrt(epsilon))-1,
    F(W)=sum_i h_i f(W_i),
    V(W)=(1/2)sum_i h_i^2[1-f(W_i)^2].

Conditional on W the sign coordinates are independent, with means
f(W_i). The scalar log MGF of a +/-1 variable of arbitrary mean has
its third derivative uniformly bounded for sources in a fixed compact
interval. Consequently, POINTWISE in W,

 log E_(Z') exp<h,sign G>
      =F(W)+V(W)+E_h(W),
    |E_h(W)|<=C sum_i |h_i|^3
                  <=C ||h||infinity||h||^2.            (7)

Put U(Z)=F(W(Z))+V(W(Z)). Since f and its first two derivatives are
bounded by constants depending only on epsilon, and
||R-epsilon I||op<=K,

    sup||grad U||<=C_(epsilon,K)||h||,
    sup||Hess U||op<=C_(epsilon,K)||h||infinity         (8)

when ||h||infinity<=1. Specifically, before whitening the Hessians of
F and V are diagonal, bounded respectively by C_epsilon max|h_i| and
C_epsilon max h_i^2; whitening multiplies their operator bound by at
most K. Gradient bounds follow by summing the squared coordinate
derivatives. Thus the smallness constant in (1) makes H<1/2 in (2).

Equations (2), (7), (8) give

    log E exp<h,sign G>
      =EU+(1/2)Var U+O_(epsilon,K)(||h||infinity||h||^2). (9)

## 4. Identifying EXACTLY the covariance in the leading term

Under the centered Gaussian W, F is odd and V is even. Hence

    EF=0,       Cov(F,V)=0.

Gaussian Poincare, or (5) at t=0 after whitening, gives

    Var V<=C_(epsilon,K)||h||infinity^2||h||^2.

If B=<h,sign G>, conditional expectation and variance are exactly

    E[B|W]=F(W),       Var(B|W)=2V(W).

The law of total variance therefore gives

    Var B=Var F+2EV.

It follows that

    EU+(1/2)Var U=(1/2)Var B
                  +O_(epsilon,K)(||h||infinity^2||h||^2).

Substitute in (9), absorbing the smaller error, to obtain (1). This
also explains why no independent residual or first-chaos approximation
is made: the leading variance is the covariance of the ACTUAL signs.
The arcsine identity for S follows from planar Gaussian rotation for
each pair of coordinates, including its diagonal value one.

### Explicit spectral-gap dependence (adversary audit refinement)

The smoothing derivatives satisfy |f'|<=C epsilon^(-1/2) and
|f''|<=C epsilon^(-1). Tracking (8) gives

    L<=C sqrt(K/epsilon)||h||,
    H<=C(K/epsilon)||h||infinity.

Thus the entire error in (1) is bounded by

    C K^2 epsilon^(-2)||h||infinity||h||^2,            (1a)

with an absolute C, provided ||h||infinity<=c epsilon/K. The extra
Var V error is only C K epsilon^(-1)||h||infinity^2||h||^2 and is
absorbed in (1a); the conditional third-order error is universal.
This supplies an explicitly paid moving-gap version rather than
silently assuming fixed-gap constants remain bounded.

## 5. Opposite-spectral full-sign bridge: exact flat specialization

Let A,D be hollow symmetric full-sign matrices of order n satisfying

    A^2=D^2=(n-1)I.

For a fixed 0<=rho<1, sample a Gaussian matrix G whose vector covariance
is

    R_rho=I-rho A tensor D/(n-1),

and put C_ij=sign G_ij. This is an actual full-sign bridge. The spectrum
of R_rho lies in [1-rho,1+rho], so Section 1 applies with fixed
epsilon=1-rho, K=1+rho. Since every nonzero off-diagonal correlation
equals +/-rho/(n-1), the sign covariance is EXACTLY

    S_rho=I-kappa_n A tensor D,
    kappa_n=(2/pi)arcsin(rho/(n-1)).                  (10)

For e_A=H_A(x)/n^(3/2), e_D=H_D(y)/n^(3/2), define

    v_n(e_A,e_D)=1-4n kappa_n e_A e_D.

Then, uniformly in ALL Boolean x,y and bounded t,

 log E exp[t x^T C y/sqrt(n)]
       =(nt^2/2)v_n(e_A,e_D)+O_(rho,t)(sqrt(n)).       (11)

Here

    v_n(e_A,e_D)
       =1-(8rho/pi)e_Ae_D+O(1/n),                    (12)

uniformly over the child energies allowed by the flat spectral bound.
In particular v_n is bounded below by 1-2rho/pi+o(1)>0.

For bounded b>0, Chernoff with t=b/v_n and its negative gives

 P(|x^T C y|>=b n^(3/2))
   <=2exp[-n b^2/(2v_n(e_A,e_D))+O_(rho,b)(sqrt(n))]. (13)

The order of limits is essential: rho<1 is fixed, then n tends to
infinity. Afterwards rho may approach one. No uniform bound as
rho->1 is asserted, and (1) has not been proved for singular R_1 by
simply setting epsilon=0.

Alternatively, (1a) permits the concrete sequence rho_n=1-n^(-1/6).
For bounded t the MGF error is O_t(n^(5/6)); replacing rho_n by one
in its leading proxy also costs O_t(n^(5/6)). Thus the limiting rho=1
proxy is realized by actual nonsingular Gaussian-sign laws at every
n, with a fully paid o(n) error. This does not claim the same uniform
estimate directly for the singular endpoint law R_1.

The exact covariance improvement rewards SAME-sign internal energies;
this dependence survives in a proved exponential tail, not merely the
second moment. It is thus an actual correlated operation with a paid
joint-energy mechanism.

## 6. General children: normalized law, and what remains unpaid

### A direct law valid for EVERY child pair (adversary simplification)

Put p=||A||op||D||op, which is at least n-1 by the fixed Frobenius
norms. Then for any fixed rho<1,

    R=I-rho A tensor D/p

is a correlation matrix with spectrum in [1-rho,1+rho]. Its Gaussian
sign bridge has the EXACT covariance

    S=I-kappa A tensor D,
    kappa=(2/pi)arcsin(rho/p).

Consequently the proved MGF and tail hold with actual variance proxy

    v_n(e_A,e_D)=1-4n kappa e_A e_D.

There is no positive-carrier response and no bounded child-operator
assumption in this version. Its explicit price is attenuation:
4n kappa is approximately 8rho n/(pi p), rather than 8rho/pi.
The cap estimate ||A||op^2<=4Q(A) permits p=O(n^(3/2)), so this
reward can vanish in the worst allowed case. The flat case p=n-1
maximizes the coefficient. Fixed-gap or the paid moving-gap version
above applies without any further normalization hypothesis.

### The constant-diagonal carrier alternative

The director's constant-diagonal normalization works for arbitrary
hollow sign children. With s=sqrt(n-1), set

    K_A=|A|+diag(s-diag|A|),

and similarly K_D. The correction is PSD because |A|_ii<=s. Thus
K_A>=+/-A and diag K_A=s. The covariance

    R_0=(K_A tensor K_D-A tensor D)/(n-1)

is PSD and has diagonal one. For fixed rho<1 use

    R_rho=(1-rho)I+rho R_0.

If its operator norm is bounded by a fixed K, (1) applies directly
to the exact arcsine covariance of the actual bridge signs. This
preserves actual, unweighted child energies in its first Gaussian
covariance term. However, bounded K, the positive K_A,K_D responses,
and any approximation of higher arcsine powers all need separate
payment. They are not consequences of the cap bound alone in this
note. The theorem must not be specialized to arbitrary optimizing
children by silently declaring them spectrally flat.

## 7. Fully paid shell criterion and its cap-only limitation

For fixed children the original parent objective is exactly

    max_(x,y) [ |H_A(x)+H_D(y)|+|x^T C y| ].

Take a finite energy-bin partition before n grows. Suppose a bin pair
has child energy values in intervals I,J, state counts at most
exp[n(s_I+o(1))], exp[n(s_J+o(1))], and a common tail variance ceiling
v_IJ for all its spin pairs. Put

    q_IJ=max{|u+v|:u in I,v in J}.

At parent target L n^(3/2), this bin is paid whenever

    L>q_IJ,       (L-q_IJ)^2/(2v_IJ)>s_I+s_J.          (14)

Strict margins pay the O(sqrt(n)) errors, the two tail polarities,
and the finite union of bins. If every bin pair satisfies (14), a
single actual bridge controls the ORIGINAL parent cap. Counts here
are counts of actual child states, not typical-noise assertions.

The exact-covariance proxy alone, even in the ideal flat case and with
rho approaching one, cannot close (14) from the cap bound only. In
fact every hollow sign child has E_x H_A(x)=0 and
E_x H_A(x)^2=binom(n,2). Therefore for every fixed delta>0,

    #{x: |H_A(x)|<=delta n^(3/2)}=2^n[1-O(1/(delta^2 n))].

So both near-zero energy bins have entropy log 2+o(1), and their exact
covariance proxy tends to one as their widths delta decrease. At the
desired balanced target L=2sqrt(2)c with c<=1/2, their limiting
first-moment exponent is at most

    L^2/2<=1 < 2log 2.

Thus this fully paid pair-by-pair first-moment shell certificate fails
already on the exponentially large near-zero-energy sector. This is
not a lower bound on the actual selected parent cap: it pinpoints the
remaining requirement for overlap-sensitive supremum control or a
different selection argument. No recurrence or convergence follows
from (1)--(14) alone.
