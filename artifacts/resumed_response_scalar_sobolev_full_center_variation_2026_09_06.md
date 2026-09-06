# An attained scalar endogenous full-center optimum and its variational equation

Date: 2026-09-06. Status: exact derivation independently reconstructed by the
director through Section 10, including criticality and the scalar class gap.
This is a bounded variational class, not a proposed solution of the full
matrix optimization or a numerical parameter-tuning exercise.

## 1. The exact endogenous scalar domain

Let gamma denote standard Gaussian measure on R and N the nonnegative
Gaussian Ornstein--Uhlenbeck number operator, N h_r=r h_r for normalized
Hermites. Consider real jointly even functions g with

    ||g||_(L2(gamma))^2=1,
    D(g)=sum_(r even) r <g,h_r>^2<=1.              (1)

The derivative energy is allowed to equal one. The scalar recovery
theorem in `fresh_tree_fixed_point_critical_independent_2026_09_05.md`
gives a UNIQUE standard first-chaos Gaussian V in the canonical space
such that V=U g(V). This is an actual same-space frame, not a formal
distributional fixed point. The matrix and cyclic-anchor extensions are
in `resumed_response_matrix_recovery_independent_audit_2026_09_06.md`.

There is no zero-mean exception hidden in (1). Evenness gives

    D(g)>=2(1-(Eg)^2),   hence |Eg|>=1/sqrt(2).     (2)

The deterministic one-child obstruction possible for an ANCHORED scalar
innovation cannot occur for an unanchored even scalar g. Replacing g by
-g will leave the value below unchanged, so one may choose either
orientation component, for example Eg>=1/sqrt(2).

Every function in (1) has a continuous representative: Gaussian Sobolev
regularity gives ordinary one-dimensional H1 regularity on each compact
interval, since the Gaussian density is bounded above and below there.

## 2. Exact one-dimensional full-response functional

For alpha>0 put H(v)=1{|v|<=alpha}, and define

    p=E H=2Phi(alpha)-1,  q=1-p,
    c=E[g(V)H(V)],  d=sqrt(p-c^2)>0.

The strict inequality is automatic. Equality in Cauchy--Schwarz would
force g=+/-H/sqrt(p), which is not in Gaussian H1 for any finite positive
alpha because of its nonzero jumps.

Let W=U H(V). Then (V,W) is Gaussian with variances (1,p) and covariance
c, so W=cV+dZ for an independent standard Gaussian Z. For the actual
bounded odd response F=sign(W)1{|V|>alpha}, Gaussian integration gives

    a=E VF=lambda+c gamma,
    b=E ZF=d gamma,

where the coefficient gamma in the following formulas is a SCALAR,
not Gaussian measure:

    lambda=2phi(alpha)[2Phi(c alpha/d)-1],
    gamma=4phi(0)/sqrt(p) * Phi(-alpha sqrt(p)/d).  (3)

The covariance identities U*V=g and U*Z=(H-cg)/d yield

    K=U*F=lambda g+gamma H,
    t=||F-P1F||_2^2
      =q-lambda^2-2c lambda gamma-p gamma^2>0.     (4)

Consequently the independently proved full-center response theorem gives
the ACTUAL matrix lower certificate

    liminf_n M_n/n^(3/2) >= T_alpha(g),
    T_alpha(g)=integral_(-alpha)^alpha
          Psi(lambda g(v)+gamma,t) phi(v) dv,
    Psi(k,t)=E|k+sqrt(t)N|.                        (5)

This is an exact one-dimensional integral, not a conditional-Jensen
surrogate for a higher-dimensional inverse feature. In particular

    E HK=lambda c+gamma p

is the marked value of this mask, and (5) is its full Gaussian-noise
center improvement. The symmetry g->-g sends c,lambda to their negatives
and leaves gamma,t,lambda g, and T unchanged.

The inequality t>0 follows independently from the Gaussian moment body:
the bounded F has support probability q, so

    a^2+b^2<=m(q)^2<q,
    m(q)=2phi(Phi^(-1)(1-q/2)).                    (6)

Thus none of the normal-density derivatives below is taken at a
zero-noise boundary.

## 3. Strong compactness and an actually attained scalar optimum

Let S be the set in (1). It is strongly compact in Gaussian L2. Indeed
the Hermite tail beyond degree N is bounded by 1/(N+2) for even N, and
the coefficients in each finite initial segment lie in a compact ball.
Any strong limit retains norm one and has derivative energy at most
one by lower semicontinuity. Evenness is closed.

For every pair in (5), Cauchy--Schwarz gives

    0<=T_alpha(g)<=sqrt(pq)<=1/2.                 (7)

To see the first upper bound, Psi(K,t)^2<=K^2+t and
||K||_2^2+t=q; restrict the expectation to H and apply Cauchy--Schwarz.
Thus values tend uniformly to zero as alpha tends to zero or infinity.

On any compact alpha interval inside (0,infinity), the maps
(alpha,g)->p,c are continuous in strong L2, since the indicator H_alpha
is L2-continuous in alpha. The strict inequality d^2>0 from Section 2,
compactness of S, and compactness of the alpha interval imply a uniform
positive d margin. Likewise (6) gives a uniform positive t margin.
The coefficients in (3)--(4) are then continuous, and

    |Psi(k,t)-Psi(k',t')|
      <=|k-k'|+E|N| |sqrt(t)-sqrt(t')|

proves continuity of T, including its moving integration interval by
Cauchy--Schwarz. Hence

    C_scalar=sup_(alpha>0,g in S) T_alpha(g)       (8)

is attained by some finite alpha_* and g_* in S, provided it is positive.
Positivity is immediate, for example from g=1 and any nontrivial mask.
The optimizing g_* is an actual endogenous scalar inverse by Section 1.
The construction uses a two-Gaussian frame (V,W), though its canonical
tree-coordinate support can be infinite.

For a specified certified lower reference c_0>0, the superlevel set
T>=c_0 lies in the explicit alpha window determined by

    p_-<=p<=p_+,
    p_+/-=(1+/-sqrt(1-4c_0^2))/2.                 (9)

This makes the compactness argument quantitative once its d margin is
bounded. No assertion is made that C_scalar equals the full response
supremum, the GFOM ceiling, or the original matrix asymptotic value.

## 4. Quantitative finite-Hermite approximation

Let P_N g retain even Hermite degrees at most the even integer N, and
let e=||g-P_Ng||_2^2. Normalize

    g_N=P_Ng/sqrt(1-e).

The denominator is nonzero by (2). This normalization does NOT spoil
the derivative constraint. In fact

    D(P_Ng)<=1-(N+2)e,
    D(g_N)<= [1-(N+2)e]/(1-e)<=1,                 (10)
    ||g-g_N||_2<=2sqrt(e)<=2/sqrt(N+2).

If e>0 the normalized truncation is strictly subcritical. If e=0 and
D(g)=1, the exact critical recovery theorem still applies.

Let C_N be (8) restricted to even Hermite polynomials of degree at most
N. The classes are nested and every C_N is attained. Therefore

    C_N increases to C_scalar.                    (11)

There is a quantitative bound. Restrict alpha to the compact window
(9) at a fixed positive reference below the maximum. On this window
and S, let d>=d_0>0 and sqrt(t)>=tau_0>0. The coefficients (3)--(4)
have uniformly bounded derivatives in c, and the L2 gradient in
Section 6 is uniformly bounded by some explicit finite L_0. The segment
between g and g_N has norm at most one, and c stays between its endpoint
values, so the same coefficient bounds apply. Consequently

    0<=C_scalar-C_N<=2L_0/sqrt(N+2)               (12)

for all sufficiently large N for which C_N exceeds the chosen reference.
This is a convergence theorem for a genuine attainable certificate
class, NOT a convergence theorem for M_n/n^(3/2).

Here is one deliberately rough explicit choice of the fixed-alpha
Lipschitz constant. Write alpha<=alpha_max, p>=p_min, and put

    Lbar=4alpha_max phi(0)^2/d_0^3,
    L_0=2phi(0)+2Lbar
          +(4phi(0)^2/tau_0)[2phi(0)/sqrt(p_min)+Lbar].

The derivative identities in Section 6 imply this bound using
|lambda|<=2phi(0), |gamma|<=2phi(0)/sqrt(p_min), ||g||_2<=1,
and B<=phi(0)/tau_0. Sharper constants are unnecessary for (12).

## 5. The d margin is certifiable, not only existential

Let beta_r(alpha)=<H_alpha,h_r>. The exact coefficients are

    beta_0=p,
    beta_r=-2phi(alpha)h_(r-1)(alpha)/sqrt(r),
                                       r>=2 even.

For any a>0, weighted Cauchy--Schwarz and (1) give

    c^2 <= (a+1) sum_(r even) beta_r^2/(a+r).      (13)

At fixed nontrivial alpha, the weighted threshold derivative energy is
infinite. Moreover

    a[p-(a+1)sum beta_r^2/(a+r)]
      =-p^2+sum_(r>=2 even)(r-1)a beta_r^2/(a+r)

increases to infinity as a tends to infinity. On a compact alpha window,
the nested-open-set/finite-subcover argument therefore supplies one
finite a for which (13) gives a strictly positive uniform d^2 gap.

For a completely finite certificate, bound the sum in (13) above by

    sum_(r<=L even) beta_r^2/(a+r)
        +[p-sum_(r<=L even) beta_r^2]/(a+L+2).     (14)

The threshold family is compact in L2 on the alpha window, so these
upper bounds converge uniformly to the sum. A sufficiently large finite
L retains a positive uniform margin. Formula (14), the elementary
Hermite endpoint formulas, and rational interval subdivision in alpha
give a reproducible way to certify a numerical d_0 without evaluating
an infinite series. The t margin follows explicitly from (6).

Similarly, finite-dimensional coefficient and alpha nets together with
interval evaluation of the one-dimensional integral in (5) provide
two-sided certificates for C_N. Feasible polynomial points are dense:
slightly shrink all nonconstant coefficients, put the remaining norm
in the constant coefficient with fixed sign, then approximate the
nonconstant coefficients rationally inside the resulting strict energy
margin. This avoids accidentally violating (1) when discretizing a
critical boundary. Formula (12) supplies the tail error for C_scalar.

## 6. Exact first variation in g and the weak Euler--Lagrange equation

For fixed alpha, write

    z=c alpha/d,
    L=4alpha phi(alpha)phi(z)/d^3.

Direct differentiation of (3), using
phi(0)phi(alpha sqrt(p)/d)=phi(alpha)phi(c alpha/d), gives

    lambda_c=L p,
    gamma_c=-L c,
    t_c=-2lambda(gamma+L d^2).                    (15)

In particular the marked value lambda c+gamma p has derivative exactly
lambda; the apparently extra terms cancel. This is a useful independent
check on the coefficients.

Put K_c(v)=lambda g(v)+gamma on the center, let
s_c(v)=2Phi(K_c(v)/sqrt(t))-1, and define

    B=E[H phi(K_c/sqrt(t))]/sqrt(t),
    A=lambda_c E[H g s_c]+gamma_c E[H s_c]+B t_c.

Then the L2 derivative is

    D_g T[delta g]
      =<H[lambda s_c+A],delta g>.                 (16)

At a maximizing pair, the Hilbert-space multiplier conditions apply to
the norm equality and derivative-energy inequality. The constraint
qualification at D(g)=1 is explicit. If b_0=Eg, the even direction

    delta g=b_0-b_0^2 g

is tangent to the norm sphere and has strictly negative energy derivative
-2b_0^2D(g). Equation (2) ensures b_0!=0. Thus there are multipliers
kappa in R and eta>=0 such that

    kappa g+eta N g=H[lambda s(lambda g+gamma)+A], (17)
    eta(D(g)-1)=0,

in the weak Gaussian H1 sense. These equations are necessary, not
sufficient for a global maximum. No sign of kappa is assumed, and the
case eta=0 is not discarded without an additional argument.

If eta>0, the bounded right side makes g belong to the Gaussian
number-operator domain. In ordinary one-dimensional notation it solves

    eta(-g''+v g')+kappa g
       =lambda s(lambda g+gamma)+A,  |v|<alpha,
    eta(-g''+v g')+kappa g=0,          |v|>alpha.  (18)

The solution is even, has g'(0)=0, and g,g' match continuously at the
two threshold endpoints. The forcing can jump there, but there is no
delta term requiring a jump in g'. Gaussian L2/H1 behavior supplies
the conditions at infinity. This semilinear interface problem replaces
the linear resolvent equation obtained when optimizing only correlation
c; the full Gaussian center reward changes the forcing itself.

## 7. Threshold stationarity and the scope of the optimum

If alpha is optimized at an interior value, the envelope derivative at
fixed g supplies one more scalar condition. Let partial_alpha denote
differentiation of the coefficient formulas with c held fixed (but
p=2Phi(alpha)-1 varying), and put p'=2phi(alpha). Then

    0=p' Psi(lambda g(alpha)+gamma,t)
       +(partial_alpha lambda) E[H g s_c]
       +(partial_alpha gamma) E[H s_c]
       +B partial_alpha t+A p' g(alpha).          (19)

The last term comes from c'=p'g(alpha); the first is the moving-mask
boundary contribution. The Sobolev representative is continuous, so
g(alpha) is defined and this derivative is legitimate. Equations
(17)--(19), normalization, and the energy constraint describe the
interior scalar candidates. They are not a uniqueness theorem and do
not replace global certification by (12)--(14).

The maximizing scalar frame has a continuous inverse g_*. Appending
W=UH(V) is a measurable causal addition. Therefore, if its value exceeds
phi(0), the continuous-core theorem supplies a further full-response
improvement outside this scalar subfamily. Thus an attained C_scalar
above phi(0) is a class ceiling which is itself strictly improvable in
the larger finite response family. It is not a barrier for actual
signings. Conversely every value within the class is an actual signing
lower certificate by (5), with finite canonical approximation before
the matrix-dimension limit. The independently established ceiling for
the applicable fixed-depth algorithm family still applies; compactness
and stationarity here do not bypass it or prove a larger asymptotic
matrix value.

## 8. Uniform full-response escape from the entire scalar class

The attained scalar ceiling is separated strictly from the larger
response family, uniformly over all sufficiently good scalar states.
Fix c_0>phi(0) below C_scalar. The set of scalar parameters satisfying
T_alpha(g)>=c_0 is compact as in Section 3. On this entire set, not
merely on one fixed frame, the full unrestricted gradient gap has a
strictly positive minimum

    delta_scalar(c_0)>0.                          (20)

Here is an exact low-information formula for this gap. Regard all
functions below as functions of the scalar standard Gaussian v, and put

    K(v)=lambda g(v)+gamma H(v),
    B=E[H phi(K/sqrt(t))]/sqrt(t),
    w(v)=H(v)[2Phi(K(v)/sqrt(t))-1]-2BK(v),
    A_1=<w,g>,   nu=||w-A_1 g||_2.

The full gradient U w has covariance A_1 with the original fixed-point
Gaussian V and total variance ||w||_2^2. Conditional on V=v it is
distributed as A_1v+nu Z. This conditional identity concerns a Gaussian
first-chaos gradient, not an unproved law of an actual matrix output.
Since K depends only on V, the full gap is exactly

    gap=E_(V,Z) max{Psi(K(V),t), |A_1V+nu Z|+B}
          +qB-2T_alpha(g).                        (21)

Indeed <K,w>=T-2qB, using ||K||_2^2+t=q, so all terms involving the
old response F reduce to the displayed scalar moments. The Z integral
in (21) has a closed Gaussian-tail expression; only one outer Gaussian
integral remains. For b>=0,

    E(|m+nu Z|-b)_+
      =nu[phi((b-m)/nu)+phi((b+m)/nu)]
         +(m-b)Phi((m-b)/nu)
         +(-m-b)Phi((-m-b)/nu),

with its continuous zero-variance interpretation. For b<0 the same
expectation is Psi(m,nu^2)-b.

Every function g in S is continuous. The actual two-coordinate frame
consists of the continuous scalar core V and the causal addition
Z=(W-cV)/d. Its response F is genuinely ternary. The continuous-core
stationarity theorem therefore makes (21) strictly positive whenever
T>phi(0).

The gap is continuous across the compact scalar parameter set. The
inverse pairs (g,(H-cg)/d) vary strongly in L2 because d is uniformly
positive; the residual t is uniformly positive as well. Therefore K,B,w,
A_1 and nu vary continuously in their L2/scalar norms. The Lipschitz
absolute-value/positive-part bounds give continuity of (21), without
requiring nu>0. Positivity and compactness prove (20).

A uniformly bounded Hessian for a subsequent half-line toward the full
linearized optimum now gives a uniform actual objective gain. With M_0
chosen using the scalar class's positive t floor, set

    theta=min(1/2,delta_scalar/M_0),
    epsilon_scalar=theta delta_scalar/2>0.

If C_full denotes the supremum of the broader feasible full-response
family, then

    liminf_n M_n/n^(3/2)>=C_full
        >=C_scalar+epsilon_scalar.                (22)

The archived scalar polynomial certificate already exceeds 0.426
(`fresh_tree_energy_independent_audit_2026_09_05.md`, exact certificate
section), so c_0>phi(0) can indeed be chosen. Equation (22) is a
separation theorem, not a new evaluated decimal bound. Its positive
constant can in principle be certified by finite-Hermite approximation
and compact subdivision of the continuous gap formula (21). It is not
uniform over arbitrary growing Gaussian response frames.

## 9. Uniform actual fixed-point recovery and finite canonical banks

The director found a uniform critical height estimate, independently
checked here. Write g=sum_(d even) b_d h_d. Its scalar kernel obeys,
for every -1<=r<=1,

    K_g(r)=sum_d b_d^2 r^d >=(1+r^2)/2.            (23)

Indeed for even d>=2,
(1-r^d)/d <=(1-r^2)/2, and D(g)<=1. For the forced height covariance
p_(h+1)=K_g(p_h), p_0=0, the deficit delta_h=1-p_h therefore satisfies

    delta_(h+1)<=delta_h-delta_h^2/2,
    delta_h<=2/(h+2).                             (24)

The second inequality follows by induction from delta_0=1; the map
x-x^2/2 is increasing on [0,1], and
2/(h+2)-2/(h+2)^2<=2/(h+3). This bound holds uniformly over the ENTIRE
Sobolev class, including critical inverse functions. Thus

    ||V_g-Pi_h V_g||_2<=sqrt(2/(h+2)).             (25)

There is also an explicit global continuity estimate for the recovered
same-space frames. Let g,h belong to S, epsilon=||g-h||_2, and
r=E[V_g V_h]. The two actual creation equations give

    r=K_(g,h)(r)=sum_(d even) g_d h_d r^d.

Polarization and (23) imply

    K_(g,h)(r)
      =(K_g(r)+K_h(r))/2
         -(1/2)sum_(d even)(g_d-h_d)^2 r^d
      >=(1+r^2)/2-epsilon^2/2.

Every r^d here is nonnegative and at most one, even when r<0, so no
orientation or positive-correlation assumption is hidden. Rearranging
the fixed-point identity gives

    (1-r)^2<=epsilon^2,
    ||V_g-V_h||_2<=sqrt(2||g-h||_2).               (26)

This global one-half-Hoelder recovery bound shows directly that the
compact Sobolev class maps continuously into actual canonical first
chaos. Combining (26) with normalized degree-N truncation and (25),

    ||V_g-Pi_h V_(g_N)||_2
       <=2/(N+2)^(1/4)+sqrt(2/(h+2)).              (27)

For fixed N and h, every variable Pi_h V_(g_N) uses the SAME finite
canonical tree bank: height at most h and even offspring count at
most N. The number of those trees may be huge, but it is finite and
does not depend on g. Thus (27) gives uniform finite-coordinate
approximation for the whole scalar class, not just for a selected
maximizer. Subsequent bounded mask and response approximation still
precedes the matrix-dimension limit as required by the response theorem.
Neither a six-hour numerical search nor an unbounded-depth algorithm
at a fixed matrix dimension is being substituted for those limits.

## 10. High scalar extrema necessarily use the full derivative budget

The multiplier eta in (17) is strictly positive at EVERY local maximum
with c!=0. Consequently D(g)=1 there. Here local maximality refers to
the scalar Sobolev feasible set, and alpha can be held fixed in this
argument.

Suppose eta=0. If kappa!=0, the Euler equation forces g=0 outside the
central interval. Inside, g takes values in the zero set of the scalar
real-analytic function

    u -> kappa u-lambda s(lambda u+gamma)-A.

Since c!=0, lambda!=0, and this analytic function is not identically
zero: its derivative is kappa minus a nonconstant Gaussian density
multiple. Its real zero set is discrete. Continuity of g on the
connected central interval makes g constant there. Continuity at the
endpoints then makes that constant zero, contradicting ||g||=1.

If instead kappa=0, the full L2 gradient in (16) is zero. The strict
monotonicity of u -> lambda s(lambda u+gamma) again makes g constant
on the center. Choose a nonzero smooth even function delta supported
strictly inside the center with Gaussian mean zero. Then

    <H,delta>=0, <g,delta>=0, <g',delta'>=0.

Along the unnormalized straight perturbation g+theta delta, the moment
c, and hence lambda,gamma,t, stay fixed. Its objective has strictly
positive second variation

    lambda^2 E[H Psi_kk(lambda g+gamma,t) delta^2]>0.

The norm and energy can be corrected at order theta^2 without changing
this second variation, because the entire first derivative is zero.
Explicitly put e=b_0-b_0^2g, where b_0=Eg, and use

    g_theta=(g+theta delta+L theta^2 e)
                /||g+theta delta+L theta^2 e||_2.

If D(g)<1, choose L=0. If D(g)=1, its energy expansion is

    D(g_theta)=1+theta^2[D(delta)-||delta||_2^2
                                      -2L b_0^2]+O(theta^3).

By (2), b_0!=0, so a sufficiently large fixed L makes this strictly
feasible for all sufficiently small nonzero theta. The objective still
has the strictly positive displayed second variation. This contradicts
local maximality and proves eta>0.

The case c=0 cannot contain the global high-value scalar optimum. To
see this with an elementary exact bound, c=0 makes g orthogonal to H,
so (Eg)^2<=1-p. Equation (2) forces p<=1/2. Put C=2/pi<2/3. Formulas
(3)--(4) at c=0 give gamma=sqrt(C)q/sqrt(p) and t=q-Cq^2. Jensen then
gives

    T^2<=p^2(gamma^2+t)=p^2q+C p q^3.

The right side increases on 0<=p<=1/2. Its derivative is
2p-3p^2+C(1-p)^2(1-4p), manifestly positive on [0,1/4]. On [1/4,1/2],
using C<2/3 bounds the derivative below by

    [2-6p+9p^2-8p^3]/3>=1/12>0.

The polynomial in brackets decreases on that interval and equals 1/4
at 1/2. Hence

    T^2<=(2+C)/16<1/6,
    T<1/sqrt(6)<0.409.                            (28)

The archived actual scalar certificate exceeds 0.426. Therefore every
global scalar maximizer in (8) has c!=0, eta>0, and D(g)=1. Its inverse
is in the number-operator domain and obeys the genuine semilinear
interface equation (18); a strictly subcritical polynomial is an
approximation to this critical optimum, not an exact global maximizer.
