# Causal inverse-shape optimization above a fixed rich Gaussian core

Date: 2026-09-06. Status: exact mechanism and fixed-covariance dual independently
reconstructed. The director read and replayed the separate rational certificate,
proving liminf M_n/n^(3/2)>=.4314713871136104. Numerical shape-search outputs
remain diagnostics; they are not certificates of a global shape optimum.

## 1. The one-dimensional conditional inverse is sufficient information

Let V be an already actual standard first-chaos Gaussian, with inverse
g=U*V of norm one. The inverse g may depend on a rich finite creation-closed
Gaussian frame, including a cyclic core. Assume the conditional projection

    m(v)=E[g|V=v]

is known. In the current 21-anchor construction it is the exact even
degree-200 polynomial saved by the conditional-V certificate.

Choose ANY even u in L2 of the scalar Gaussian V with ||u||_2=1. Set
W_u=U u(V). This is a causal addition: u is a function of the already
constructed V, and is not a function of its own newly created coordinate.
There is NO derivative-energy restriction on u. The pair (V,W_u) is
standard jointly Gaussian with covariance

    c=<g,u(V)>=<m,u>,  d=sqrt(1-c^2).

Fix a center threshold alpha>0, let H=1{|V|<=alpha}, p=EH, q=1-p,
and take F=sign(W_u)(1-H). For d>0 define

    lambda=2phi(alpha)[2Phi(c alpha/d)-1],
    gamma=4phi(0)Phi(-alpha/d),
    t=q-lambda^2-2c lambda gamma-gamma^2.

The exact first-chaos calculation gives

    U*F=lambda g+gamma u(V),   ||F-P1F||_2^2=t.

The full-center theorem, followed by conditioning only on V and Jensen,
therefore gives the ACTUAL signing lower bound

    liminf_n M_n/n^(3/2)
      >=E[H Psi(lambda m(V)+gamma u(V),t)].        (1)

This does not replace g by m inside its own fixed-point equation. In
particular ||m|| can be less than one and its scalar derivative energy
can exceed one. Neither fact affects the causal construction of W_u.
For the current polynomial m, ||m|| is about .99504, so all unit u have
|c|<=||m||<1 and there is a uniform nondegenerate output covariance.

## 2. A concrete exactly normalized three-shape candidate

For a chosen outer mask alpha define

    b=E[Hm],  s_2=E[Hm^2],  r^2=E m^2,
    v_1=s_2-b^2/p,  v_2=r^2-s_2.

When v_1,v_2>0, the three even functions

    e_0=H/sqrt(p),
    e_1=H(m-b/p)/sqrt(v_1),
    e_2=(1-H)m/sqrt(v_2)

are EXACTLY orthonormal. Their correlations with m are
b/sqrt(p), sqrt(v_1), sqrt(v_2), respectively. Thus

    u=x_0 e_0+x_1 e_1+x_2 e_2,  sum x_i^2=1,

is an actual unit causal inverse with explicitly known covariance c.
On the center, the conditional inverse in (1) is simply

    lambda m+gamma u=A m+B,
    A=lambda+gamma x_1/sqrt(v_1),
    B=gamma[x_0/sqrt(p)-x_1 b/(p sqrt(v_1))].       (2)

Consequently every norm, covariance, and bin mean needed for a lower
certificate is a one-dimensional polynomial moment. No high-dimensional
integration of the rich inverse g is needed.

The original rich core is fixed at its already certified parameters
alpha_core=3623/5000 and resolvent=3479/1000. The new candidate uses

    alpha=91/125,
    x_1=229/20000, x_2=-647/50000,
    x_0=sqrt(1-x_1^2-x_2^2).

The numerical search also tried an orthogonalized half-center indicator.
Its observed gain over these three shapes was only about 4e-8, so it
was omitted from the exact candidate. That omission does not assert
optimality of the three-dimensional shape bank.

## 3. Exact lower-dimensional dual at fixed covariance

Fix alpha and a covariance c with |c|<r=||m||_2. In this section lambda,
gamma,t are fixed scalars determined by that c. Put

    Phi_c(u)=E[H Psi(lambda m+gamma u,t)].

The fixed-covariance problem is

    P(c)=sup{Phi_c(u): u even, ||u||_2=1, <m,u>=c}.

The norm equality can be replaced by ||u||_2<=1 without changing the
supremum. Indeed outside H the objective is zero. Given any smaller-norm
u, choose an even tail-supported L2 direction orthogonal to both m and
u, then add exactly enough of it to fill the norm. The tail Gaussian
space is infinite dimensional, so such a direction exists.

For eta>0 and beta real define

    D_c(eta,beta)=eta+beta c
       +beta^2 ||(1-H)m||_2^2/(4eta)
       +E[H sup_(y real){Psi(lambda m+gamma y,t)
                                         -eta y^2-beta m y}].       (3)

Then

    P(c)=min_(eta>0,beta real) D_c(eta,beta).       (4)

Thus the moment dual has only TWO scalar multipliers and a pointwise
one-dimensional maximization, regardless of the degree or rich-core
dimension. It is an upper certificate for this fixed covariance slice,
not an upper bound on actual signings or on other response families.

### Weak duality and pointwise maximizers

Add eta(1-||u||^2)+beta(c-<m,u>) to the objective. Outside the center,
the pointwise supremum is beta^2m^2/(4eta), attained at
u=-beta m/(2eta). Inside, maximizing the coercive scalar function in
(3) gives the claimed upper bound. Every inside maximizer obeys

    2eta y+beta m=gamma[2Phi((lambda m+gamma y)/sqrt(t))-1],          (5)
    |y|<=(gamma+|beta m|)/(2eta).

On the compact center the polynomial m is bounded, so the pointwise
maximizers are uniformly bounded. Measurable grid approximations to the
maximizer show that the supremum of the integral equals the displayed
integral of suprema. No unbounded-action interchange is left unjustified.

### Why strong duality holds despite the nonconcave pointwise reward

An elementary atomless-mixing argument supplies the needed convexity.
For two even L2 functions u_1,u_2, select u_1 or u_2 on rapidly alternating
sets of the uniform variable 2Phi(|V|)-1. The indicators converge weak-*
to any desired constant mixing fraction. Testing against the three L1
differences of squared norm, covariance, and objective shows that the
closure of their finite-dimensional moment/objective image is convex.
The assertion follows first for continuous functions of the uniform
variable and then by L1 density; exact pointwise convexity of Psi is not
being assumed.

Taking closure does not increase the supremum at an interior covariance
|c|<r. The function u_0=c m/r^2 has the correct covariance and norm
strictly below one. Mix any approximately feasible sequence with a
vanishing fraction of u_0 to create an energy margin larger than its
moment errors. Then correct the remaining covariance error by adding a
vanishing multiple of m. Cauchy--Schwarz controls the energy change, and
the gamma-Lipschitz bound in u controls the objective change. This yields
exactly feasible functions approaching the same score.

More explicitly, for the objective with its parameters frozen at c,
the optimal closed-image value as a function of an energy budget e and
a covariance z is finite and concave on e>z^2/r^2. The point (1,c) is
interior. A finite supporting plane there gives multipliers eta>=0,beta,
and hence an attained dual value equal to the primal supremum. Monotonicity
in the energy budget gives eta>=0. In fact eta cannot be zero: on a
positive-measure central set, Psi(lambda m+gamma y,t)-beta m y has an
unbounded supremum in at least one sign of y, since gamma>0. Thus eta>0,
proving (4). At |c|=r the only unit feasible function is +/-m/r; the
interior dual-attainment claim is not asserted at those endpoints.

This proof identifies a supremum and its dual, not automatic attainment
by a single pointwise maximizing branch. When branches tie, one can
use atomless mixtures and feasible approximation as above. A selected
branch satisfying the two required moments is, of course, a genuine
attaining primal certificate.

## 4. A particularly tractable unique-branch regime

If

    eta>gamma^2 phi(0)/sqrt(t),                    (6)

the scalar maximization in (3) is strictly concave: its second derivative
is at most 2gamma^2phi(0)/sqrt(t)-2eta<0. Equation (5) then has a unique
solution for every m value, and its fixed-point form is a contraction.
The outside solution remains linear in m. In this regime, the norm and
covariance equations are just two one-dimensional Gaussian integrals
in eta,beta. If they equal one and c, the resulting measurable u
simultaneously attains (3) and the primal problem, giving a GLOBAL
fixed-covariance optimum, not merely a stationary shape.

The full optimization must still compare different c and outer thresholds
alpha. Neither (4) nor a successful local shape search supplies that
comparison automatically. The exact norm and moment dual nevertheless
provides an independently checkable way to distinguish a true shape
improvement from a spurious unconstrained inverse or a self-feedback
construction that the canonical space cannot realize.

## 5. Exact all-shape upper certificate on the selected covariance slice

The separately saved program
`computations/resumed_response_feedforward_fixed_c_dual_certificate_2026_09_06.py`
uses the exact covariance of Section 2 and rational multipliers

    eta=852051/400000000, beta=1259877/5000000000.

It rigorously upper-bounds (3), hence EVERY even unit shape u at this
same alpha and c, by

    .431472206961858089947095440100843388061138810230859493559171.

The exact three-shape lower certificate is

    .431471387113610406780573650606255124749081686370408082604706,

so the certified fixed-slice gap is less than 8.2e-7. Floating solution
of both moment equations gives primal and dual .431471950413875; that
decimal is diagnostic, not an additional certified bound.

Here is the entire upper-bound mechanism. The pointwise dual value
f(m)=sup_y{Psi(lambda*m+gamma*y,t)-eta*y^2-beta*m*y} is convex in m.
The negative second derivative in y is at least
kappa=2eta-2gamma^2phi(0)/sqrt(t)>.00415784817. At any rational test y0,

    f(m)<=objective(m,y0)+|partial_y objective(m,y0)|^2/(2kappa).

Four interval-CDF fixed-point iterations select a rational y0 at each
endpoint of an enclosed m range; this inequality, not convergence of
those iterations, proves the upper bound. No floating root solver enters
the certificate.

The center is partitioned into 128 positive-half bins. An exact midpoint
value and derivative of m, with a Taylor remainder, enclose its range on
each bin. The uniform second-derivative bound is

    |m''|<=sqrt(192 sum_k k(k-1)m_k^2)<62.81484260.

Indeed Mehler's formula at r=.99 gives sum_(j<=198)h_j(v)^2<192 for
|v|<=1: r^(-198)<8, (1-r^2)^(-1/2)<8, and
exp(r*v^2/(1+r))<3. Cauchy--Schwarz proves the displayed derivative
bound. Convexity bounds f by the chord through the two upper endpoint
values. Each chord integral is then enclosed using exact Gaussian bin
mass and the Hermite-endpoint integral of m. The tail contribution is
exactly beta^2*v_2/(4eta).

All operations use rational outward intervals. The JSON field
`dual_upper_enclosure` encloses this constructed chord upper bound;
its lower endpoint is NOT a lower bound on the optimal dual or primal
value. Only its upper endpoint is used to bound P(c). The certificate
does not optimize alpha or c and does not bound any larger family of
actual response constructions.
