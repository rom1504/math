# An actual fixed-frame variational limit, stationarity, and strict non-sharpness

Date: 2026-09-06. Status: independently reconstructed by the director,
including the zero-variance boundary, uniform partition errors, restricted
stationarity, and the containing-frame escape implication. The two-variable
optimized rectangle certificate is separate.
This theorem describes what fixed-frame policy ascent optimizes toward,
without asserting that its local numerical optimizer finds the global
optimum or that the original sequence M_n/n^(3/2) converges.

## 1. The actual finite-frame functional

Let X=(X_1,...,X_d) be any fixed actual standard orthonormal Gaussian
frame in the canonical first chaos. Its inverse features h_j=U*X_j are
orthonormal even L2 functions on the full canonical space; they need NOT
be measurable in X. Put

    phi_j(X)=E[h_j|X], Gamma=E[phi phi^T]<=I_d.

For odd f and even q, measurable in X with |f|<=q<=1, define

    p=E q, a=E[X f], t=p-|a|^2,
    Kbar=a.phi,
    J_X(f,q)=E[(1-q)Psi(Kbar,t)].                   (1)

The norm in t is the TRUE inverse norm: ||sum a_j h_j||^2=|a|^2.
It is not a^T Gamma a. Finite canonical approximation, then even-gate
conditional ternary realization, the full-response theorem and Jensen
conditional on X make (1) an actual original-signing lower functional.
For a completely arbitrary finite X its inverse features may depend on
the entire countable Gaussian space; one must approximate first rather
than assert an exact unused gate independent of all those inverse
features. The present rich core has a finite closed containing frame,
so exact independent gates are available there. A value attained
in (1) is a guaranteed lower value, not a statement that the actual
output energy equals this conditional lower estimate.

For the present rich-core application, d=2, X=(V,Z),

    phi_1=G(V,Z):=E[g|V,Z],
    phi_2=[f_birth(V)-A G(V,Z)]/nu.

The actual inverse residuals are both multiples of g-G. Therefore

    I_2-Gamma=delta vv^T,
    delta=||g-G||^2, v=(1,-A/nu).

This rank-one projection loss explains why replacing the true inverse
norm by the conditional norm is not licensed. The polynomial G and the
piecewise-polynomial f_birth are known exactly. The numerical cell
ascent changed only the policy; these actual inverse features remained
fixed throughout.

## 2. A compact global optimum exists

The feasible pair domain is weak-* compact in L-infinity(X)^2. Parity
and |f|<=q<=1 are weak-* closed convex constraints. The moments p and a
are continuous because 1 and X_j are integrable. Thus t is continuous.
The Gaussian moment body gives

    |a|<=m(p):=2varphi(Phi^(-1)(1-p/2)),            (2)

where varphi denotes the standard Gaussian density. It follows that
t>=0, and t>0 whenever p>0. Indeed m(p)^2<p for every p>0, by strict
Cauchy--Schwarz in the extremizing two-sided Gaussian tail.

For clarity, the elementary moment-body argument is
|a|=E[(e.X)f]<=E[|e.X|q] for e=a/|a|. Among masks of mass p, the last
integral is maximized by the two Gaussian tails. This proves (2) without
any smoothness assumption on the policy.

The bound

    |Psi(k,t)-Psi(k',t')|
       <=|k-k'|+sqrt(2/pi)*|sqrt(t)-sqrt(t')|       (3)

shows strong L1 continuity of the Psi integrand as (p,a) varies,
including t=0. Weak-* convergence of q then gives continuity of J_X.
Therefore

    C_X=max_(|f|<=q<=1, parity) J_X(f,q)            (4)

exists. No positive t margin is assumed on the WHOLE domain. The only
t=0 boundary point has p=0 and score zero.

On a superlevel J_X>=c>0, there is a uniform positive t margin. Since
Gamma<=I, Cauchy--Schwarz gives J_X<=sqrt(p), so p>=c^2. Hence

    t>=min_(p in [c^2,1]) [p-m(p)^2]=:t_c>0.       (5)

The stronger elementary bound J_X<=sqrt(p(1-p)) is also valid: apply
Cauchy--Schwarz to (1-q)|Kbar+sqrt(t)N| and use 0<=1-q<=1.

## 3. Exact finite-partition approximation, uniformly over all policies

Let P be a finite partition invariant under X->-X. Write

    X_P=E[X|P], phi_P=E[phi|P],
    epsilon_X=||X-X_P||_(L2 vector),
    epsilon_phi=||phi-phi_P||_(L2 vector).

For any feasible pair set f_P=E[f|P], q_P=E[q|P]. Feasibility and parity
are preserved, p_P=p, and

    |a-a_P|=|E[(X-X_P)f]|<=epsilon_X,
    |t-t_P|<=2epsilon_X.                           (6)

Here |a|,|a_P|<=1. Because Gamma<=I, the mean-field difference satisfies

    ||a.phi-a_P.phi_P||_1<=epsilon_X+epsilon_phi.   (7)

Define the finite-cell functional

    J_P(f_P,q_P)=E[(1-q_P)Psi(a_P.phi_P,t_P)].      (8)

Its integrand is P-measurable, so q_P can be replaced by q in this
expectation. Equations (3), (6), and (7) yield the GLOBAL uniform bound

    |J_X(f,q)-J_P(f_P,q_P)|
       <=epsilon_phi+epsilon_X
                          +sqrt(2/pi)*sqrt(2epsilon_X).           (9)

This bound includes p=0 and arbitrarily small t. On a common region
t,t_P>=t_min>0, the sharper bound is

    |J_X-J_P|<=epsilon_phi
                   +(1+2varphi(0)/sqrt(t_min))*epsilon_X.        (10)

To use (10) uniformly on J_X>=c, first make the right side of (9)
smaller than c/2. Then J_P>=c/2, and the moment body supplies a common
positive t margin. One must NOT silently apply (10) to the whole domain.

Let C_P be the global maximum of (8) over its finite cell feasible
polytope. Conditional Jensen gives J_P(f_P,q_P)<=J_X(f_P,q_P), hence

    C_P<=C_X<=C_P+the right side of (9).            (11)

For nested symmetric partitions with dense union, both epsilons tend
to zero. Further refinement can only increase C_P: retain a coarse
pair and apply Jensen to the finer conditional feature vector. Thus

    C_P increases to C_X.                          (12)

This is a genuine finite-class variational limit theorem. The finite
optimization is generally nonconvex; (12) concerns its GLOBAL optimum,
not the output of a local policy iteration. Exhaustive finite-dimensional
optimization with interval bounds is possible in principle and is not
claimed efficient.

In the current core every required cell coefficient is a Gaussian
polynomial/rectangle moment. The scalar principal-angle theorem bounds
the omitted Z-Hermite levels. In fact the vector inverse projection
error from discarding levels >=K is at most

    sqrt(1+(A/nu)^2)*||b_core||^K.

Consequently the uniform approximations in (11) can be bounded with
finite arithmetic. For actual signings, choose the partition and
Hermite cutoff first, approximate the fixed actual core coordinates by
finite canonical tree coordinates next, and only then take matrix size
to infinity. Equation (12) does not run a growing mesh or an infinite
response depth at a fixed matrix size.

## 4. Exact stationary threshold within the fixed frame

At p>0 set H=1-q, k=a.phi, and

    s(k)=2Phi(k/sqrt(t))-1,
    B=E[H varphi(k/sqrt(t))]/sqrt(t),
    A_j=E[H s(k)phi_j]-2B a_j.

The exact derivative is

    dJ_X=E[(A.X)delta f+(B-Psi(k,t))delta q].       (13)

Therefore every local maximum obeys the pointwise threshold rule

    f=sign(A.X)q when A.X!=0,
    q=1 when |A.X|+B>Psi(k,t),
    q=0 when |A.X|+B<Psi(k,t),                    (14)

with ties allowed. These are necessary conditions, not sufficient ones.
The vector A is the gradient projected to the fixed Gaussian frame;
it is not the full creation-space gradient. In particular A=0 does NOT
mean that the full gradient is zero, and the universal full-gradient
zero-state cap must not be imported into this restricted problem.

For finite cells, replace X and phi in (13)--(14) by their cell means.
This is precisely the oracle used in the saved same-frame ascent code.
Its exact line search preserves the conditional feasible polytope and
increases the stated cell objective. A vanishing cell oracle gap proves
only first-order stationarity on that cell polytope; it is neither a
global upper certificate nor unrestricted response stationarity.

In the present rich core, phi is polynomial separately on |V|<alpha
and |V|>alpha. If A!=0, the tie set in (14) has Gaussian measure zero.
Indeed on each V-piece and each sign halfspace of A.X, the tie equation
is real analytic. An identity on any open set would analytically extend
Psi(k,t)=+/-A.X+B to all of R^2 using that piece's polynomial extension,
contradicting positivity of Psi where the affine right side is negative.
Thus nonzero-A stationary policies here are genuinely ternary threshold
policies. No such conclusion is asserted on the A=0 branch.

## 5. Why this fixed class cannot be globally sharp for response theory

There is a separate STRICT gap beyond (4) for the actual rich core.
The full core E=(21 anchors,Z_*) has polynomial inverse features and is
therefore an analytic cyclic core. Although Z is not independent of E,
write

    Z=b_core.E+sqrt(1-||b_core||^2)T.

The certificate proves ||b_core||^2<.31. Thus T is an actual standard
Gaussian independent of E, and

    U*T=[h_Z-sum_j b_core,j h_(E_j)]/sqrt(1-||b_core||^2)

is measurable in E. The enlarged standard frame (E,T) is exactly a
fixed analytic cyclic core plus one arbitrary measurable causal
appendage. It contains (V,Z).

The independently audited uniform fixed-frame escape theorem in
`resumed_response_analytic_core_causal_extension_stationarity_2026_09_06.md`
therefore applies to EVERY feasible pair in this fixed frame above
varphi(0). Let C_frame be its full, unconditioned response maximum.
Its compactness proof gives attainment, and the uniform escape theorem
gives some epsilon_frame>0 with

    C_X<=C_frame<C_frame+epsilon_frame<=C_full.     (15)

The present lower certificate is above varphi(0), so the relevant
superlevel is nonempty. Equation (15) proves that neither the optimized
two-variable policy class nor even the entire containing fixed frame
can be sharp for the unrestricted full-response functional.

This epsilon is not evaluated and is not uniform over an expanding
sequence of frames. The conclusion does not compel ascent to the
terminal functional's universal ceiling, much less to 1/2. A finite
class can have a perfectly genuine stationary optimum while still
omitting a useful creation-space direction. The new global terminal-
functional ceiling, proved separately, shows why indefinitely enriching
this same terminal objective cannot by itself solve the original
1/2/convergence question.
