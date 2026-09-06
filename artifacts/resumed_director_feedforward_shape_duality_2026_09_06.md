# Exact two-parameter dual for a causal feedforward response shape

Date: 2026-09-06. Status: proof reconstructed by the director and independently
verified in `resumed_bound_audit_feedforward_shape_dual_2026_09_06.md`.
This theorem optimizes one declared response class. It is
not a formula for the signing optimum or the full Gaussian response class.

## 1. The feasible class and the scalar objective

Work on a nonatomic probability space; for Gaussian even functions use
the nonatomic space of absolute Gaussian values. Let H be an indicator
of mass p>0, let h be real L2 with rho=||h||2^2>0, and fix

    |c|<sqrt(rho), gamma>0, t>0, lambda real.

Among real functions u satisfying ||u||2=1 and <h,u>=c, consider

    J*=sup_u E H Psi(lambda h+gamma u,t),
    Psi(k,t)=E|k+sqrt(t)Z|, Z standard normal.                 (1)

In the actual response application h=E[g|V], where V=Ug is already an
actual canonical creation coordinate. Any such u(V) is an admissible
feedforward inverse: its output Uu(V) is Gaussian and has covariance c
with V. No scalar derivative-energy restriction is imposed on u, and no
new self-feedback fixed point is being assumed. Conditional Jensen gives
(1) as a lower certificate for that actual two-coordinate construction.
The constants lambda,gamma,t depend on the chosen mask and covariance c;
they are held fixed throughout the shape optimization below.

Define

    A=lambda+gamma*c/rho,
    d=gamma*sqrt(1-c^2/rho)>0,
    C_t(s)=2sqrt(t) phi(Phi^(-1)((1+s)/2)), -1<=s<=1,

with C_t(-1)=C_t(1)=0. The exact theorem is

    J*=inf_(eta>0,b real) [ d*eta/2+eta*rho*b^2/(2d)
       + E H max_(-1<=s<=1)
             { (A-b)h*s+C_t(s)+d*s^2/(2eta) } ].              (2)

Thus an infinite-dimensional feedforward shape search has a two-scalar
global dual and a one-dimensional pointwise maximization. Multiple
pointwise maximizing branches are allowed. Ignoring those branches can
give a false optimization claim.

## 2. Gaussian absolute-value conjugacy and the Hilbert sphere

Differentiating Psi gives Psi_k(k,t)=2Phi(k/sqrt(t))-1. Hence direct
Legendre inversion, including the endpoint limits, proves

    Psi(k,t)=max_|s|<=1 {sk+C_t(s)}.                           (3)

For measurable s on H, write

    m=E H h s, z=E H s^2, e=E H C_t(s), v=z-m^2/rho>=0.

The two suprema over u and s may be interchanged: both are suprema, not
a minimax exchange. The support function of the Hilbert sphere slice is

    sup_(||u||=1,<h,u>=c) <Hs,u>
      =c*m/rho+sqrt(1-c^2/rho)*sqrt(z-m^2/rho).                (4)

When v>0 it is attained at

    u=(c/rho)h+sqrt(1-c^2/rho)*(Hs-(m/rho)h)/sqrt(v).          (5)

For v=0 use any orthogonal unit direction; this exists in the nonatomic
L2 space. Therefore

    J*=sup_s { A*m+e+d*sqrt(v) }.                             (6)

Equation (5) also constructs a genuine feedforward inverse from a
near-optimal dual/primal witness. No full Boolean landscape is encoded.

## 3. Why scalar moments can be convexified without a response loss

Let K be the closure of the triples (m,z,e) arising from deterministic
measurable s in [-1,1]. Then K is compact and convex. A direct finite
approximation proof suffices: truncate/approximate h in L1 by a simple
function, approximate the compact control interval by a finite mesh
(C_t is uniformly continuous), and partition each resulting nonatomic
cell into arbitrary prescribed proportions. This implements any finite
randomized mixture of controls to arbitrary error in all three moments.
The same approximation shows K is the set of integrals of randomized
control kernels, or equivalently their compact closed moment set.
The m-coordinate is uniformly integrable since |Hs|<=1 and h is L1;
the other two coordinates are bounded. This also proves closedness by
the finite-cell approximation or by diagonal weak compactness of the
randomized kernels. No pointwise concavity in s is asserted.

Consequently the continuous concave function

    F(m,z,e)=A*m+e+d*sqrt(z-m^2/rho)

attains its maximum on K, and deterministic controls approach that
maximum. Passing to K does not change the supremum in (6). A claim of a
particular deterministic maximizing shape requires additional selection
or purification; the theorem needs only approximating shapes.

Every maximizer has v>0. If v=0 and m!=0, mix its triple with the s=0
triple. At mixing fraction theta, the new v is
theta(1-theta)m^2/rho, whose square-root gain dominates the O(theta)
change in the linear terms. If v=m=0, then z=0 and the triple is s=0.
Take a randomized control equal to +/-delta with equal probabilities
on H. It has m=0,z=p delta^2 and
e=p C_t(delta)=p C_t(0)-O(delta^2). Its d sqrt(p)delta gain is positive
for small delta. The convexification argument supplies deterministic
approximations to the same strictly improved value. Thus the boundary
v=0 cannot maximize (6).

## 4. Supporting-plane proof of the exact dual

For eta>0,b real, completing a square gives, for every feasible triple,

    d sqrt(v) <= d eta/2 + d z/(2eta)
                  +eta*rho*b^2/(2d)-b*m.                    (7)

Equality holds precisely when eta=sqrt(v) and
b=d*m/(eta*rho). Equivalently first use
sqrt(v)<=eta/2+v/(2eta), then represent the negative quadratic
-d*m^2/(2eta*rho) as the infimum of
eta*rho*b^2/(2d)-b*m.

This proves that the right side of (2) is at least J*. Conversely let
(m*,z*,e*) maximize F over K. Its v*>0, so F is differentiable there.
Its supporting gradient maximizes the associated linear functional
on K: for every k in K, grad F(k*) dot (k-k*)<=0. Set

    eta*=sqrt(v*), b*=d*m*/(eta* rho).

The affine majorant in (7), with A*m+e added, is exactly this supporting
functional plus a constant, touches F at k*, and has maximum J* on K.
Maximizing that affine expression over control kernels, or deterministic
controls, separates pointwise and gives the integrand in (2).
This proves equality and actual attainment of the two-scalar infimum.
No joint convexity of the displayed (eta,b) parameterization and no
unjustified exchange of infimum and supremum is used.

For interior pointwise stationary controls the equation is

    (A-b)h - sqrt(t) Phi^(-1)((1+s)/2)+(d/eta)s=0.             (8)

All roots and endpoints must be compared. A certified search using (2)
must enclose this one-dimensional maximum, not merely solve (8) on one
locally selected branch.

## 5. Consequences and limits

This is a strict finite-dimensional optimization reduction for a new
causal inverse above a fixed actual scalar observation, at fixed output
covariance. The reduction can produce both rigorous upper bounds on
this declared shape class and constructive lower witnesses. The parent
Gaussian inverse may be high-dimensional; conditioning uses only h.

It does not optimize the original inverse g, the mask H, or the full
unconditioned response. Additional optimization over c and the mask is
still required, and conditional Jensen may be strict. It supplies no
cross-order inequality and no convergence conclusion by itself.
