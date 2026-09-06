# Independent audit: exact fixed-covariance feedforward shape dual

Date: 2026-09-06. Status: proof reconstructed independently; no open
step found under the hypotheses below. Audited files:

- `resumed_director_feedforward_shape_duality_2026_09_06.md`;
- Section 3 of `resumed_response_causal_inverse_shape_dual_2026_09_06.md`.

This is an optimization theorem for one fixed causal response class.
It is not an upper bound on the signing optimum.

## 1. Precise domain and reduction

Work on an atomless probability space. For scalar Gaussian even
functions, use the atomless quotient represented by |G|. Let H be an
INDICATOR with p=EH>0, h in L2, rho=E h^2>0, and assume

    gamma>0, t>0, |c|<sqrt(rho), lambda real.

All controls and h,H are measurable on the same space. Define

    Psi(y,t)=E|y+sqrt(t)N|,
    C_t(s)=2sqrt(t) phi(Phi^(-1)((1+s)/2)), |s|<=1,
    A=lambda+gamma c/rho,
    d=gamma sqrt(1-c^2/rho)>0.

At s=+/-1, set C_t(s)=0. Direct differentiation and endpoint limits give

    Psi(y,t)=max_|s|<=1 {sy+C_t(s)}.

For a measurable control s on H, put

    m=E H h s, z=E H s^2, e=E H C_t(s), v=z-m^2/rho.

Cauchy--Schwarz gives v>=0. Interchanging TWO suprema is legitimate,
and the support function of the sphere slice is exactly

    sup_(||u||=1,<h,u>=c) <Hs,u>
      =cm/rho+sqrt(1-c^2/rho)*sqrt(v).

For v>0, the maximizing u is

    u=(c/rho)h+sqrt(1-c^2/rho)*(Hs-(m/rho)h)/sqrt(v).

For v=0, any orthogonal unit direction completes the prescribed norm;
the atomless L2 space supplies one. Consequently the original supremum
is exactly

    J*=sup_s [A m+e+d sqrt(v)].                     (1)

The INDICATOR hypothesis matters: for a general weight H, the norm
moment would be E[H^2 s^2], not E[H s^2]. No such weighted extension is
silently used here.

## 2. The moment closure has no relaxation gap

Let K be the closure in R^3 of attainable triples (m,z,e). It is bounded,
hence compact. It is convex by atomless mixing: for two controls, the
three coordinate differences are integrable. Their vector-measure
range permits the prescribed common mixing fraction, or finite simple
approximations and arbitrarily fine atomless subdivisions give it in
the closure. Controls may remain even by doing this on |G|.

The function F(m,z,e)=A m+e+d sqrt(z-m^2/rho) is continuous on K and
concave on v>=0. By the DEFINITION of K, its maximum is the supremum
of actual deterministic controls in (1). Thus randomized moments add
no value. A deterministic maximizing control is not required or claimed.

Every maximum has v>0. If v=0 but m!=0, mix a putative maximizing triple
with the s=0 triple, retaining a fraction 1-theta of the former. The
new v is theta(1-theta)m^2/rho. Its positive sqrt(theta) gain dominates
the O(theta) linear loss for sufficiently small theta.

If v=m=0, then z=0 and necessarily e=p C_t(0). The last assertion also
holds in the closure: z->0 forces s->0 in H-weighted probability, and
C_t is bounded and continuous. A randomized equal choice +/-delta on
H has m=0, z=p delta^2, e=p C_t(delta). These moments belong to K by
the same atomless approximation. Since C_t(delta)=C_t(0)-O(delta^2),
the gain d sqrt(p)delta dominates the quadratic loss. This excludes the
remaining boundary. In particular the dual below attains eta>0; an
eta->0 limiting convention is unnecessary under these hypotheses.

## 3. Supporting-plane dual, with signs and factors checked

For eta>0 and b real, completing squares gives

    d sqrt(v) <= d eta/2+d z/(2eta)
                    +eta rho b^2/(2d)-b m.         (2)

The two nonnegative errors are

    d(eta-sqrt(v))^2/(2eta),
    eta rho [b-dm/(eta rho)]^2/(2d).

Thus equality holds at eta=sqrt(v), b=dm/(eta rho). At a maximizing
triple with v>0, the gradient of the concave F is a supporting linear
functional of K: its dot product with every feasible displacement is
nonpositive. Taking these touching eta,b in (2) therefore gives

    J*=min_(eta>0,b real) [d eta/2+eta rho b^2/(2d)
       +E H max_|s|<=1 {(A-b)h s+C_t(s)+d s^2/(2eta)}].       (3)

Pointwise separation is valid: a continuous maximization over compact
[-1,1] admits measurable grid approximations, bounded by a constant
plus |h|. There is no infimum/supremum swap hidden in this proof.

As a fallback that also handles related degenerate variants, replace
sqrt(v) by sqrt(v+epsilon). Its differentiable tangent gives the same
dual with an additional d epsilon/(2eta). Uniform convergence within
d sqrt(epsilon) then proves equality in the limit. Here the strict
interiority proof above yields the stronger attained minimum.

## 4. Exact equivalence with the energy/covariance dual

Let kappa>0 and beta be the response agent's multipliers. Taking the
Gaussian conjugate before maximizing the scalar u gives its dual as

    kappa+beta c+beta^2 rho/(4kappa)
      +E H max_|s|<=1 {
          [lambda-gamma beta/(2kappa)]h s
             +C_t(s)+gamma^2 s^2/(4kappa)}.         (4)

The bijective parameter change

    kappa=gamma^2 eta/(2d),
    beta=(gamma eta/d)(b-gamma c/rho)

turns (4) EXACTLY into (3): the s^2 coefficient becomes d/(2eta), the
linear coefficient becomes A-b, and the constant becomes
d eta/2+eta rho b^2/(2d). This independently checks both duals' factors
and signs. It also proves the energy dual's exact value without relying
on an unbounded-control minimax theorem.

Both duals optimize a fixed covariance c with lambda,gamma,t frozen.
Changing c requires recomputing those response parameters. Pointwise
maximizers can have multiple branches; solving only one stationarity
equation does not certify the maximum appearing in (3) or (4).
