# Full-center exact gain and purification of the variational problem

Date: 2026-09-06. Status: mathematical derivations and exact rational replay
complete; the director independently reconstructed and replayed the endpoint
in `resumed_full_center_numerical_director_audit_2026_09_06.md`.
It depends on the independently reconstructed full
nonlinear center theorem, not merely on the earlier restricted edge
channel. No claim of convergence or iterative closure is made.

The full theorem and this track's independent audit are in
`resumed_response_full_channel_independent_audit_2026_09_06.md`.

## 1. Signed weighted Jensen improves the simpler norm bound

Write

    Gamma_tau(k)=E_N |k+tau N|,
    g_tau(x)=Gamma_tau(x)-x
             =2[tau phi(x/tau)-x Phi(-x/tau)],  x>=0.

For every feasible old pair F,H, put K=U*F, tau=||F-P_1F||_2,
mu=EH, and J=EHK. The full theorem gives the lower bound
E H Gamma_tau(K). If mu>0, convexity of Gamma_tau on ALL real k gives

    E H Gamma_tau(K)
       >= mu Gamma_tau(EHK/mu)
       = J+mu g_tau(J/mu),                         (1)

when J>=0. This is simply Jensen under the probability density H/mu.
It is stronger than first replacing K by |K| and then bounding its
weighted absolute mean by its L2 norm. No distributional information
about the inverse creation polynomial K is required.

## 2. Exact two-dimensional formulas for the canonical mask

Take the canonical 21-anchor degree-200 construction with alpha=361/500
and resolvent parameter 17/5. Its Gaussian V,W satisfy

    EV^2=1,  EW^2=p,  EVW=c,
    H=1_{|V|<=alpha},  F=sign(W)1_{|V|>alpha},
    p=EH=2Phi(alpha)-1.

The earlier finite certificate reconstructs the exact p,c values and
checks the strict creation contraction. Put d=sqrt(p-c^2)>0 and
write W=cV+dN_1 with V,N_1 independent standard Gaussians. Since F
depends only on these two first-chaos directions, its full first
Gaussian projection is aV+bN_1, where

    a=2phi(alpha)[2Phi(c alpha/d)-1]
        +4c/sqrt(p) phi(0) Phi(-alpha sqrt(p)/d),
    b=4d/sqrt(p) phi(0) Phi(-alpha sqrt(p)/d).       (2)

For completeness,

    a=2 integral_alpha^infinity
                 v phi(v)[2Phi(cv/d)-1] dv,
    b=4 integral_alpha^infinity phi(v)phi(cv/d) dv.

Integrating the first expression by parts gives its boundary term
and (c/d) times the second; the second Gaussian integral gives (2).
These formulas retain the exact signs and do not replace a by its
coarser upper bound 2phi(alpha).

Orthogonal projection and F^2=1-H now give

    sigma^2=||P_1F||_2^2=a^2+b^2,
    tau^2=||F-P_1F||_2^2=1-p-a^2-b^2,
    J=E F W=ca+db=E H U*F.                        (3)

Thus (1) is entirely computable from p,c and elementary one-dimensional
Gaussian functions, despite the complicated distribution of U*F.

## 3. Exact rational result

The script
`computations/resumed_response_full_center_certificate_2026_09_06.py`
reconstructs the canonical finite coefficient sums directly, using the
audited outward Fraction interval class, Machin pi enclosure, finite
exponential/Gaussian series, and explicit tails. It does not invoke the
canonical script's file-writing main, alter old certificates, or use
floating point, optimization, or numerical quadrature. Its output is
`computations/results/resumed_response_full_center_certificate_2026_09_06.json`.

The decisive intervals are

    a = .6143777828536417947547684299175914845799620478288993982286...
    b = .0017113479929010940154837256928531638842582718755096041833...
    tau^2 = .0928314887467733136762183059015820433517977600035876826565...
    tau = .3046826032886901038338515340995397481304193380589933961978...
    J = .4306581794055286027240538046347110263271732383363251904555...

The exact rational LOWER endpoints are

    p g_tau(J/p)
      > .000379279323472719152390980507389095679106795070238950997741,

    J+p g_tau(J/p)
      > .431037458729001321876444785142100122006280033406564141453321.

The latter is strictly greater than 431/1000, which the script checks
as an exact rational assertion. After independent theorem reconstruction
and numerical replay, a safely rounded original-problem lower endpoint is

    liminf M_n/n^(3/2) >= .4310374587290013.

An independent floating-point calculation using the direct Gaussian
formulas agrees; it is a diagnostic, not part of the certificate.

## 4. Purification: the full variational supremum may use ternary F

Define the full functional on bounded odd F and even feasible H by

    C(F,H)=E H Gamma_{||F-P_1F||_2}(U*F),
    0<=H<=1,  |F|+H<=1.

Since Gamma is nonnegative, first increase H to 1-|F|. This only
increases the functional. This step concerns the FULL feasible-pair
class; it need not preserve the more restrictive original relationship
F=sign(UH)(1-H).

Suppose F depends on a fixed finite ancestor-closed old Gaussian
family. Choose a tree coordinate G_* outside this family, adjoining its
finitely many ancestors if needed. It is independent of the old
coordinates on which F and K=U*F depend. Define

    U_*=2Phi(|G_*|)-1,
    F'=sign(F) 1_{U_*<|F|},
    H'=1-|F'|.

U_* is uniform on (0,1) and independent of the original finite
family. The gate is even in G_*, so F' remains jointly odd and H'
jointly even. They are pointwise feasible and take values
F' in {0,+1,-1}, H' in {0,1}. Conditional on the original fields,

    E[F'|old]=F,
    E[|F'||old]=|F|,
    E[H'|old]=1-|F|.                              (4)

Every old first-chaos coefficient is preserved by (4). The new G_*
coefficient is zero because its gate is even; all other newly adjoined
coordinate coefficients vanish by independence. Therefore

    P_1F'=P_1F,  U*F'=U*F=K,
    tau'^2=E|F|-||P_1F||_2^2
             >= EF^2-||P_1F||_2^2=tau^2.

The inverse creation K is still measurable with respect to the
original ancestor-closed family. Conditioning (4) and using that
Gamma_tau(k) is nondecreasing in tau give

    C(F',H')
      =E(1-|F|) Gamma_tau'(K)
      >=E(1-|F|) Gamma_tau(K)
      >=C(F,H).                                  (5)

If the original F is Gaussian-a.e.-continuous, so is F': conditional
on the old family, the threshold equality U_*=|F| has probability
zero. At F=0 the gate is eventually closed in a neighborhood whenever
U_*>0. Hence the finite bounded-response theorem applies directly.

For a general infinite-coordinate F there need not be an unused old
coordinate. We do NOT silently append one while claiming to preserve
the given creation operator. Instead first approximate F,H jointly
by finite ancestor-closed smooth feasible pairs. The full functional
is L2-continuous by the preceding audit. Purify each finite approximant
using its own unused coordinate. Consequently

    sup_all feasible countable (F,H) C(F,H)
      = sup_finite ternary F, H=1-|F| C(F,H).       (6)

This is equality of variational SUPREMA, not a direct purification of
every infinite-coordinate response while keeping all of its old data.

## 5. What purification does not give

The gate is an even function of an unused Gaussian tree coordinate,
not a new product-input spin vector after applying a matrix response.
It does not solve the twin-dependence obstruction and does not yield
an iterative fresh-center theorem. Its use is static: it reduces the
full variational optimization to ternary odd responses with binary
center masks, while preserving the first-chaos adjoint and increasing
the usable nonlinear variance.
