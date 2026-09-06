# Independent mechanism audit: a rich-core response birth

Date: 2026-09-06. Status: structural proof reconstructed completely,
exact rectangle implementation audited, and independent rational replay
passed. The new certified lower endpoint is .4320510727290484005169....
The larger floating .43205335 diagnostic is still NOT a certificate.

Read in full:

- `resumed_response_rich_core_response_birth_2026_09_06.md`;
- `resumed_response_conditional_v_full_center_2026_09_06.md`;
- `fresh_finite_anchor_fixed_point_2026_09_05.md`;
- `computations/resumed_response_rich_core_birth_2026_09_06.py`;
- the conditional-V exact coefficient certificate implementation.

## 1. Actual creation variables, not a formal scalar replacement

The old core supplies an ACTUAL standard first-chaos Gaussian V=Ug,
where g is an even degree-200 polynomial in the independent finite
Gaussian frame E=(G_0,...,G_20,Z_*), with ||g||2=1. Write

    V=w.E, w=(rho,s), ||w||2=1, m(V)=E[g|V].

The old causal inverse u(V) has norm one and gives W=Uu. It is allowed
to be discontinuous at the old mask boundary; no self-feedback equation
for u is being solved. Let c=<g,u>=<m,u>. The old ternary response has
first-chaos inverse K0=lambda g+gamma u and support probability p0.

Choose any even scalar L2 function f(V), define

    L=Uf, A=<g,f>=<m,f>, nu^2=||f||2^2-A^2>0,
    Z=(L-A V)/nu.

Then V,Z are independent standard Gaussians, because they are actual
orthogonal first-chaos coordinates. The exact inverse of Z is

    U*Z=(f-A g)/nu.

The norm of g is ONE throughout; it is never replaced by ||m||<1.
This matters in nu and every subsequent first-chaos norm.

For an odd ternary response F1(V,Z), put p1=E F1^2,
aV=E[V F1], aZ=E[Z F1]. Gaussian independence shows that there are no
other first-chaos coefficients, even when the core has many more
coordinates. Therefore

    K1=(aV-aZ A/nu)g+(aZ/nu)f,
    ||K1||2^2=aV^2+aZ^2.

For the threshold diagnostic F1=sign(A V+nu Z) times its outer-threshold
indicator, the conditional endpoint formulas in the code are exact:
the conditional signing mean is Phi(-upper)-Phi(lower), and its first
Z moment is phi(upper)+phi(lower). The cutoff is an even function of V,
so F1 is jointly odd under (V,Z)->(-V,-Z).

## 2. The interpolated support variance is genuinely realizable

Take an independent standard Gaussian gate G_gate and an EVEN event
of probability theta. Choose F1 on this event and F0 off it. The
resulting response is exactly ternary, has second moment

    p_theta=(1-theta)p0+theta p1,

and has first-chaos inverse

    K_theta=(1-theta)K0+theta K1.

The gate creates no first-chaos coefficient: its indicator is even,
and it is independent of all old variables. Nor do unused Gaussian
coordinates create coefficients. Thus the exact residual variance is

    t_theta=p_theta-||K_theta||2^2.

This supplies a concrete realization of the conditional feasible-pair
purification. It is NOT the residual variance of the raw convex-average
function (1-theta)F0+theta F1, whose squared norm would be different.

The old/new first-chaos inner product used in the code is exactly

    <K0,K1>=lambda aV
       +gamma[aV c+aZ(<u,f>-A c)/nu].

Expanding K1 in g,f gives precisely the four terms in the implementation.
All objects are fixed before finite canonical approximation and the
matrix-order limit, so no growing-depth algorithm is being presumed.

## 3. Core covariance and the bivariate conditional polynomial

For each anchor inverse h_j=U*G_j, its projection on V is

    E[h_j|V]=v_j h_(d_j)(V),
    v_j=sqrt(d_j!/product_l m_l!) product_l rho_l^(m_l).

Hence Cov(L,G_j)=v_j E[f h_(d_j)(V)]. The remaining innovation covariance
is determined by Cov(L,V)=A:

    Cov(L,Z_*)=(A-sum_j rho_j Cov(L,G_j))/s.

Set b=(Cov(L,E)-A w)/nu. It is the ACTUAL vector Cov(Z,E), so b.w=0
and ||b||<=1. These properties do not assume Z is independent of the
core. Numerically its core projection is substantial, ||b|| about .55.

The old inverse has exact form

    g=sum_j [rho_j-s beta_(d_j)v_j/(a_res N)]h_j(E)
       +(s/N) sum_(d even<=200) beta_d R_(a_res)h_d(V),

where N is the norm of the anchor-deleted unnormalized resolvent.
Partial OU in Z_* has covariance with (V,Z)

    A_r=R^2+s^2 r, B_r=s(r-1)b_*.

Pure-chaos Gaussian projection then gives

    E[h_d(V_r)|V,Z]
      =sum_k sqrt(binomial(d,k)) A_r^(d-k) B_r^k
                         h_(d-k)(V)h_k(Z).

This follows equally by the Hermite generating function; the independent
OU noise preserves the unit variance of V_r. Integrating against
r^(a_res-1) gives the resolvent coefficients.

For an anchor monomial, expand product_j(rho_j x+b_j y)^(m_j).
Its x^(d-k)y^k coefficient is multiplied by
sqrt((d-k)! k!/product_j m_j!). This agrees with the implementation's
normalization. The k=0 column reproduces the old m coefficients, and
the d=0 resolvent term cancels the removed constant exactly.

For outward arithmetic, the scalar resolvent integral has the useful
nonalternating formula

    integral_0^1 r^(a-1)(R^2+s^2 r)^j[s b_*(r-1)]^k dr
      =(-s b_*)^k sum_(ell=0)^j binomial(j,ell)
          (R^2)^(j-ell)s^(2ell)
             k!/product_(v=0)^k(a+ell+v).

After the single extracted sign, all terms are nonnegative. For rational
a and R^2,s^2 the scalar sum is rational. This avoids cancellation from
expanding (r-1)^k and removes the need for numerical quadrature.

## 4. Weighted conditional Jensen is the correct handoff

After integrating the independent gate, the effective mask is

    H_theta=(1-theta)H0+theta H1.

Write K_theta=a_theta g+k_theta(V). On each V-bin, let M be the
conditional mask mass and I the conditional integral of H_theta K_theta.
Convexity gives the legitimate lower score

    E H_theta Psi(K_theta,t_theta)
       >= sum_bins Prob(bin) M Psi(I/M,t_theta),

with contribution zero when M=0. Here

    I=a_theta E[H_theta g|bin]+E[H_theta k_theta(V)|bin].

The first expectation MUST use the bivariate projection E[g|V,Z].
Although H0 is a function of V alone, H1 depends on Z. Therefore
E[H_theta g|V] generally differs from E[H_theta|V]m(V). The code keeps
this correlation correctly; making the substitution would falsely
increase the diagnostic by about .00089 at the selected interpolation.

For a piecewise-rectangle policy in (V,Z), every needed term is a finite
product of Gaussian masses and Hermite endpoint integrals. The old mask
boundary must be included among the V-bin boundaries. Setting the new
response to sign(V) outside a finite cutoff makes its center empty there;
the tail support and first moments are still included exactly.

## 5. The Z-degree truncation penalty is valid

Let r_b=||b||<1. If r_b>0, put B=(b.E)/r_b. Orthogonality b.w=0 gives
independent V,B, and Gaussian regression gives

    Z=r_b B+sqrt(1-r_b^2)N,

with N independent of the entire core. Conditional projection from B
to Z multiplies Hermite degree k by r_b^k. Consequently the part of
E[g|V,Z] with Z-degree at least K has L2 norm at most r_b^K||g||2=r_b^K.
For r_b=0, every positive Z-degree is zero and the claim is immediate.

Take K>=1. The omitted tail has zero conditional expectation given V,
so the old-mask term has EXACTLY zero error. Only theta H1 remains.
The perspective function M Psi(I/M,t) is 1-Lipschitz in I at fixed M,t.
Summing bins, applying conditional Jensen to the absolute error, then
Cauchy--Schwarz gives the total score penalty

    |a_theta| theta r_b^K.

This remains true for arbitrary rational rectangle masks. It is a norm
bound on the full inverse tail, not an independence approximation or an
unjustified use of the smaller conditional-projection norm.

The director's generalization in
`resumed_director_conditional_birth_tail_theorem_2026_09_06.md` was also
read completely and passes. For finitely many new orthonormal Gaussian
coordinates, the largest singular value of their old-space projection
replaces r_b. Singular-value decomposition and the product Hermite
channel prove the same bound by TOTAL new degree, which is unchanged
by orthogonal rotations. If degrees <=K are retained instead, the
tail exponent is K+1, not K. With a mask in [0,1], Cauchy--Schwarz gives
the optional sharper multiplier sqrt(EH_new), so the mixed-mask error
is at most |a_theta| theta sqrt(EH_new) r_b^(K+1) in that convention.

## 6. Replay and remaining numerical obligation

Independent replay at 128 Gauss-Legendre nodes per V interval with the
degree-six polynomial surrogate gave

    old value .4314718440696391,
    best value .43205334796517925,
    theta .43694747000407147,
    wrong-independence value .43294485377968217.

These are diagnostics only. The proof of a new decimal requires freezing
rational rectangle boundaries and theta, enclosing the actual moments,
using a positive lower variance, evaluating the finite Jensen score,
and subtracting the certified truncation penalty. The Gaussian mechanism
and its limit order leave no additional structural gap found here.

## 7. Completed exact rectangle-certificate audit and replay

The complete source
`computations/resumed_response_rich_core_birth_certificate_2026_09_06.py`
and its imported interval, Hermite, and central-moment routines were read.
The policy file contributes only rational rectangle endpoints. The exact
surrogate is defined by its four displayed rational degree-six power
coefficients and the tail coefficient -9/250; it need not equal the
slightly different floating optimizer fit. Theta is exactly 437/1000.

The following normalization-sensitive implementation steps pass:

- all support and first moments include the exact |V|>2 sign(V) tails;
- doubled positive rectangles correctly use j+k even under joint parity;
- the old mask boundary is included among V-bin boundaries;
- scalar integrals factor against the new mask only because its Z
  thresholds are constant inside each declared V bin;
- the u/f cross moment is computed from the actual normalized u and f;
- degrees k=0,...,19 are retained, so the penalty uses ||b||^20;
- bins discarded to avoid division by a tiny mass are legitimately
  discarded as NONNEGATIVE true contributions;
- the final outward subtraction uses an upper bound on |a_theta| and
  on the entire conditional-tail penalty.

The independent invocation wrote a distinct output at
`tmp/resumed_bound_audit_rich_birth_replay_2026_09_06.json`. It exited
successfully, and byte comparison with the response agent's saved
certificate JSON found the files IDENTICAL. In particular:

    residual variance > .09304242935883662111075117987894229248,
    core covariance norm squared < .30275113007298663878577518855657,
    Jensen bin sum > .43205252656165672224647469271233223033,
    tail penalty < .000001453832608321729546536145959765576.

After the penalty, the exact outward interval for the certified
lower-score expression is

    [.432051072729048400516928156566372464760895594160009607711726,
     .432051072729048401335353378850585427211923740055161746488190].

Thus the safely rounded ORIGINAL signing consequence is

    liminf_n M_n/n^(3/2) >= .4320510727290484.

This improves the universal lower bound without asserting convergence,
an optimal response, or validity of the larger uncorrected conditional-
independence diagnostic.
