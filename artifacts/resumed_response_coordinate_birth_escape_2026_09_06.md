# Exact coordinate-birth escape for finite canonical-tree responses

Date: 2026-09-06. Status: new structural theorem, proved below and submitted
for independent audit. It uses the full nonlinear center functional and
the original canonical Gaussian creation dictionary. It does NOT apply
merely because a response uses finitely many arbitrary rotated Gaussian
directions with nonpolynomial inverse features.

## 1. Scope and notation

Let A be a FINITE ANCESTOR-CLOSED set of canonical old Gaussian tree
coordinates, including the edge coordinate. Write X=(G_T:T in A).
Let F(X) be a bounded jointly odd, ternary response and set

    F in {0,+1,-1},  H=1-|F| in {0,1},
    mu=EH,  q=EF^2=1-mu,
    a_T=E[F G_T],
    K=U*F=sum_{T in A} a_T h_T(X),
    tau^2=q-sum_{T in A}a_T^2,
    C=E H Gamma_tau(K),  Gamma_tau(k)=E_N|k+tau N|.

Each h_T is the normalized even Hermite child monomial whose creation
coordinate is G_T. Ancestor closure ensures K is measurable with
respect to the same finite X. The full center theorem makes C a valid
original-problem lower certificate. Assume Gaussian-a.e.-continuity if
one wants direct use of its finite-response version; otherwise use its
bounded measurable extension.

The distinction from an arbitrary finite Gaussian frame is essential.
A linear combination of countably many canonical coordinates can be
a single Gaussian direction, while its inverse creation feature has
infinitely many Hermite coefficients. The polynomial argument in
Section 2 is for the finite CANONICAL bank above, not such a frame.

## 2. A new independent creation coordinate must be available

Assume

    C > sqrt(8/(27*pi)).                              (1)

This threshold is approximately .3071, so every response relevant to
the current .43-level lower bound satisfies it. Necessarily 0<mu<1.
Also tau>0: otherwise F=P_1F is a bounded Gaussian linear functional,
which must be zero, contradicting C>0.

Define the even bounded function

    R_0(X)=H(X) Gamma_tau'(K(X))
          =H(X)[2Phi(K(X)/tau)-1].

R_0 is not identically zero. Indeed, if it vanished, K=0 on H, hence

    C=mu tau sqrt(2/pi)
      <=mu sqrt(1-mu) sqrt(2/pi)
      <=sqrt(8/(27*pi)),

where mu sqrt(1-mu) has maximum 2/(3sqrt(3)). This contradicts (1).
Positive C alone would not have justified this step: a pure-noise
state can have R_0=0 and positive value.

Expand R_0 in the complete normalized even Hermite basis of X. Each
such basis monomial h has a canonical creation coordinate G_T=U h,
whose children lie in A. Suppose every nonzero coefficient occurred
at a tree T already in A. There are only finitely many such T, so
R_0 would be a polynomial of X. An essentially bounded Gaussian
polynomial must be constant. But R_0 vanishes on {H=0}, an event of
positive probability q, so the constant would be zero, contradicting
the preceding paragraph.

Consequently there exists a normalized even Hermite monomial h_T of
X such that

    T notin A,  c_T=E[R_0 h_T] != 0.                     (2)

All children of T belong to A, so adjoining T preserves finite
ancestor closure. The canonical Gaussian G_T is independent of X.
This is the legitimate new independent coordinate; it is NOT a
matrix best-response spin declared to be fresh.

## 3. Exact nonperturbative update

Fix a finite L>=1, let G denote G_T, and put

    p=P(|G|>L),  z=E[|G|1_{|G|>L}]=2phi(L),
    eta=sign(c_T).

Define

    F'=F+H eta sign(G) 1_{|G|>L},
    H'=H 1_{|G|<=L}.                                  (3)

On {H=0}, F was already +1 or -1 and is unchanged. On {H=1}, F was
zero and is filled only in the new Gaussian tail. Thus F' remains
ternary, H'=1-|F'|, and the full cube feasibility condition is exact.
Parity is also preserved: H is even and the added sign-tail response
is odd in the new coordinate. Gaussian-a.e.-continuity is preserved
because the new threshold has zero probability of equality.

Independence and oddness of the new tail response imply that every
old first-chaos coefficient is unchanged. The new coefficient is

    a'_T=E[F'G]=eta mu z.

Other unused coefficients remain zero. It follows EXACTLY that

    mu'=mu(1-p),             q'=q+mu p,
    K'=K+eta mu z h_T,
    tau'^2=tau^2+mu p-mu^2 z^2.                       (4)

By Cauchy--Schwarz,

    z^2 <= p E[G^2 1_{|G|>L}] <= p,

so tau'^2>=tau^2. In particular the nonlinear noise does not decrease
under this finite step. No first-order approximation is being used.

K' is measurable with respect to the OLD X: the inverse creation of
the new G_T is precisely h_T(X). Consequently the new gate factors
out of its full functional,

    C'=(1-p) E H Gamma_{tau'}(K+eta mu z h_T),         (5)

Equivalently C'=(1-p)E H E_N|K+eta mu z h_T+tau' N|.

Gaussian absolute expectation is nondecreasing in its standard
deviation and convex in its shift. Therefore (2), (4), and (5) give

    C' >= (1-p)[C+mu z |c_T|],
    C'-C >= p[(1-p)mu(z/p)|c_T|-C].                  (6)

For L>=1, p<1/2, and z/p=E[|G| ||G|>L]>L. Choose any finite L with

    L>=1,  L>2C/(mu |c_T|).

Then the right side of (6) is strictly positive. A convenient stronger
choice L>=max(1,4C/(mu|c_T|)) gives the explicit lower gain C'-C>p C.
The Gaussian tail p may be extremely small, but it is strictly positive.
All new parameters and the new tree remain fixed before the matrix
order tends to infinity.

Thus EVERY finite canonical ternary response satisfying (1) admits a
strictly improving, exactly feasible, one-coordinate static expansion.

## 4. A legitimate finite-step iteration, with a countable limit

The update can be repeated: each new response remains finite,
ancestor-closed, ternary, and its certificate strictly increases, so
condition (1) persists. No fresh-product-spin assumption is made at
any step. At every fixed stage the full response theorem constructs
actual feasible sign means and an original-problem lower bound.

For an infinite sequence of such deterministic parameter choices,
the supports {|F_t|=1} increase and existing nonzero signs are never
changed. Hence F_t converges pointwise to a ternary F_infinity and

    ||F_infinity-F_t||_2^2=mu_t-mu_infinity -> 0.

The adjoints K_t converge in L2 by contraction of U*, and the nonlinear
variances converge by norm continuity (they are also monotone by (4)).
The full functional is L2-continuous for jointly feasible pairs, so

    C_t increases to C_infinity=C(F_infinity,1-|F_infinity|).

It is bounded above by the original-problem liminf and therefore by
the known 1/2 upper endpoint. This constructs a convergent sequence
of static variational certificates. It is not a convergence theorem
for M_n/n^(3/2).

The new directions may have very small Hermite coefficients, very
high canonical degree, or very large tail thresholds. Formula (6)
provides no uniform positive step size. The increasing sequence can
therefore converge to a value strictly below 1/2. Its countable limit
need not have an unused independent creation coordinate; the finite
polynomial-support argument no longer applies. No stationarity or
optimality assertion about that limit is made.

## 5. Consequences and explicit nonconsequences

There is no maximizer of the full functional among finite canonical
ternary response banks at any value above sqrt(8/(27*pi)). Through
finite approximation and purification, such banks still have the same
supremum as the entire full feasible response class. Thus that
supremum, if attained at all, requires genuinely infinite canonical
structure (or a finite rotated frame whose inverse features have such
infinite structure).

The no-finite-maximizer conclusion also applies to general finite
canonical feasible pairs at those values: first maximize H to 1-|F|
and purify using one unused even Gaussian gate. This gives a finite
canonical ternary pair with no smaller value. The strict coordinate-
birth step then gives a strictly larger value than the original pair.

The theorem identifies a nonlocal, distribution-specific escape which
the Hadamard coordinatewise plateau examples do not exclude. It also
avoids the high-value zero-strip obstruction: the perturbation uses
the existing entire zero-response mask and a new independent Gaussian
tail, not slack near zero of UH.

It does not assert that the actual matrix best-response spins are
independent, that a fixed finite bank reaches a universal sharp bound,
that the variational supremum is 1/2, or that the original sequence
converges. Those are separate questions.

## 6. Exact criterion for finite creation-closed Gaussian frames

The update has a useful extension beyond finite canonical banks. Let
E=span(G_1,...,G_d) be an orthonormal finite first-Gaussian-chaos frame
inside the countable canonical space. Assume each inverse feature

    g_j=U*G_j

is a jointly even measurable function of THIS SAME frame. This is the
finite creation-closure hypothesis; arbitrary Gaussian frames need not
satisfy it. The g_j are orthonormal by isometry. For a ternary F
measurable in the frame, its entire first projection lies in E and
K=U*F belongs to M=span(g_1,...,g_d), hence is frame-measurable.

Let R_0=H Gamma_tau'(K). If

    r=||R_0-Proj_M R_0||_2 > 0,

take h=(R_0-Proj_M R_0)/r and G_new=U h. The normalized h is even
and measurable in the old frame. Orthogonality to every g_j makes
G_new orthogonal, and therefore independent, of the entire old
Gaussian frame. Its inverse feature is h. The enlarged frame is
again creation-closed.

Every identity and inequality in Section 3 now applies with h in
place of h_T and c_T=r. There is no need for h to be a polynomial:
the full countable-field response theorem and its measurable
approximation apply. Thus a nonzero residual gives the same exact
nonperturbative strict escape.

If all inverse features are polynomials of the finite frame, the
boundedness argument of Section 2 ensures r>0 at the high values in
(1). For nonpolynomial inverse features r can vanish. This is the
precise obstruction, and is why finite Gaussian dimension alone
does not prove nonstationarity.

After ANY successful coordinate birth, however, another residual is
available at high value. All inverse features of the enlarged frame
are measurable in the PREVIOUS frame, whereas its new response

    R'_0=1_{|G_new|<=L} D(old),
    D(old)=H Gamma_{tau'}'(K+eta mu z h),

has independent new-coordinate gating. Its squared distance to ANY
old-frame-measurable function is at least

    p(1-p) E[D(old)^2] > 0.

The last strict inequality follows again from the pure-noise cap
and the fact that the new value remains above (1). Thus once a
finite creation-closed frame admits one birth, this particular
construction admits all subsequent births, including when its
inverse features become nonpolynomial. This is still not a uniform
gain statement.

## 7. The canonical finite-resolvent frame: genuine cyclic closure

The exact construction in `fresh_finite_anchor_fixed_point_2026_09_05.md`,
Sections 1--2, uses E_0=(21 canonical anchors,Z_*), where Z_* is a
same-space Gaussian fixed point. Although Z_* has countably many
canonical coefficients, its inverse feature h(G,Z_*) is explicitly
a FINITE degree-200 polynomial of this 22-dimensional Gaussian frame.
Consequently

    V=rho.G+sZ_*,
    g=U*V=sum rho_T h_T+s h(G,Z_*)

is also a finite polynomial in E_0. This does not make the frame
causal or triangular: the inverse of Z_* depends on Z_* itself.
The positive innovation derivative energy in the exact certificate
confirms that this self-dependence is genuine. It is a contractive
cyclic closure, not a finite canonical tree-depth construction.

For the FULL center response one must additionally include W=UH,
where H=1_{|V|<=alpha}. This introduces the NONPOLYNOMIAL inverse
feature H(V). After orthogonalizing W against E_0, the full frame
has inverse-feature span

    M=M_0+span(H(V)),

where every element of M_0 is a polynomial of E_0. The extra Gaussian
direction is nondegenerate: H cannot belong to a finite polynomial
space, so its inverse-isometric image W cannot belong to E_0.
All inverse features of the full frame remain E_0-measurable.

The explicit first-projection calculation gives

    K=lambda g+gamma H,
    lambda=a-bc/d > 0,  gamma=b/d,

with the notation of the exact full-center certificate. Thus the
full-frame residual is nonzero. If R_0=H Gamma_tau'(K) belonged to M,
write R_0=P+kappa H with polynomial P of E_0. On the open set
{|V|>alpha}, P vanishes almost everywhere, hence P is identically
zero. On the open strip {|V|<alpha}, this would force
Gamma_tau'(lambda g+gamma)=kappa. Since tau>0 and Gamma_tau' is
strictly increasing, g would be constant on that strip and therefore
constant everywhere as a polynomial. This is impossible: ||g||_2=1
and Eg=rho_edge=.8108, so its variance is positive.

Therefore the actual canonical full-response frame, despite its
cyclic core and nonpolynomial added H feature, satisfies the precise
residual criterion of Section 6. It admits a direct static coordinate
birth and then an infinite succession of such births. This conclusion
does not rely on incorrectly treating the canonical fixed point as
a finite causal frame, and does not require a preliminary finite-
canonical approximation to obtain the structural strict escape.
