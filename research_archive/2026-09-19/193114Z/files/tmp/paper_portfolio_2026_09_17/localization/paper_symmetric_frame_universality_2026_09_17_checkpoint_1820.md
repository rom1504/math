# Symmetric scalar frame universality with a common Gibbs bound

2026-09-17. Self-contained proof jointly derived by the director and
localization track, independently reconstructed by the discrepancy track.
The Bernoulli track has an independent product-mixture proof. External
novelty is NOT established. This is a sharper all-offset transfer for the
balanced mode construction, not a comparison to iid physical-edge Gaussian
disorder and not a proof of the original minimizer convergence.

## 1. Finite all-offset theorem

Let c_1,...,c_K be deterministic coefficient vectors in R^d, K>=2,
and let h_1,...,h_K be arbitrary deterministic real offsets. Define

    B=max_(j,b)|c_(j,b)|,
    W_4=max_j sum_b c_(j,b)^4,
    V=max_j sum_b c_(j,b)^2,

so W_4<=B^2 V. Let epsilon have independent fair sign coordinates and
g have independent standard Gaussian coordinates. For every tau>0,

    |E max_j(h_j+c_j dot epsilon)-E max_j(h_j+c_j dot g)|
      <= 2log(K)/tau+(19/12)tau^3 exp(4tau B) W_4.      (1)

If B=0 the two quantities are identical. If tau B is bounded and tau
is chosen to balance the terms, this gives in particular

    error = O((B^2 V)^(1/4) [log K]^(3/4)),            (2)

provided B^2 log K<=C_0 V for a fixed C_0; the constant can depend on
C_0. The more precise W_4 in (1) can give a stronger bound.

The important distinction from separate-coordinate Taylor estimates
is that the fourth-power coefficients are summed at ONE common Gibbs
measure. The bound uses max_j SUM_b c_(j,b)^4, not
SUM_b max_j c_(j,b)^4.

### 1.1 Unconditional fourth-budget bound by querywise truncation

The director observed that the condition involving B in the optimized
rate can be removed completely. The localization track independently
checked the reduction and its constants. For L=log(2K),

    |E max_j(h_j+c_j dot epsilon)-E max_j(h_j+c_j dot g)|
                              <=8 W_4^(1/4) L^(3/4).   (2a)

This holds for every finite nonempty query family and every offset,
without a coefficient cap. When W_4=0, or when K=1, the assertion is
immediate. Otherwise put B0=(W_4/L)^(1/4), and split EACH QUERY into
the coordinates of magnitude at most B0 and those larger than B0.
For every deleted tail r_j,

    ||r_j||_1 <=W_4/B0^3,
    ||r_j||_2^2<=W_4/B0^2.

Deleting the tail therefore changes the Rademacher maximum pointwise
by at most W_4/B0^3. The expected Gaussian maximum changes by at most

    E max_j |r_j dot g|<=sqrt(2L W_4/B0^2).

These are deterministic maximum/triangle inequalities; neither deleted
part has to be independent of the retained part. Apply (1) to the
truncated queries with tau=1/(2B0). Their coefficient bound is B0 and
their fourth budget is at most W_4. Writing E0=W_4^(1/4)L^(3/4), the
two deletion costs and the two terms of (1) total at most

    [1+sqrt(2)+4+19e^2/96] E0 <8 E0.

This proves (2a). The truncation is a proof device for the affine query
family, not an edit of the physical signing or a purported new common
parent matrix. Both endpoints of (2a) are the ORIGINAL query maxima.

The split has the same functional roles as the Bernoulli-process
ell1/Gaussian decomposition: the large coefficients are paid in ell1
on the sign side and in Gaussian width on the Gaussian side. No
decomposition theorem is imported; the threshold calculation is explicit.

## 2. Exact interpolation identity

Put

    F(z)=tau^(-1)log sum_j exp[tau(h_j+c_j dot z)],
    Z_t=sqrt(t)epsilon+sqrt(1-t)g,
    Phi(t)=E F(Z_t).

Let u=sqrt(t), and condition on all signs except epsilon_b and on
the Gaussian vector. Write z^(b) for the resulting vector whose b-th
coordinate is sqrt(1-t)g_b, while all other coordinates agree with Z_t.
Gaussian integration by parts and the exact two-point sign average give

    Phi'(t)=(1/2)sum_b E_[except epsilon_b]
       [(1/(2u))integral_(-u)^u F_bb(z^(b)+s e_b) ds
          -(F_bb(z^(b)+u e_b)+F_bb(z^(b)-u e_b))/2].    (3)

This has no approximate integration-by-parts premise. For a twice
differentiable scalar function f,

    (1/(2u))integral_(-u)^u f(s)ds-[f(-u)+f(u)]/2
      =-(1/(4u))integral_(-u)^u (u^2-s^2)f''(s)ds.     (4)

Its absolute value is at most (u^2/3)sup|f''|; the outer factor 1/2
in (3) makes the derivative contribution at most
`(t/6)sup|F_bbbb|`. Formula (4) follows directly by two integrations
by parts; it is also exact for quadratics, fixing its normalization.

## 3. Fourth cumulants at one common environment

Let p_z be the Gibbs weights at z. The fourth derivative is

    F_bbbb(z)=tau^3 [E_(p_z)(c_b-mu_b)^4
                              -3 Var_(p_z)(c_b)^2].

The elementary inequalities

    E(c_b-mu_b)^4<=16 E c_b^4,
    Var(c_b)^2<=(E c_b^2)^2<=E c_b^4

give

    |F_bbbb(z)|<=19tau^3 E_(p_z)c_b^4.                 (5)

Every point z^(b)+s e_b in (3) differs from the ACTUAL Z_t only in
coordinate b, by at most 2sqrt(t)<=2. Therefore every Gibbs weight
at the shifted point is at most exp(4tau B) times its weight at Z_t.
This holds for either value of epsilon_b. Consequently

    |F_bbbb(z^(b)+s e_b)|
       <=19tau^3 exp(4tau B) E_(p_(Z_t)) c_b^4.         (6)

Average the right side over epsilon_b if necessary. It now uses the
SAME random Gibbs measure for every b. Summation is legitimate before
taking a maximum over queries:

    sum_b E_(p_(Z_t)) c_b^4
           =E_(p_(Z_t)) sum_b c_b^4 <= W_4.

Combining (3)--(6) proves

    |Phi'(t)|<=(19/6)t tau^3 exp(4tau B) W_4.

Integrating t from zero to one gives the 19/12 term in (1). The
maximum/softmax discrepancy is at most log K/tau at each endpoint.
All differentiation can first be carried out on compact subintervals
of (0,1); the displayed uniform derivative bound and finite Gaussian
moments permit passage to the endpoints. Arbitrary offsets disappear
from the estimates but remain in the Gibbs weights throughout.

## 4. Balanced-mode specialization

Use the sign bridge from
[the balanced-mode construction](paper_localization_balanced_modes_2026_09_17.md).
For a full-parent witness (sigma,x,y), the scalar coefficient of the
bundle sign g_(j,t) is

    c_(j,t)=sigma y_j sum_b h_(a(j),b) x_(b,t).

Hence |c_(j,t)|<=k. The leftover coordinates have coefficients of
absolute value one. Deterministic balance gives, for EVERY query,

    sum_b c_b^2 <=(q+k)n,
    sum_b c_b^4 <=k^2(q+k)n.                           (7)

The children contribute arbitrary offsets
sigma[H_A(x)+H_D(y)]. The absolute cap has at most 2^(n+q+1) witnesses.
Taking q=floor(epsilon n) for fixed epsilon>0, and choosing
tau=n^(-1/4)k^(-1/2), yields

    |E Q([A,C;C^T,D])-E Q([A,C^G;(C^G)^T,D])|
                              <=C_epsilon n^(5/4) k^(1/2)           (8)

whenever k<=c sqrt(n), with constants depending on a fixed c.
Thus k=o(sqrt(n)) gives o(n^(3/2)) loss, uniformly in BOTH children.
For k<=n^(1/2-delta), the loss is O_epsilon(n^(3/2-delta/2)).

The Gaussian bridge uses the SAME deterministic mode orientations and
replaces scalar signs by independent Gaussians. Its physical-coordinate
covariance is singular and anisotropic. Replacing it by iid edge
Gaussians would discard essential information and is not claimed.

The earlier separate-block third-order estimate
O(n^(4/3)k^(2/3)) and separate-scalar fourth-order estimate
O(n^(5/4)k^(3/4)) remain valid, but (8) is stronger. The improvement
comes from deterministic balance plus the common Gibbs measure, not
from unproved independence between query coefficients and the Gibbs
weights.

The unconditional bound (2a) also applies to any prescribed subset of
the full-parent witnesses with arbitrary offsets. In particular, at a
critical scaling where a selected sector has W_4<=theta n^3 and
log(2K)=O(n), its expected sign/Gaussian maxima differ by at most

    C theta^(1/4) n^(3/2).                              (9)

Thus a vanishing fourth-mass sector transfers even when an individual
coefficient is too large for a naive common-temperature argument. The
complementary coherent sector must still be bounded separately. An
upper bound for one sector is not an upper bound for the original cap
until every witness and all child offsets have been included.

## 5. Archive and verification scope

The archived `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`
uses conditional independent replacement followed by a two-star Stein
kernel contraction for bounded-spectrum Gaussian signs. Its displayed
replacement is third order and its covariance is a different object.
The archived fixed-frame theorem also uses third-order separate-block
replacement. These are adjacent mechanisms, not an already established
version of (1) identified in this audit. No external novelty claim follows.

All claims above have self-contained proofs. A finite replay is supplied
separately to check the trapezoid constants, balanced fourth-power sum,
and numerical sign/Gaussian softmax differences; it is not a premise.

Both peer researchers independently audited the complete frozen proof
and returned PASS; the director independently reconstructed the same
smart interpolation. The Bernoulli track's product-mixture proof is
retained in its Section 19 as an independent verification.

The replay is
`computations/paper_localization_2026_09_17_symmetric_frame.py`, with
output `tmp/paper_portfolio_2026_09_17/localization/symmetric_frame_audit.json`,
seed 2026091705. It checks the exact trapezoid identity on degrees
zero through eight; enumerates 4,672 balanced-frame queries for exact
second/fourth coefficient sums; and performs 1,000 cumulant diagnostics
and deterministic Gauss--Hermite softmax comparisons. The Gaussian
quadrature is not certified integration and is not used in the proof.

