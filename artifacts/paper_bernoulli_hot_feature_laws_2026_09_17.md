# Hot quadratic feature tilts: actual signs, uniform response, covariance error

2026-09-17. Independent reconstruction of the director's hot-feature
proposal. This is a physical sign-law theorem, not a statement about the
geometry of unknown minimizing ground codes. The paired cold/hot repair
and any parent construction are separate obligations.

**Sharpened range:** Section8 proves the stronger covariance error
O_F(sqrt(n) p_max), without a partition-smallness hypothesis. Under
p_max=O(r/n), combining it with the already uniform scalar estimate
extends the HOT realization to r=o(sqrt(n)). The original entrywise
bound and its r=o(n^(1/3)) corollary are retained as valid weaker proofs.
Section9 further proves the HOT scalar response and operator covariance
conclusions whenever p_max tends to zero, including diffuse r=o(n).
It does not assert high-dimensional Gaussian total-variation convergence.

## 1. Finite theorem

Let U be an n-by-r real matrix with orthonormal columns, let P=UU^T,
write u_i for its rows and p_i=||u_i||^2, and set

    d_2=sum_i p_i^2,       p_max=max_i p_i.

Fix v>1 and a=1-1/v. On ACTUAL physical signs h in {+-1}^n define

    nu(h)=2^(-n) exp[a||U^T h||^2/2]/Z.

The law is centered by global sign symmetry and has full support.
For d_2 sufficiently small, with constants depending only on v,

    Z/v^(r/2)=1+O_v(d_2),                           (1)
    Cov(nu)=I+(v-1)(P-diag P)+E,
    |E_ij|<=C_v sqrt(p_i p_j)(d_2+p_i+p_j), i!=j,
    E_ii=0,
    ||E||_F<=C_v r(d_2+p_max).                     (2)

Uniformly over EVERY Boolean query x, with t_x=||U^T x||^2/n in[0,1],

    E_nu|h.x|/sqrt(n)
       =kappa sqrt[1+(v-1)t_x]
         +O_v(n^(-1/2)+r/n+sum_i p_i^(3/2)/sqrt(n)+d_2),   (3)

where kappa=sqrt(2/pi). There is no entropy, cardinality, or adaptivity
qualification on this query-uniform conclusion: U is fixed, and the
same numerical bound holds at every x.

In particular, if p_max<=Lr/n for fixed L and r=o(n^(1/3)), then

    Z/v^(r/2)=1+O_(v,L)(r^2/n),
    ||E||_F=O_(v,L)(r^3/n)=o(1),
    sup_x response error in(3)=O_(v,L)(r^2/n+n^(-1/2))=o(1). (4)

The scalar conclusion alone needs only r=o(sqrt(n)) under this leverage
bound. The smaller rank range in(4) pays the stated Frobenius covariance
error. The constants are for a fixed hot temperature v, not uniformly
as v tends to infinity.

## 2. Hubbard representation and the relative partition estimate

Introduce a standard r-dimensional Gaussian G. The elementary Gaussian
MGF identity gives

    exp(a||U^Th||^2/2)=E_G exp(sqrt(a)G.U^Th),
    Z=E_G product_i cosh(sqrt(a)u_i.G).

Set z_i=sqrt(a)u_i.g and

    R(g)=sum_i [z_i^2/2-log cosh(z_i)].

For ALL real z,

    0<=z^2/2-log cosh z<=z^4/12.                    (5)

The upper bound follows by integrating tanh z>=z-z^3/3 for z>=0;
the latter follows from z^2-tanh^2 z>=0 after differentiation.
Since sum_i z_i^2=a||g||^2, changing the Gaussian reference measure to
gamma_v=N(0,vI_r) gives the EXACT formulas

    Z=v^(r/2)c,       c=E_(gamma_v) exp(-R),
    d rho/d gamma_v=exp(-R)/c.                     (6)

Under rho, conditionally on g, the physical h_i are independent signs
with means tanh z_i. Furthermore

    1-c<=E_(gamma_v)R
         <=a^2 v^2 d_2/4.                          (7)

Consequently c>=1/2 eventually and the auxiliary mixing law is within
O_v(d_2) in total variation of gamma_v. The later unbounded observable
comparison does NOT merely multiply this TV estimate by a nonexistent
uniform bound.

## 3. Covariance, including the entrywise error

For i!=j, Cov(nu)_ij=E_rho[tanh z_i tanh z_j]. Under gamma_v,

    E z_i z_j=av P_ij=(v-1)P_ij.

The pointwise inequalities |tanh z|<=|z| and
|tanh z-z|<=|z|^3/3 show that replacing the tanh product by z_i z_j
costs at most

    C_v sqrt(p_i p_j)(p_i+p_j)

in expectation. Here rho<=2 gamma_v, and ordinary Gaussian moments
give E|z_i|^3|z_j|<=C_v p_i^(3/2)p_j^(1/2).

For the change of auxiliary law, use(6), |1-exp(-R)|<=R and(7):

    |E_rho z_i z_j-E_(gamma_v) z_i z_j|
      <=2 E_(gamma_v)|z_i z_j|R
          +2(1-c) E_(gamma_v)|z_i z_j|
      <=C_v sqrt(p_i p_j)d_2.                      (8)

For the final inequality, sum(5) and apply Cauchy--Schwarz to each
E|z_i z_j|z_l^4. This is bounded by C_v sqrt(p_i p_j)p_l^2,
with no independence assumption between these Gaussian linear forms.
Combining proves(2). Its Frobenius bound follows from sum_i p_i=r:

    sum_(i,j) p_i p_j(d_2+p_i+p_j)^2
        <=r^2(d_2+2p_max)^2.

The diagonal is EXACTLY one, explaining why diag P must be subtracted
from the proposed leading covariance. No diagonal approximation is used.

## 4. A conditional scalar comparison valid even at tiny variance

We use the following elementary Wasserstein form of scalar Stein's
normal approximation. If independent centered variables Y_i satisfy
|Y_i|<=b and sum_i E Y_i^2=sigma^2>0, then for every1-Lipschitz test H,

    |E H(sum_i Y_i)-E H(sigma G)|
       <=3 sum_i E|Y_i|^3/sigma^2<=3b.             (9)

For completeness, solve sigma^2 f'-yf=H-EH(sigma G); the standard
one-dimensional Stein solution has ||f''||_infinity<=2/sigma^2.
Writing S_i=S-Y_i, insert sigma_i^2 E f'(S_i) in each summand and
Taylor-expand f(S_i+Y_i). The error is at most

    (2/sigma^2)sum_i[variance(Y_i)E|Y_i|+E|Y_i|^3/2].

The product of the first and second absolute moments is at most the
third absolute moment, proving(9). The usual Stein derivative bound
is the same one reconstructed in the campaign's
[shifted-absolute proof](paper_bernoulli_shifted_absolute_2026_09_17.md).
Approximation handles nonsmooth Lipschitz H. In the zero-variance case
the comparison is exact. Thus(9) has no lower-variance hypothesis.

Conditionally on g, let

    m_x(g)=sum_i x_i tanh(z_i)/sqrt(n),
    sigma(g)^2=sum_i sech^2(z_i)/n.

The centered normalized sign summands have absolute value at most
2/sqrt(n). Applying(9) with H(y)=|m_x+y| yields, UNIFORMLY in g,x,

    |E(|h.x|/sqrt(n) | g)-E_G|m_x(g)+sigma(g)G||<=6/sqrt(n). (10)

Moreover

    0<=1-sigma(g)^2<=sum_i z_i^2/n=a||g||^2/n,

so replacing sigma G by G costs O_v(r/n) after rho averaging. The
linearized conditional mean is

    l_x(g)=sqrt(a)(U^Tx/sqrt(n)).g,

and its mean absolute error is bounded by

    E_rho|m_x-l_x|
       <=sum_i E_rho|z_i|^3/(3sqrt(n))
       <=C_v sum_i p_i^(3/2)/sqrt(n).              (11)

This controls EVERY query by the same deterministic error.

## 5. The unbounded auxiliary comparison and the final scalar limit

Put F_x(g)=E_G|G+l_x(g)|. Since ||U^Tx/sqrt(n)||<=1,

    F_x(g)<=kappa+|l_x(g)|,
    E_(gamma_v) l_x(g)^2<=v-1.

Using(6)--(7) directly,

    |E_rho F_x-E_(gamma_v)F_x|
       <=2 E_(gamma_v)F_x R+2(1-c)E_(gamma_v)F_x
       <=C_v d_2.                                 (12)

Indeed E|l_x| z_i^4<=C_v p_i^2 by Gaussian moments and the uniform
variance bound for l_x. This is the required uniform integrability
calculation; it avoids a spurious sqrt(r) loss from bounding by ||g||.

Finally, under gamma_v and an independent standard scalar G,
G+l_x(g) is Gaussian with variance1+(v-1)t_x. Equations(10)--(12)
give(3), and the leverage bounds give(4).

## 6. Scope and checks

This theorem supplies genuine cube measures and uniform scalar
response approximations with quantified covariance errors. It does
NOT claim exact isotropy of a temperature mixture without an additional
finite covariance repair, a dimension-free subGaussian proxy for an
arbitrary such repair, or a favorable complete-parent cap. Nor does
it establish that an actual minimizing nearcode has large projection
onto the selected rank-r feature space.

Finite normalization checks are in
`computations/paper_bernoulli_2026_09_17_hot_feature_laws.py`.
Their numerical quadrature is explicitly a replay, not the proof of
the uniform asymptotic theorem.

The discrepancy track independently read and reconstructed the entire
proof: PASS, including the tiny-variance Stein step and the unbounded
auxiliary comparison. The replay checks six frames, including diffuse
Hadamard and nonconstant-leverage frames, and all2,112 Boolean queries
at their dimensions. Exact-cube and Gaussian-quadrature partitions and
covariances agree within2e-7 (observed errors about1e-15). Rank-one
binomial replays at orders16 through1,024 show the stated limiting
normalization and response. Compilation also passes. The output is
`tmp/paper_portfolio_2026_09_17/bernoulli/hot_feature_laws.json`.

## 7. Exact hot subGaussianity and entropy, without a rank restriction

The discrepancy track supplied the following concentration strengthening;
the director and this track independently reconstructed its algebra.
For EVERY n,r,U and fixed v>1, the hot law satisfies

    E_nu exp(theta.h)<=exp(v||theta||^2/2).         (13)

No small leverage or approximation is required for(13). Under an
arbitrary external field theta, the Hubbard mixing potential is

    Phi_theta(g)=||g||^2/2-sum_i log cosh(theta_i+sqrt(a)u_i.g).

Its Hessian lies between I/v and I. The classical strong-log-concavity
Poincare inequality therefore has constant v. For clarity, this
particular consequence follows directly from the Langevin semigroup:
the derivative flow contracts by e^(-t/v), so
||grad P_t F||^2<=e^(-2t/v)P_t||grad F||^2. Invariance and the variance
identity 2 integral_0^infinity E||grad P_tF||^2 dt give
Var(F)<=v E||grad F||^2. The potential here is smooth with globally
bounded Hessian and a uniform positive lower Hessian, so the derivative
flow, invariant density and limiting variance identity apply by ordinary
smooth approximation. This is the strong-convexity case of the
[classical Brascamp--Lieb variance principle](https://doi.org/10.1016/0022-1236(76)90004-5),
not a new concentration theorem.

Let D=diag(sech^2(theta_i+sqrt(a)u_i.g)). Conditional variance plus
this Poincare inequality gives, for any fixed vector b,

    Var_(nu,theta)(b.h)
      =E b^T D b+Var_rho(sum_i b_i tanh(theta_i+sqrt(a)u_i.g))
      <=E b^T D b+av E||U^T D b||^2
      <=(1+av)||b||^2=v||b||^2.

Thus the log moment-generating function has Hessian at most vI at
EVERY external field. Its value and gradient vanish at zero, proving
(13) by integration along a line.

There is also a sharper hot entropy calculation than one obtains by
tracing the Frobenius covariance error. Under the unshifted Hubbard law,
Gaussian integration by parts in the reference gamma_v gives

    E_rho||g||^2
      =vr-v E_rho sum_i z_i(z_i-tanh z_i).          (14)

The summands in the correction are nonnegative and at most z_i^4/3.
When d_2 is small, rho<=2 gamma_v therefore bounds the correction by
O_v(d_2). Differentiating the finite cube partition, or equivalently
integrating by parts in the original standard-Gaussian representation,
gives the EXACT identity

    E_nu||U^Th||^2=(E_rho||g||^2-r)/a.

Consequently

    E_nu||U^Th||^2<=vr,
    E_nu||U^Th||^2=vr+O_v(d_2),
    D(nu||Uniform)=(r/2)(v-1-log v)+O_v(d_2).       (15)

The first inequality in(15) holds without the small-leverage condition;
the error bound uses it. This does not differentiate an unspecified
partition-error term, and its error is O(r^2/n) under diffuse leverage.

Finally there is a safe exact concentration reduction for a cold
temperature1/2<v<1. Write beta=1/v-1<1. On the PHYSICAL cube,

    exp[-beta h^TPh/2]
      =exp[-beta n/2] exp[beta h^T(I-P)h/2].

Apply(13) to an orthonormal basis of the complementary projection.
The cold law is therefore v/(2v-1)-subGaussian, regardless of rank or
leverage. The discrepancy track's stronger Fourier proof treats ALL
fixed cold v>0 with asymptotic proxy1+o(1) in the diffuse growing-rank
regime; the complement argument is only an exact simpler fallback.

## 8. Sharpened covariance via latent Fisher information

The director proposed the following stronger bound, reconstructed here
in full. For EVERY n,r,U, with fixed v>1,

    Cov(nu)=I+(v-1)(P-diag P)+E,
    ||E||op<=C_v p_max,
    ||E||F<=C_v sqrt(n) p_max.                     (16)

There is no d_2-smallness requirement in(16). In particular, diffuse
leverage p_max<=Lr/n gives Frobenius error O_(v,L)(r/sqrt(n)); combining
with(1),(3) and(15) gives all three hot realization estimates whenever
r=o(sqrt(n)). The covariance improvement does not by itself improve a
separate cold argument or an exact-isotropy repair.

### 8.1 The latent law is nearly Gaussian at the covariance level

In the representation(6), let Phi(g)=||g||^2/(2v)+R(g). Direct
differentiation gives

    Hess R=a U^T diag(tanh^2 z_i) U >=0.           (17)

Thus rho is centered and1/v-strongly logconcave, irrespective of rank
or leverage. The Poincare argument in Section7 gives Cov_rho g<=vI.
Applied after arbitrary linear tilts of g, it also gives

    E_rho exp(t.g)<=exp(v||t||^2/2).

Consequently E z_i^2<=av p_i and E z_i^4<=C_v p_i^2. The latter needs
only the elementary subGaussian tail integral; an equality with a
Gaussian fourth moment is not assumed.

Let J=E_rho Hess Phi. From(17) and tanh^2 z<=z^2,

    I/v<=J<=I/v+a^2 v p_max I.                    (18)

The needed matrix Cramer--Rao bound is elementary here. Integration by
parts gives

    E[g (grad Phi)^T]=I,
    E[(grad Phi)(grad Phi)^T]=E Hess Phi=J.

Both vectors are centered. Positivity of their joint covariance, or
minimizing E||c.g-d.grad Phi||^2 over d, yields

    Cov_rho g>=J^(-1).

Combining with the upper Poincare bound,

    v/(1+a^2 v^2 p_max) I<=Cov_rho g<=vI,
    ||Cov_rho g-vI||op<=a^2 v^3 p_max.            (19)

The potential is smooth, has bounded Hessian and a positive quadratic
lower bound; all these integration-by-parts identities have Gaussian
tail justification without an auxiliary truncation limit.

### 8.2 The nonlinear conditional mean has a small operator error

Write z=sqrt(a)Ug, m=tanh z, and e=m-z, all centered under rho.
For any physical coefficient vector b,

    grad_g(b.e)=-sqrt(a)U^T diag(tanh^2 z_i)b.

Poincare and ||U^T||op=1 imply

    Var_rho(b.e)
      <=va E||U^T diag(tanh^2 z_i)b||^2
      <=va sum_i b_i^2 E tanh^4 z_i
      <=C_v sum_i b_i^2 p_i^2
      <=C_v p_max^2 ||b||^2.                      (20)

Therefore ||Cov e||op<=C_v p_max^2. Meanwhile
Cov z=aU(Cov g)U^T has bounded operator norm. Cauchy--Schwarz for the
cross covariance gives

    ||Cov(z,e)||op
       <=sqrt(||Cov z||op ||Cov e||op)<=C_v p_max.

Using(19),

    Cov m=(v-1)P+O_op,v(p_max).                   (21)

The physical conditional variance is diagonal:

    Cov_nu h=E_rho diag(sech^2 z_i)+Cov_rho m.

Since E tanh^2 z_i<=av p_i, the diagonal term is I+O_op,v(p_max).
Together with(21) this gives I+(v-1)P+O_op,v(p_max). Subtracting the
explicit diagonal (v-1)diag P costs at most(v-1)p_max in operator
norm and makes the leading term's diagonal exactly one. This proves
the operator bound in(16); its Frobenius bound follows from
||E||F<=sqrt(n)||E||op.

This argument exploits strong logconcavity, a matrix score identity,
and conditional-mean Poincare control jointly. It does not replace a
random covariance by its expectation inside a moment-generating
function, nor infer a distributional Gaussian approximation merely
from covariance convergence. The scalar theorem remains proved by
the separate uniform comparison in Sections4--5.

## 9. Directional Stein comparison: hot response for diffuse r=o(n)

The director's further proposed extension also follows from the same
latent strong logconcavity, without comparing the entire r-dimensional
auxiliary law to a Gaussian. For fixed v>1 and EVERY n,r,U,

    sup_(x Boolean) |E_nu|h.x|/sqrt(n)
        -kappa sqrt[1+(v-1)||U^Tx||^2/n]|
       <=C_v(p_max+n^(-1/2)).                     (22)

Together with(16), this proves operator-covariance convergence and the
uniform scalar response whenever p_max->0; under diffuse leverage this
permits r=o(n). The physical hot law remains exactly v-subGaussian by
(13). Small Frobenius error and exact isotropic mixtures are separate
requirements and must not be inferred from the larger range in(22).

The same error is UNIFORM in every real affine offset s:

    sup_(x Boolean,s real) |E_nu|s+h.x/sqrt(n)|
       -E_G|s+sqrt(1+(v-1)||U^Tx||^2/n)G||
       <=C_v(p_max+n^(-1/2)).                     (22a)

The conditional scalar test is still1-Lipschitz, and the final
directional test y->E_G|s+G+y| has the same Lipschitz constant for all
s. Thus the proof below changes neither its constants nor its order
of quantifiers. It does not multiply a constant depending on s by a
total-variation error.

### 9.1 Every auxiliary linear projection has a small Wasserstein error

Fix b in R^r, put L=b.g and sigma^2=v||b||^2. Integration by parts
under rho, whose potential is ||g||^2/(2v)+R, gives

    E[L f(L)]
       =v||b||^2 E f'(L)-v E[f(L)b.grad R].        (23)

For a1-Lipschitz scalar test H, the solution of
sigma^2 f'-Lf=H-EH(sigma G) has ||f||infinity<=1, independently of
sigma. One direct formula is f(y)=-integral_0^infinity e^(-t)
E H'(e^(-t)y+sigma sqrt(1-e^(-2t))G) dt, with smooth approximation
for nonsmooth H. For b=0 the comparison is exact.

The function F_b=b.grad R is centered because it is odd under the
symmetric law rho, and

    grad F_b=a U^T diag(tanh^2 z_i) U b.

Poincare and the fourth-moment bounds from Section8 show

    E F_b^2
       <=va^2 sum_i (Ub)_i^2 E tanh^4 z_i
       <=C_v p_max^2 ||b||^2.

Equation(23) and the bounded Stein solution therefore prove

    W_1(L,N(0,v||b||^2))<=C_v p_max||b||.          (24)

There is NO division by a small ||b||, and no unbounded-test
total-variation comparison is involved.

### 9.2 Uniform physical response

Put c=x/sqrt(n), so ||c||=1. The conditional mean is c.tanh z;
the linearized mean is c.z=sqrt(a)(U^Tc).g. Equation(20) gives

    E_rho|c.(tanh z-z)|<=C_v p_max.                (25)

The conditional variance is the same as in Section4 and obeys

    E_rho(1-sigma(g)^2)
       <=a v r/n<=a v p_max.

Use the variance-uniform conditional scalar estimate(10), then replace
sigma(g)G by G and the mean by c.z. Apply(24) with
b=sqrt(a)U^Tc to the1-Lipschitz test y->E_G|G+y|. The resulting
Gaussian variance is1+(v-1)||U^Tx||^2/n. This proves(22), uniformly
over all queries, with no code cardinality condition.

### 9.3 Entropy at the larger rank range

The additive log-partition estimate does not need d_2 small. From(6),
Jensen and(5),

    -C_v d_2<=log E_(gamma_v)e^(-R)<=0.

Hence logZ=(r/2)logv+O_v(d_2) for ALLr. Likewise the exact identity
(14), using the strong-logconcavity fourth moments instead of
rho<=2gamma_v, gives E_nu||U^Th||^2=vr+O_v(d_2) for ALLr. Therefore

    D(nu||Uniform)=(r/2)(v-1-logv)+O_v(d_2).        (26)

Since d_2<=r p_max, its error is o(r) whenever p_max->0. This is an
additive free-energy estimate, not relative partition convergence
when r^2/n is large. That distinction matters in the diffuse r=o(n)
range and is explicitly retained here.

## 10. Trace-sensitive Frobenius refinement

The discrepancy track observed a sharper Frobenius consequence of the
same operator proof; this track independently reconstructed it:

    ||Cov(nu)-I-(v-1)(P-diag P)||F
       <=C_v sqrt(r) p_max.                       (27)

This is valid at every rank and order, with fixed v. Under diffuse
leverage it is O_(v,L)(r^(3/2)/n). It improves the HOT covariance-repair
scale but does not on its own enlarge the separate cold Fourier range.

Retain z=sqrt(a)Ug and e=tanh z-z. The coordinate version of(20) gives

    tr Cov e<=C_v sum_i p_i^3<=C_v r p_max^2,
    ||Cov e||op<=C_v p_max^2.

Thus ||Cov e||F<=sqrt(||Cov e||op tr Cov e)
<=C_v sqrt(r)p_max^2. Positivity of the covariance block matrix of
(z,e) gives the standard contraction factorization of the cross
covariance; hence

    ||Cov(z,e)||F^2<=||Cov z||op tr Cov e
       <=C_v r p_max^2.

Finally Cov z-(v-1)P has rank at most r and operator norm at most
C_v p_max by(19). Adding the three terms shows

    ||Cov(tanh z)-(v-1)P||F<=C_v sqrt(r)p_max,

using p_max<=1. Physical conditional variances restore the diagonal
EXACTLY to one. Therefore the actual covariance error is precisely
the hollow projection of this last matrix, and hollow projection is
an orthogonal contraction for the Frobenius norm. This proves(27),
without a spurious sqrt(n) diagonal cost.

The discrepancy researcher independently audited all of Sections7--9
in full, including the Fisher lower bound, variance-uniform directional
Stein solution and larger-rank entropy identity: PASS. Section10 is
the jointly reconstructed trace refinement of those same estimates.

## 11. All-offset scalar comparison retains an arbitrary new-spin child

Here is a precise conditional composition consequence, consistent with
the director's earlier codewise field-universality theorem. It states a
comparison, not a bound on the target Gaussian-mixture value.

Suppose an actual centered sign law nu has linear subGaussian proxy K,
and a centered Gaussian variance-mixture column G has the same common
linear proxy K. Assume uniformly over Boolean x and real s that

    |E_nu|s+h.x/sqrt(n)|-E|s+G.x/sqrt(n)||<=delta.   (28)

For the feature realization, G conditional on the variance label v_j
has covariance I+(v_j-1)P, which is positive definite. The mixture is
NOT silently replaced by a single covariance-I Gaussian. Independent
physical columns and independent mixture-Gaussian columns are used.

Let I be any finite declared old-query index set. Each index i has a
Boolean word x_i, any deterministic offset b_i, and ANY deterministic
function d_i on the q new signs. Define

    F_i(h_1,...,h_q)
       =b_i+max_(y in{+-1}^q)[d_i(y)+sum_j(h_j.x_i)y_j],

and define F_i(G_1,...,G_q) identically. Absolute-parent queries are
included by indexing old words and both global polarities and taking
d_i to be the corresponding signed child quadratic. No new-spin
configuration is discarded.

For a fixed index i and all other fields frozen, the maximum as a
function of the jth normalized scalar field u is EXACTLY

    constant+sqrt(n)|u-s|

for some real s: separate the maxima with y_j=+1 and y_j=-1.
Replacing the q independent fields successively and applying(28)
therefore gives

    |E F_i(h)-E F_i(G)|<=q sqrt(n) delta.          (29)

There is also a fully paid maximum-over-old-code estimate. If V_j are
independent centered K-subGaussian scalar fields, a separately
sqrt(n)-Lipschitz function of them has centered MGF proxy2Kqn.
To check this without assuming concentration of an arbitrary vector
law, use the Doob martingale and an independent scalar copy: for each
conditional1-coordinate function g,

    E exp(lambda(g(V)-E g(V)))
      <=E cosh(lambda[g(V)-g(V')])
      <=E cosh(lambda sqrt(n)(V-V'))<=exp(Kn lambda^2).

Tensorizing these conditional bounds proves the claimed proxy. Each
F_i has exactly the required separate Lipschitz constants. Thus its
expected maximum exceeds its maximum mean by a number in
[0,sqrt(4Kqn log|I|)], for BOTH ensembles. Combining with(29),

    |E max_i F_i(h)-E max_i F_i(G)|
       <=q sqrt(n) delta+sqrt(4Kqn log|I|).        (30)

The factor in front of the concentration term is one, not two: both
expected maxima lie above their maximum means and within the same
one-sided error interval. For |I|=1 that term is correctly zero.

If q=O(n), delta=o(1), and log|I|=o(n), the comparison error is
o(n^(3/2)), retaining arbitrary child couplings and offsets. Choosing
or preparing such an old-query code, controlling excluded old words,
and bounding the target mixture-Gaussian parent remain separate
obligations. In particular(30) does not claim that matching covariance
or reducing per-column means already solves the original signing
minimization problem.
