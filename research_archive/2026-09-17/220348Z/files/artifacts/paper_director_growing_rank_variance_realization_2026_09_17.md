# Growing-rank variance allocation in actual sign vectors

2026-09-17. **Proved; independent cold, hot, concentration and combined-proof reconstructions PASS.** This extends the
fixed-rank realization theorem, not the original signing bound. Constants
depend on the fixed variances and leverage constant, not on the rank.
No external novelty claim is made. The Gaussian transform and high-temperature
expansion are classical; the exact physical covariance repair is essential
to the particular conclusion below.

## 1. The theorem and limit order

Let U=U_n have r=r_n orthonormal columns, with row vectors u_i, and write
P=UU^T, p_i=P_ii. Assume

    max_i p_i <= Lr/n,              r=o(n^(1/3)).             (1)

Fix finitely many variances v_j>0 and probabilities pi_j with
sum_j pi_j v_j=1. They are fixed BEFORE n grows. There exist centered,
exactly isotropic laws nu_n on actual physical signs {+-1}^n, of full
support, such that uniformly over ALL Boolean queries x,

    E_nu |h.x|/sqrt(n)
      = kappa sum_j pi_j sqrt(1+(v_j-1)t_x)
          + O_L,v((r^2/n)^(1/4)+r^3/n),
    t_x=||U^T x/sqrt(n)||^2,       kappa=sqrt(2/pi).         (2)

Moreover

    D(nu_n || uniform signs)
       <= -(r/2) sum_j pi_j log(v_j) + o(r).                (3)

The exact-isotropy repair has mass O_L,v(r^3/n). The laws also have a
common dimension-free subGaussian proxy: any fixed constant larger than
max{3/2,max_j v_j} works at all sufficiently large n. This conclusion
requires the separate Fourier-modulus/determinant argument below; it does
not follow from covariance control. All formulas concern actual signs,
not Gaussian vectors rounded without accounting for changed covariance.

## 2. Physical component laws

For a fixed v define a=1-1/v and

    w_v(h)=v^(-r/2) exp(a ||U^T h||^2/2),
    Z_v=E_uniform w_v,      d nu_v / d uniform=w_v/Z_v.

All these finite-cube densities are strictly positive and even. The two
Gaussian-transform reconstructions prove, with S2=sum_i p_i^2,

    Z_v=1+O_v,L(r^2/n)+negligible,
    Cov(nu_v)=I+(v-1)(P-diag(P))+E_v,
    ||E_v||F <= C_v,L r^3/n+negligible,                    (4)

and uniformly in x,

    E_nu_v |h.x|/sqrt(n)
       =kappa sqrt(1+(v-1)t_x)+O_L,v((r^2/n)^(1/4)).       (5)

Complete independent branch proofs are in
`paper_localization_growing_rank_cold_tilt_2026_09_17.md` and
`paper_bernoulli_hot_feature_laws_2026_09_17.md`. The cold proof obtains its
rate from the absolute-value Fourier integral, split at delta^(1/4) and
delta^(-1/4), where delta=O(r^2/n). The hot proof in fact gives
O(r^2/n+n^(-1/2)).

The negligible error in the cold case is at most a fixed polynomial in
n,r times exp(C_v r-c_v n/(Lr)); it tends to zero faster than powers under
(1). The hot case needs no exclusion of bad Gaussian coordinates.
The case v=1 is independent signs.

### Cold branch, 0<v<1

Set beta=1/v-1 and z_i=sqrt(beta) u_i.g for standard Gaussian g in R^r.
The exact Fourier identity is

    E_h exp(-beta ||U^T h||^2/2)=E_g product_i cos(z_i).

On max|z_i|<=c<pi/2, compare the product with exp(-sum_i z_i^2/2).
The relative error is bounded by C sum_i z_i^4. After Gaussian tilting,
g has covariance vI, and its expectation is O_v(S2). Outside that set,
the union bound is <=2n exp(-c_v/pmax). Dividing by the Gaussian
normalization v^(r/2) still gives the stated negligible error.

For i!=j, the covariance numerator is exactly

    -E_g sin(z_i)sin(z_j) product_(k!=i,j) cos(z_k).

Use tan(z_i)tan(z_j) product cos ONLY on the good set. Taylor expansion
there, Gaussian moment bounds, and the original sine product outside give
entry error at most

    C_v sqrt(p_i p_j)[S2+p_i+p_j+negligible].

The reference entry is (v-1)P_ij. Frobenius summation gives
O(r S2+r pmax)=O_L(r^3/n); the physical diagonal stays exactly one.

For the query V=x.h/sqrt(n), the characteristic numerator is
E_g product_i cos(z_i+t x_i/sqrt(n)). Completing the Gaussian square
gives the reference characteristic function with variance
1+(v-1)t_x. The same fourth-moment comparison works uniformly in x for
each fixed real t. Its fourth-moment sum is
O_v,t(S2+pmax+1/n). Covariance control in (4) supplies uniformly bounded
second moments, so uniform convergence of absolute first moments follows.
The order is: fixed variances, fixed characteristic argument, then n;
tail truncation is removed last. No joint inverse covariance is used.

### Hot branch, v>1

The Hubbard--Stratonovich identity gives

    E_h exp(a||U^T h||^2/2)=E_g product_i cosh(sqrt(a)u_i.g).

Since a<1, completing the square produces the Gaussian reference with
covariance vI. Globally,

    0 <= z^2/2-log cosh z <= C z^4.

This proves the partition estimate by a fourth-moment bound. Covariance
uses tanh(z_i)tanh(z_j) times the cosh product and
|tanh z-z|<=C|z|^3, giving the same entrywise error as the cold case.

Conditioned on the auxiliary Gaussian, physical signs are independent
with means tanh(sqrt(a)u_i.g). The auxiliary law is asymptotically close
to N(0,vI); its density relative to that law is bounded by 1+o(1).
The conditional query mean differs in mean from sqrt(a)(U^Tx/sqrt(n)).g
by O_v,L(r^(3/2)/n). Its variance tends uniformly to one, and the bounded
summand scalar CLT gives (5). Uniform second moments justify removing
truncations rather than multiplying a total-variation error by sqrt(n).

## 3. Exact isotropy, with vanishing repair mass

Let tilde_nu=sum_j pi_j nu_vj. The leading hollow covariance terms in
(4) CANCEL EXACTLY, since sum pi_j(v_j-1)=0. Thus

    Cov(tilde_nu)=I+E,       f=||E||F=O(r^3/n)=o(1).

If f>0, set theta=4f/(1+4f), D=-E/(4f), and
R=I+sin((pi/2)D), with sine entrywise. Then R is a correlation matrix
with spectrum in [1/2,3/2], and the actual sign of N(0,R) has covariance
I+D by the arcsine formula. The mixture with tilde_nu at weight theta
has covariance EXACTLY I. If f=0, do not repair.

Both component covariances are bounded in operator norm, so the repair
changes every normalized absolute query by O(theta), not O(theta sqrt(n)).
This proves (2). Full physical support already holds for the component
laws; it also holds after mixing. The repair law has relative entropy at
most pi^2/64 by Gaussian data processing and the Frobenius bound.

## 4. Information cost, not a free high-dimensional miracle

The exact branch divergence equals

    D(nu_v||uniform)
       =((1-1/v)/2) Tr[P Cov(nu_v)]-(r/2)log v-log Z_v.

Formula (4) implies

    Tr[P Cov(nu_v)]
       =vr-(v-1)S2+O(sqrt(r)||E_v||F)=vr+o(r).

Consequently the branch cost is
(r/2)(v-1-log v)+o(r). Entropy convexity, sum pi_j v_j=1,
and the bounded-cost vanishing repair prove (3).

Two variances suffice to obtain arbitrarily cheap responses ON the feature
space. For v0<1<v1 use pi1=(1-v0)/(v1-v0), pi0=1-pi1. At t_x=1,
the response is kappa(pi0 sqrt(v0)+pi1 sqrt(v1))+o(1), which can be
made arbitrarily small by taking v0 small and v1 large, both still fixed.
The information cost is O(r); constants and required orders deteriorate.

## 5. Exact relation to the original problem

This is a growing-dimensional physical realization theorem and a sharp
distinction from bounded-band Gaussian response mixtures. It combines a
deterministic mean-field Gaussian transform with exact sign-covariance
repair. It does NOT show that high-energy states of actual minimizing
signings are captured by a diffuse o(n^(1/3))-dimensional subspace.

In fact a query distribution with covariance norm at most L0 satisfies
E t_x<=L0 r/n, and (2) then gives average response kappa-o(1) when r=o(n).
The low-query-covariance branch in the earlier actual-low-cap theorem
requires macroscopic rank. Thus this result cannot be substituted for
that missing geometry, nor does it close any original-value recurrence.

Classical neighboring literature includes the Hubbard--Stratonovich
method and low-storage Hopfield fluctuation theorems (Bovier--Gayrard,
WIAS Preprint283,1996). Those random-pattern results are not imported as
proofs of this deterministic all-query, exactly-isotropic statement.
External novelty remains unestablished.

## 6. The concentration obligation is also closed

For v>1, the hot law is EXACTLY v-subGaussian, for every rank and order.
Under any external field its Hubbard auxiliary potential has Hessian at
least I/v. Conditional product-spin variance plus the Gaussian
Brascamp--Lieb variance bound gives the physical covariance at most vI
under every field; integrating the log-partition Hessian proves the MGF
claim. The hot-component proof includes this argument.

For 0<v<1 the separate proof in
`paper_discrepancy_cold_feature_subgaussian_2026_09_17.md` gives

    K_cold=1+O_v,L(sqrt(r/n)(r+log n))=1+o(1).

Its decisive step is to take the MODULUS of the biased-product Fourier
integrand before a Gaussian determinant comparison. The latent ball
radius squared is O_v(r+log n), paying the exponentially small partition
normalization. Uniformly for every external test vector t this gives

    log E exp(t.h) <=(1+beta pmax)||t||^2/2+b_n,
    b_n=O_v,L(r(r+log n)^2/n)=o(1).

An even-cosh scaling inequality removes the additive defect at arbitrarily
small t: evaluate first at norm b_n^(1/4), then interpolate the positive
even power series toward zero. The resulting proxy is
1+beta pmax+O(sqrt(b_n)). No all-field covariance theorem or Gaussian
surrogate is substituted for the physical law.

The corrected mixture has proxy max{3/2,max_j v_j,1+o(1)}. Each step was
independently reconstructed by all three researchers and the director.
The same proof is uniform for the retained-rank budget used below.

## Appendix A. Adaptive covariance capture, with a uniform information budget

Independent derivation and audit by the discrepancy track. This corollary
removes the need to prescribe one diffuse subspace in advance. It does
not assert that original minimizing codes satisfy its geometric premise.

Let C_n be any nonempty Boolean code, and let 1<=r_n=o(n^(1/3)). Define

    s_r(C_n)=inf_(mu supported on C_n)
                 [sum_(j=1)^r lambda_j(E_mu xx^T)]/n,       (A1)

with eigenvalues in decreasing order. For EVERY fixed theta>0 there is
a finite constant B_theta, independent of n, the code and r_n, such that
for all sufficiently large n there is a centered, exactly isotropic,
full-support physical sign law nu_n satisfying

    max_(x in C_n) E_nu |h.x|/sqrt(n)
       <=kappa sqrt(1-s_r(C_n))+theta,
    D(nu_n || uniform signs)<=B_theta r_n.                  (A2)

The order threshold may depend on theta and on the rate r_n^3/n->0.
The law can ALSO be required to have a common finite subGaussian proxy
K_theta, independent of n, the code and the rank budget. Section6 supplies
this additional property; it is not inferred from the information bound.
The infimum in (A1) is over ALL laws on the full code, not only uniform
or extremal-energy laws.

### A.1. Uniformity over all retained ranks and all dual laws

Fix a small eta>0 and a dual law mu. Put Sigma=E_mu xx^T. Among the
top r covariance eigenvectors, retain those whose eigenvalues are at
least eta n/r. Let U contain these q<=r orthonormal vectors, and P=UU^T.
The retained mass and leverage obey

    Tr(P Sigma)/n >= s_r(C_n)-eta,
    P_ii <= Sigma_ii/(eta n/r)=r/(eta n).                   (A3)

The second inequality uses Sigma_ii=1 and the spectral inequality
Sigma >= (eta n/r)P. It does not assume that the original top-r
eigenvectors were diffuse. If no eigenvector is retained, use the
independent sign law; the estimates below also have their evident
rank-zero interpretation.

The component proofs remain uniform with this rank BUDGET r, even when
the actual retained rank q is much smaller. Indeed

    S2=sum_i P_ii^2 <= q max_i P_ii <= r^2/(eta n),
    ||E_v||F <= C_v q(S2+max_i P_ii)+negligible
               <=C_(v,eta) r^3/n+negligible.               (A4)

The cold tail is bounded by a polynomial times
exp(C_v q-c_v eta n/r), uniformly by the same expression with q replaced
by r. Its negligibility is therefore uniform over mu and q. The hot
error bounds likewise use only q<=r, (A3), and (A4). The uniform scalar
error tends to zero at rate at most
O_(v,eta)((r^2/n)^(1/4)+r^3/n).

Choose fixed cold/hot variances and weights with sum pi_j v_j=1 and

    m_theta=kappa sum_j pi_j sqrt(v_j)<=theta/4.

For t_x=||U^T x||^2/n, the response profile satisfies

    kappa sum_j pi_j sqrt(1-t_x+v_j t_x)
       <=kappa sqrt(1-t_x)+m_theta sqrt(t_x)
       <=kappa sqrt(1-t_x)+m_theta.

Jensen's inequality and (A3) bound its mu average by

    kappa sqrt(1-(s_r(C_n)-eta)_+)+m_theta
       <=kappa sqrt(1-s_r(C_n))+kappa sqrt(eta)+m_theta.     (A5)

Choose eta so that kappa sqrt(eta)<=theta/4 and then take n large enough
for the uniform realization error to be at most theta/4. This constructs
an admissible physical law for EVERY dual mu with average payoff at
most the right side of (A2), with a remaining theta/4 margin.

### A.2. The entropy budget survives adaptive choice and minimax

For the q-dimensional retained feature space, the branch entropy
calculation in Section 4 is uniform. Specifically its trace remainder
is at most sqrt(q)||E_v||F=O_(v,eta)(r^(7/2)/n), while its diagonal
correction is O_(v,eta)(r^2/n). Both are o(r). The partition logarithm
has the same uniform control. Thus the repaired mixture has divergence

    <=-(q/2)sum_j pi_j log(v_j)+o(r)<=B_theta r             (A6)

for one fixed B_theta and all sufficiently large n, uniformly over mu.
The exact repair has bounded divergence and vanishing mass, as in
Section 3. In particular the law is not selected from an unbounded
information class.

For fixed n take the compact convex set of probability laws on the cube
that are centered, have covariance exactly I, satisfy
D(nu||uniform)<=B_theta r, and obey the common K_theta MGF bound.
The MGF conditions are closed linear inequalities, preserved by mixing.
Relative entropy is continuous on this finite
simplex, with 0 log 0=0, and its sublevel set is convex. The payoff
E_mu E_nu|h.x|/sqrt(n) is bilinear. Finite-dimensional minimax applied
to this set and the full probability simplex on C_n therefore converts
(A5) into ONE law satisfying the uniform query bound.

Minimax need not preserve strict positivity of every atom. Add a mixture
of weight n^(-1) of uniform independent signs afterward. This preserves
centering and exact isotropy, does not increase the entropy budget, and
changes each normalized response by at most n^(-1), since both laws have
query second moment one. The spare margin above absorbs that change.
This proves (A2), including full physical support.

### A.3. Geometric meaning and limitation

If C_n lies in a linear space of dimension at most r_n, then s_r(C_n)=1.
Such a code admits arbitrarily small FIXED uniform response with exact
isotropy and an O_theta(r_n) information budget, even if its originally
presented spanning basis had poorly localized rows. The retained large
covariance modes in (A3) are automatically diffuse.

Conversely, if C_n admits a law mu with covariance norm at most L0, then
s_r(C_n)<=L0 r_n/n. The guarantee (A2) is then asymptotically the ordinary
kappa bound. This is a limitation of this covariance-capture mechanism,
not a lower bound against arbitrary non-Gaussian sign laws; the finite
ground-code games earlier in the campaign already distinguish those
claims. No covariance-capture estimate for actual minimizing nearcodes
is supplied by this appendix.
