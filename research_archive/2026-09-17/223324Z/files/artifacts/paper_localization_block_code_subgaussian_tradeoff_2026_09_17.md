# Sharp large-rank subGaussian cost of protecting a block code

2026-09-17. Independently obtained by the localization and discrepancy
tracks. This is a declared-query benchmark, not an actual-minimizer
geometry assertion. It complements the sharp relative-entropy price in
[the block-code theorem](paper_director_block_code_response_information_2026_09_17.md).

## 1. Universal lower bound

Partition n=br coordinates into r equal blocks and let C_(b,r) consist
of all block-constant Boolean words. Suppose a physical law nu is exactly
isotropic and has linear subGaussian proxy K I, where K is fixed. Put
S_j=b^(-1/2)sum_(i in block j)h_i. Then S is exactly isotropic in R^r
and K-subGaussian. In particular K>=1. Uniformly in b and nu,

    max_(x in C_(b,r)) E_nu |h.x|/sqrt(n)
       >=kappa/sqrt(K)-o_K(1),     r->infinity.      (1)

Proof. Gaussian integration and the assumed linear MGF imply, for
0<t<1/(2K),

    E exp(t||S||^2)
      =E_g E exp(sqrt(2t)g.S)
      <=E_g exp(Kt||g||^2)=(1-2Kt)^(-r/2).

Thus for0<d<=1, with R=||S||,

    P(R^2>=Kr(1+d))
       <=exp[-r(d-log(1+d))/2]<=exp(-rd^2/8).

Using the same Chernoff parameter t=d/[2K(1+d)] beyond this threshold
and integrating gives

    E R^2 1_{R^2>Kr(1+d)}
       <=[Kr(1+d)+2K(1+d)/d]exp(-rd^2/8).

Since E R^2=r by EXACT isotropy, R^2<=sqrt(Kr(1+d))R below the
threshold yields

    E R/sqrt(r)>=1/sqrt(K(1+d))-o_K(1)

when d->0 and rd^2 grows sufficiently fast, for example d=r^(-1/4).

For any fixed realization s, a scalar third-moment Gaussian comparison
for independent fair query signs z gives

    E_z |sum_j z_j s_j|>=kappa||s||-3max_j|s_j|.    (2)

Indeed the standard normalized Wasserstein bound costs at most
3sum|s_j|^3/||s||^2<=3max|s_j|; the zero vector is immediate. The
scalar estimate was reconstructed in the
[hot-feature proof](paper_bernoulli_hot_feature_laws_2026_09_17.md), Section4.
Also subGaussianity and the log-sum-exp bound give

    E max_j |S_j|<=sqrt(2K log(2r)).

Average the block queries z, apply (2), and divide by sqrt(r). This
proves (1). No independence of the physical coordinates or block
magnetizations is assumed in this lower bound.

Consequently a large-rank block code with all normalized responses at
most epsilon requires K>=kappa^2/epsilon^2 asymptotically. This is a
different parameter from covariance, which is exactly I throughout.

## 2. Physical sharpness for K>=2

For even b choose the largest positive even integer m_b<=sqrt(Kb),
and put K_b=m_b^2/b. For b sufficiently large,1<=K_b<=K and K_b->K.
Draw one common active label with probability1/K_b. In the inactive
phase, independently sample a balanced uniform slice in every block.
In the active phase, choose independent fair block poles and sample
each block uniformly with magnetization plus or minus m_b accordingly.

Each block has mean zero and magnetization second moment b, while
different blocks have zero cross-correlation. The law is therefore
EXACTLY isotropic. The fixed-slice perpendicular MGF from
[the swap proof](paper_discrepancy_kernel_slice_repair_2026_09_17.md),
Section2, and the bound |S_j|<=sqrt(K_b) give matrix proxy

    2(I-P)+K_b P<=K I,

where P projects onto the block-constant subspace. The inequality uses
K>=2. Every block query has the SAME exact normalized response

    E|epsilon_1+...+epsilon_r|/sqrt(K_b r),          (3)

which tends to kappa/sqrt(K) as b,r both tend to infinity, with no
relation between their rates. Adding any vanishing positive uniform
cube mass gives full support, preserves exact isotropy and the same
K-proxy, and does not change this limiting response.

Thus (1) is sharp for every fixed K>=2 in the physical sign model.
The information-optimal shared-phase construction's proxy O(epsilon^-2)
has the correct dependence on the desired response, not an avoidable
dimension or block-count loss. The scalar, rank-one case is different:
a variance-one1-subGaussian three-point law can have absolute moment
1/sqrt(3), whereas (1) describes the constraint imposed by an entire
growing-rank Boolean query code.

The discrepancy researcher independently read and reconstructed both
sections in full, including the integrated norm tail and finite slice
rounding: PASS. The director's subsequent
[joint information/tail frontier](paper_director_block_response_information_tail_frontier_2026_09_17.md)
uses physical hot/cold Gibbs blocks to extend asymptotic sharpness to
EVERY K>1. It also quantifies the information needed to approach the
floor: per-feature cost has leading coefficient (1-1/K)log(1/delta).
This track independently audited that complete theorem. The explicit
slice construction above remains a separate elementary proof for K>=2.
