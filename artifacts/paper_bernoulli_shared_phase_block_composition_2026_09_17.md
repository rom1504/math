# Identical physical block marginals can have different leading responses

2026-09-17. Independent finite reconstruction of the director's
[exact shared-phase block construction](paper_director_block_code_response_information_2026_09_17.md),
Appendix A, with the product-marginal comparison and exact total
correlation. The example concerns declared block-constant queries,
not unknown actual-minimizer ground codes.

## 1. Fix the physical block size before taking a block-count limit

Fix0<epsilon<1/4 and ANY integer b>=4096 epsilon^(-2). For a uniform
b-sign word H let S=b^(-1/2)sum_i H_i. Write

    a=epsilon/4,     t=4/epsilon,
    P_0=Law(H | |S|<=a),
    P_1=Law(H | t<=|S|<=t+1),
    v_j=E_(P_j) S^2,
    pi=(1-v_0)/(v_1-v_0).

Both conditioning events have strictly positive probability for every
such integer b; the finite lattice and binomial bounds are paid in the
linked appendix. No b-to-infinity limit is needed here. The bands are
disjoint and each law is invariant under all coordinate permutations
and global sign reversal. In particular

    v_0<=a^2<1<t^2<=v_1<=(t+1)^2,
    0<pi<1,       (1-pi)v_0+pi v_1=1,
    pi<=2/t^2.

Put q=(1-pi)P_0+pi P_1. This SINGLE-BLOCK physical law is centered and
EXACTLY isotropic: by exchangeability, for i!=j,

    E_q H_i H_j=(E_q S^2-1)/(b-1)=0.              (1)

## 2. Shared phase versus independent copies of the very same marginal

For every integer r>=1, form n=br physical coordinates in r blocks.
Consider the two genuine sign laws

    rho_r=(1-pi)P_0^(tensor r)+pi P_1^(tensor r),
    eta_r=q^(tensor r).                           (2)

They have EXACTLY the same entire physical distribution q on every
block. Both are centered and have full global covariance I_n. For
rho_r, cross-block covariances vanish conditional on the shared label
because the conditional blocks are independent and centered. For
eta_r this is ordinary independence. Thus neither global covariance
nor all the block marginals distinguishes their scalar response.

Let C_(b,r) be all2^r Boolean words constant on each block. Its query
with block signs z_j has normalized field

    h.x/sqrt(br)=sum_j z_j S_j/sqrt(r).

For rho_r, conditional independence and centering give second moment
v_j in phase j. Consequently, for EVERY r and EVERY declared query,

    E_(rho_r)|h.x|/sqrt(br)
       <=(1-pi)sqrt(v_0)+pi sqrt(v_1)
       <=a+2/t=3epsilon/4.                        (3)

For eta_r, the S_j are iid, centered, symmetric, have variance one,
and satisfy |S_j|<=t+1. Symmetry means that the scalar field has EXACTLY
the same distribution for all z, not merely a union-bound estimate.
The scalar Wasserstein estimate reconstructed in the
[hot-feature proof](paper_bernoulli_hot_feature_laws_2026_09_17.md),
Section4, yields the finite uniform bound

    |E_(eta_r)|h.x|/sqrt(br)-kappa|
       <=3 E_q|S|^3/sqrt(r)
       <=3(t+1)/sqrt(r).                          (4)

In particular, keeping epsilon and b fixed and then taking r to
infinity, EVERY block-constant query has normalized response tending
to kappa under eta_r, while rho_r retains(3). This is a leading-order
physical response separation with identical block laws and covariance.

## 3. The common label has one-bit entropy but extensive total correlation

Let h_2(pi)=-pi log pi-(1-pi)log(1-pi), in natural units. Because the
two bands are disjoint, observing ANY single physical block reveals
the hot/cold label exactly. Therefore the entropy identities are

    H(q)=h_2(pi)+(1-pi)H(P_0)+pi H(P_1),
    H(rho_r)=h_2(pi)+r[(1-pi)H(P_0)+pi H(P_1)].

It follows EXACTLY that

    D(rho_r || q^(tensor r))=(r-1)h_2(pi).          (5)

Equivalently, with U_b the uniform physical block law,

    D(rho_r || U_b^(tensor r))
       =r D(q||U_b)+(r-1)h_2(pi).                 (6)

The shared latent label has entropy at most log2, independent of r.
That does not make its correlation cost vanish: its repeated imprint
across r disjoint physical blocks creates the extensive total
correlation in(5). Replacing it by independent labels removes precisely
that dependence and leads to the central-limit response in(4).

The linked finite conditioning estimates give

    D(rho_r || uniform cube)<=r log(1/epsilon)+23r.

This is a paid physical information cost, not a free hot variance
reservoir or a comparison of merely formal Gaussian covariances.

## 4. Full-support version with exactly matching repaired marginals

If full support is required, fix delta=epsilon/8 and define

    rho'_r=(1-delta)rho_r+delta U_b^(tensor r),
    q'=(1-delta)q+delta U_b,
    eta'_r=(q')^(tensor r).                        (7)

Again rho'_r and eta'_r have exactly the SAME block marginal q', are
centered, have covariance I_n, and now assign positive probability to
EVERY physical sign word. The shared law has normalized response
at most7epsilon/8, by(3) and the uniform response bound at most one.

Under q', the block magnetization remains centered, symmetric and
variance one. Also

    E_(q')|S|^3
       <=(1-delta)(t+1)+delta sqrt(3),

because E_(U_b) S^4<=3. Hence the same variance-uniform scalar estimate
shows that eta'_r has response tending to kappa at every declared
query, with error at most three times this third-moment bound/sqrt(r).
This rate is independent of b once the finite band conditions hold.

The exact formula(5) is stated BEFORE this uniform-support repair.
After repair the three phase distributions overlap, so it must not be
reused unchanged. Nevertheless total correlation is still extensive:
if J is the three-valued common label and I_1=I(J;one block)>0, then

    D(rho'_r || (q')^(tensor r))
       =r I_1-I(J;all blocks),
    r I_1-H(J)<=D(rho'_r || (q')^(tensor r))<=r I_1. (8)

All laws, mixture weights and entropies in this example are finite.
The only limit used to separate the product response is the block
count r tending to infinity after epsilon,b are fixed. The stronger
uniform bound in(4) shows there is no hidden Gaussian feature limit.

The localization track independently read and reconstructed the entire
finite construction and comparison: PASS, including both entropy
identities and the use of the repaired marginal q' in the full-support
version. Its separate conditional-slice concentration result can supply
a dimension-free subGaussian proxy at each fixed epsilon; that extra
input is not needed for the exact response separation above.

Finite replay is preserved in
`computations/paper_bernoulli_2026_09_17_shared_phase_composition.py` and
`tmp/paper_portfolio_2026_09_17/bernoulli/shared_phase_composition.json`.
Its independent b=8 zero/ferromagnetic-slice toy enumerates ALL65,536
physical order16 words and checks centering, covariance I and identical
entire block marginals by exact integers, both before and after a
full-support repair. The unrepaired normalized responses are exactly
1/4 versus15/32. A compressed scalar FFT replay through4,096 blocks
approaches kappa/sqrt(8) versus kappa. The toy checks the composition
identities, not the separate epsilon-band size condition in Section1;
its limiting FFT values are numerical diagnostics, not substituted for
the finite and asymptotic proofs above.
