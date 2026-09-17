# Active research state

Evidence cutoff: 2026-09-17 sixth paper-combination checkpoint.
Campaign ACTIVE16:42:21--22:42:21UTC; stop for assessment at its endpoint.
Read [STEERING](STEERING.md) and the [current campaign](artifacts/paper_portfolio_campaign_2026_09_17.md).
The long ledger and dated archives retain prior proofs, drafts and failures.

## Exact problem and reported frontier

For hollow symmetric full signings,
H_A(x)=sum_(i<j)a_ij x_i x_j, Q(A)=max_x|H_A(x)|, M_n=min_A Q(A).
Determine convergence OR nonconvergence of M_n/n^(3/2), without fixing its value.

    .4333221116640807 <= liminf M_n/n^(3/2)
       <= limsup M_n/n^(3/2) <= U0-zeta < .493608094.

U0=[97/20+(24/25)log2-5151/6250]/[(97/10)sqrt(24/25)];
zeta>0 is explicit but tiny. Convergence remains OPEN; convergence to1/2
is excluded by the preserved strict upper theorem. No original-bound
improvement is claimed from this paper campaign.

For d=binom(n,2), the AUGMENTED cut code is
C_n^+={(c+b_i+b_j)_(i<j)}. Then Q(a)=d-2dist(a,C_n^+),
M_n=d-2rho(C_n^+). The constant bit is essential for the absolute maximum.

Recorded M_3,...,M_14=(3,4,4,5,9,10,12,13,17,18,20,21).
Forty-two new physical-column response games on stored orders4--14 have
rational primal/dual certificates. Root replayed every physical column
and cap without a solver; global optimality of the stored signings is
imported, not newly re-proved. Finite games are not asymptotic evidence.

## Original proof dependencies

Lower: marked Gaussian response, covariance/tree estimates, nonlinear
replacement, the finite policy certificate, and spectral deletion.
The half-range extension was reread and the exact policy replayed during
this campaign. Finite approximations precede order; cutoff removal is last.
See decisive_audit_fresh_full_lower_chain_2026_09_07.md and
decisive_audit_certified_minimum_width_lower_2026_09_07.md in artifacts/.

Upper: precision Schur, Gaussian direct stopping, ternary certificate,
full-spin/both-polarity counting, terminal-uniform Fock bound, all-order
realization, then selector/local-stability refinement. Fix margin and
finite depth before order, then remove the margin. The current campaign
preserves rather than re-certifies the entire upper chain.
See artifacts/principle_director_final_synthesis_2026_09_07.md.

## Strongest new theorem about actual low-cap signings

For EVERY sufficiently large full signing with Q(A)<=.5n^(3/2), ONE
centered, exactly isotropic, (3/2)-subGaussian physical law satisfies

    sup_(|H_A(x)|>=.30n^(3/2)) E|h.x|/sqrt(n)
        <=kappa-2^(-67),       kappa=sqrt(2/pi).

No operator-norm or dual-covariance hypothesis is imposed. The proof
combines signed covariance, Grothendieck coordinate localization and the
universal full-sign half-range budget. Three independent reconstructions
check all quantifiers. This discount does not reach the roughly3c/2
parent slope and does not control the full old/new-spin maximum.
Canonical: artifacts/paper_director_low_cap_uniform_response_2026_09_17.md.

## Strongest reusable paper combinations

- Same-order cloned-block surgery regularizes ANY prescribed
  T=o(n^(3/2)) nearlevel window into exp(o(n)) witnesses at o(n^(3/2))
  cap cost. Mixed-tail chaining proves actual independent sign-flip
  stability and a random Hamming-path statement. Its certified error
  exceeds the mean contraction gain; no favorable noise theorem follows.

- Diffuse feature spaces of rank r=o(sqrt(n)) admit actual non-Gaussian
  sign laws with exact isotropy, uniform Boolean scalar responses,
  fixed subGaussian control and O(r) relative entropy. Gaussian transforms,
  positive convex auxiliary measures, operator covariance and exact
  covariance repair are all independently reconstructed. Adaptive
  covariance capture follows by compact minimax. An efficient cold
  sampler or minimax algorithm has NOT been proved.
  Canonical: paper_director_growing_rank_variance_realization_2026_09_17.md.

- On block-constant query codes, the optimal physical information cost
  is r log(1/epsilon)+O(r), for every rank at large enough block size.
  Fixed tail budget K>1 has sharp growing-rank floor kappa/sqrt(K);
  approaching it within delta costs
  (1-1/K+o(1))r log(1/delta). Physical common-phase Gibbs laws attain
  the coefficient. Independent phases with exactly the same block
  marginals and covariance fail. This is a genuine composition benchmark,
  not an assertion that original minimizing codes have block geometry.
  Canonical: paper_director_block_response_information_tail_frontier_2026_09_17.md.

External novelty is unestablished; classical ingredients and imported
hypotheses are explicitly attributed. New here means proved in this
repository, not a priority claim.

## Strongest scoped obstructions

1. Child reversal forces max_epsilon|a+b+epsilon c|=|a+b|+|c|.
   No independent cancellation argument can violate this identity.
2. The Gaussian latent band[1/2,3/2] has sharp isotropic-mixture response
   floor .786393873897..., above every useful slope below .75.
   Non-Gaussian laws escape the class; this is not a sign-law impossibility.
3. Query covariance norm<=L forces ANY physical discount Delta to cost
   KL>=Delta^2/(64L)-2. Positive Walsh densities of degree o(n) therefore
   cannot achieve fixed normalized discount on bounded-covariance duals.
   No such dual is proved for every minimizing ground code.
4. Fixed bounded-Lipschitz functions of finitely many actual central
   energies leave scalar responses at kappa+O(n^(-1/4)). Extensive
   Gibbs tilts and rare conditioning are not covered.
5. Whole nearlevel entropy is at least order eta log(1/eta); finite
   full-code entropy slope and response^2=o(eta) shortcuts are vacuous.
   Anchored centers can still be useful if their residuals are paid.
6. Original local-profile, independent-port, fixed-frame, scalar Schur
   and algebraic-family failures retain their precise quantified scopes.
   They do not rule out selectable exact children or global rewrites.
7. The earlier BH proof is valid, with a checked analytic refinement,
   but powers through o(n^(3/4)) miss bounded-cap extrema. Linear-degree
   optimized signed moments remain missing. Generic power variants are
   not revived by the current paper campaign.

## Exact remaining original-value gap and stopping condition

No theorem shows that actual minimizing high-energy codes possess the
low-rank covariance capture needed by the new realization. Their diffuse
signed-covariance component can have macroscopic rank. A useful physical
mean response also still needs full all-energy fluctuation/escape control.

One sufficient original theorem is full-sign replacement of the
row-square-normalized two-child weighted block diagonal with cap loss
O(N^(3/2-delta)); archived almost-subadditivity would then give convergence.
This is not proved or assumed to be the only possible architecture.

At22:42UTC finish preservation, synthesis and assessment; do not start
another autonomous campaign. Every counterexample above is scoped, not
a proof of nonconvergence. A real nonconvergence proof needs separated
infinite subsequences.

Never use /tmp. Preserve research-bearing scratch files, including failed
derivations, at substantive checkpoints. See research_archive/README.md
for manifests, exclusions, verified hashes and dependency restoration.
