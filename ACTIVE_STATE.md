# Active research state

Evidence cutoff: 2026-09-07 closing synthesis-to-principle checkpoint.
The six-hour campaign runs15:24--21:24UTC; final verification is in progress.
Read [STEERING](STEERING.md) and the
[current synthesis](artifacts/principle_director_final_synthesis_2026_09_07.md).
Earlier detailed states remain in Git and the archival ledger.

## Exact problem and verified frontier

For hollow symmetric full signings,
H_A(x)=sum_(i<j) a_ij x_i x_j, Q(A)=max_x |H_A(x)|, M_n=min_A Q(A).
Determine convergence OR nonconvergence of M_n/n^(3/2), without fixing a value.

    .4333221116640807 <= liminf M_n/n^(3/2)
       <= limsup M_n/n^(3/2) <= U0-zeta < .493608094.

U0=[97/20+(24/25)log2-5151/6250]/[(97/10)sqrt(24/25)].
The new selector/local-stability construction gives explicit zeta>0,
but it is extremely small; no changed displayed decimal is claimed.
Convergence remains OPEN. Convergence to1/2 is excluded.

For d=binom(n,2), the AUGMENTED cut code is
C_n^+={(c+b_i+b_j)_(i<j)}. Then Q(a)=d-2dist(a,C_n^+),
M_n=d-2rho(C_n^+). The constant bit is essential for the absolute maximum.

Recorded M_3,...,M_14=(3,4,4,5,9,10,12,13,17,18,20,21).
External witnesses give M15<=27,M16<=30; their global catalogue completeness
was not replayed. New asymptotic theorems do not depend on these finite optima.

## Proof dependencies, not inherited verdicts

Upper: precision-Schur BE<=E; Gaussian-boundary direct stopping with budget
Phi-G; ternary directed certificate E<=-5151/6250; full-spin/both-polarity
counting; terminal-uniform Fock bound; H2/H12 all-order realization.
Fix margin, then finite depth, then large order, then remove margin.
The selector uncertainty + actual local-stability refinement yields zeta.
Original certificate replay:70 digits,86,041 boxes,43,021 leaves, PASS.

Lower: marked two-Gaussian response, covariance/tree estimates, nonlinear
replacement, exact finite policy certificate and spectral deletion.
Finite approximations precede order; deletion cutoff is removed last.
See decisive_audit_fresh_full_lower_chain_2026_09_07.md.
The same c* holds for minimum half-range and fixed bounded-amplitude
row-energy-regular weighted matrices; unbounded amplitudes are separate.
Square full-sign bilinear minimum is at least2c* asymptotically.

## New positive organizing principle and actual construction

Retain physical rare-packet entropy BEFORE separating algebraically linked
response channels. A uniform rare-source rate-distortion theorem and exact
Hadamard clipping inequality make this quantitative, not just terminology.

Safe same-law model: Sylvester L32, one hole per node group, marked constant
column, independent recursive children, exact balanced-column sign repair.
For its actual sign-block bulk B, q/m->p=31/32, N=mq:

- Balanced-face cap <.497761196 N^(3/2), including selector entropy.
- Near-constant variance coefficient sqrt(p)/2+o_r(1)<.492125493.
  One common temperature covers constant fibres and densities [r^alpha,r]
  for alpha<=31/16; at coefficient .499 it covers [r^3,r].
- A verified mixed-profile entropy cone permits balanced fibres to pay
  the other rows' costs, including arbitrarily many tiny density scales.
- Stratified restricted-Gram conditioning controls EVERY outer seed and
  all biases when o(m/log^4 m) fibres are nonconstant, at relative
  coefficient sqrt(p)/2+o(1). This modifies the selector law by a uniformly
  bounded conditional density and preserves all prior strict row exponents.
  Elementary unconditioned coherence also covers o(sqrt(m/log m)).
- Actual fair-edge local stability excludes small-minority local maxima
  unless almost all relevant minority mass is coherent or nearly balanced.

These regions do NOT yet cover every mixed profile. Bulk diagonal zeros
are explicit; no full minimax upper bound or ensemble optimality is inferred.
See principle_director_joint_packet_audit_2026_09_07.md and the synthesis.

## New statement about actual liminf families

Hierarchy-aware Gram--Schmidt augmentation preserves all old edges.
If e(eta) is normalized gamma2-square of the FULL absolute near-level set,
every liminf-realizing family satisfies

    limsup_eta e(eta)h(min(e(eta),1/2))/eta >=3c_inf/(40 phi^5).

This is a real necessary condition, not just a cleaner name for convergence.
It does not imply that the near-level code has low complexity.

## Strongest scoped boundaries

1. Child reversal forces max_epsilon|a+b+epsilon c|=|a+b|+|c|.
2. Local iid profiles coexist with preserved leading cap; random restrictions
   need not transfer favorable value. Rare chosen restrictions remain open.
3. Independent-port pressure AND minimum defect are seed-universal at
   Q=o(m^2/log m), sharply. This includes bounded-cap exact seeds but not
   correlated placement or global parent rewrites. Querywise is not minmax.
4. Every specified independent fixed-k orthogonal-frame law has actual cap
   >.50717387; rare/correlated frames and growing physical fibres are distinct.
5. Exact block balance plus strict balanced-face cap does not control mixed
   words: actual macroscopic balanced scars falsify that implication.
6. Unstratified scalar E has rare coefficient1/(2sqrt p); joint packet
   stratification escapes it. This is a certificate comparison, not cap lower.
7. Naive tensor/DC-spine packet iteration has explicit temperature conflicts.
   It does not prohibit all adaptive hierarchies.
8. H24 output-permuting automorphisms do not close source-law mixing.
   H32 translations preserve each fixed channel and give the safe proof.
9. Spectral/subspace residual certificates, planted near-minimizer failures
   and particular algebraic floors must not be generalized to all selectable
   exact children or unrestricted global rewriting.

## Exact remaining original-value gap

For selectable exact children at comparable N=m+n, an actual full-sign
replacement of the row-square-normalized weighted block diagonal with
cap loss O(N^(3/2-delta)) would give summable almost-subadditivity for
u_n=M_n/sqrt(n-1), hence convergence. This is one sufficient theorem,
not a mandatory architecture and NOT proved here.

The new construction still needs complete intermediate/coherent multiscale
control AND a separate favorable-seed-value/all-order transfer theorem.
A genuine nonconvergence proof instead needs two separated infinite
subsequences. No route falsifier counts as that result.

## Preservation

Follow README: never /tmp; preserve drafts, failures, code and outputs
in tracked dated archives. Canonical proofs stay in place. Environments,
caches, credentials and reviewed reproducible builds are excluded.
See research_archive/README.md for hashes, exclusions and dependency recovery.
