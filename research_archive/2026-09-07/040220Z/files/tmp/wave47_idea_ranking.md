# Wave 47 independent idea ranking

Evidence cutoff: ledger §10.99 and STEERING.md after Wave 46.

1. **Polyhedral tight-decomposition separation.**  Rewrite failure of
   (10.1203) as a finite zero-sum/block-game dual certificate.  Test whether
   the active-state optimality conditions forced by exact signing minimality
   contradict that certificate.  This is exact, strongly finitely supported,
   and admits a clean proof or counterexample.
2. **Complement slack-or-improving-replacement dichotomy.**  At the scalar
   optimizer, assume the high-slack conditional fraction (10.1191) is too
   small.  Use the resulting concentration of active constraints near the
   threshold to seek one edge or block flip that lowers every relevant
   parent ground state, contradicting exact minimality.  This directly attacks
   the leading route without geometric Harnack.
3. **Active-witness coarea alignment.**  Couple the principal-cap subgradient
   attaining Q(A[S]) with the fixed center z in D_z(S).  Seek cancellation in
   Johnson replacements strong enough for (10.1198), or construct an actual
   minimizer fiber showing that such alignment is false.
4. **Nested maximal-principal chain.**  Test whether maximal-norm selectors
   admit a nested chain carrying one parent ground restriction.  A greedoid-
   or valuated-exchange property would prove tight decomposition for every m;
   an exact finite failure would sharply prune this strengthening.
5. **Scalar Pareto-curve interpolation.**  Vary the price lambda and track
   derivatives/subgradients of min_d[-log alpha_d+lambda R_2(d)].  Try to show
   that some Pareto transition has both project mass and enough active slack,
   or give a scalable Pareto curve that avoids the desired quadrant.
6. **Posterior-entropy bound on high-side fan-in.**  Starting from
   sum_i r_i=n-m and the exact exclusion log ratio (10.1205), replace nodewise
   degree by a weighted entropy/Fisher-information charge that counts a high
   endpoint only once.  It must beat the tropical-star obstruction using
   minimizer structure.
7. **High Johnson-level recovery.**  Instead of forcing small boundary in
   (10.1195), lower-bound the nonnegative correction Xi_ell from threshold
   geometry of the hard fibers.  Quantify exactly how much high-level mass
   would rescue the negative nearest-core examples.
8. **Larger exact falsifier search.**  Encode exact minimizer status plus
   failure of (10.1203), or failure of complement high slack, as SAT/MILP with
   switching normalization.  Target order nine or structured blow-ups.  This
   is evidence/falsification, not an asymptotic proof.
9. **Profile-conditioned bare-tail count.**  Count arbitrary favorable cuts
   after conditioning on their row norm and the entire principal-norm
   profile, using dependent quadratic-form small-ball estimates.  This avoids
   complement average incidence but must exhibit an exp{-O(H)} lower tail.
10. **Signed multiscale Green-kernel extraction.**  Solve a Poisson equation
    for clipped deficit on the Johnson graph and ask whether exact minimality
    controls the signed multiscale flux even though nonnegative kernel
    mixtures cannot improve extraction.  Positivity loss is the primary
    falsification test.

Selected independent assignments:

- Wave 47A: idea 2, complement slack-or-improving-replacement dichotomy.
- Wave 47B: idea 3, active-witness coarea alignment.
- Wave 47C: ideas 1 and 4, tight-decomposition separation and nested-chain
  strengthening.

Ideas 5--10 remain ranked alternatives; none is being treated as a directive
or as proved mathematics.
