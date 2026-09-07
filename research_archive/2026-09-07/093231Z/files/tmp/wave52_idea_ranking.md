# Wave 52 independent idea ranking

Evidence basis: README, the Wave 51 Updated frontier, and STEERING at the
Wave 52 boundary. This ranking was made before assigning agents.

1. **Level-two spectral-excess identity.** Expand `P_2-lambda_2(2)` exactly
   through pair loads of the favorable selector families. Seek a
   minimizer-specific inequality forcing saved positive excess at a low-row
   cap, or prove that the level-two formulation is still degree-tautological.
2. **Adaptive higher-core spectral averaging.** Express the family of excesses
   `P_ell-lambda_2(ell)`, `2<=ell<m`, in the Johnson harmonic weights and
   average over `ell`. Prove that some scale has saved excess unless the
   favorable family has a structural form incompatible with child grounds of
   an exact minimizer.
3. **Selector-exchange overlap from child optimality.** Starting with one
   box-witness pair `(S,y,w)`, exchange ports between `S` and `S^c` and use
   exact child-ground regret identities to force a second favorable selector
   for a common low-row completion. Quantify how many exchanges survive.
4. **Positive signed-shell budget.** Weight the exact cap increments
   `sum_{R_2(z)=R} a_z^2(p_ell(z)-lambda_2)` by a function of the row and seek
   a telescoping identity from `A^2`, child deficits, or full optimality that
   makes some controlled prefix positive at the saved scale.
5. **Row-neighborhood collision amplification.** Combine the outside-flip
   ball around a box witness with exchange maps between nearby selectors;
   test whether double counting forces distinct selector labels to collide at
   one center before the cap exceeds `C(R0+n^2)`.
6. **Rademacher K-functional far completion tail.** Apply a verified reverse
   lower-tail theorem for the cross-linear coefficients
   `beta=2 A[S^c,S]y`. At entropy budget `K=n^(3/4-c)`, the K-functional may
   reach `n^(3/2-c)` once cross variance is above the project-row scale.
   Control the outside quadratic noise separately and audit noncircularity.
7. **Coefficient-profile dichotomy.** Split the cross-linear field into
   concentrated and diffuse coefficients. Fix the largest `K` signs in the
   concentrated case and use a Gaussian-scale reverse tail in the diffuse
   case; make both branches quantitative at the saved entropy scale.
8. **Signed-spectrum far tail for the outside quadratic block.** When the
   cross-linear K-functional is too small, seek a negative Rademacher-chaos
   tail from the eigenvalue/triangle profile of `A[S^c]`. The target must
   contribute outside the `O(n^(3/2-c))` local-margin band and cannot rely on
   moments through order three alone.
9. **Partial-completion box discrepancy via variance dichotomy.** If the cross
   field is below `O(n^(9/4-c))`, combine it with child one-spin stability and
   a partial-coloring extension theorem to bound the internal and completed
   row norm; otherwise feed the high cross variance into idea 6.
10. **Exact finite spectral-shell atlas.** On all stored exact minimizers,
    compute `P_ell-lambda_2` at every row shell and correlate it with pair
    loads, exchange regrets, and cross coefficient K-functionals. Use this
    only to falsify candidate inequalities and choose invariant quantities,
    not as asymptotic evidence.

Selected independent attacks: 1--3 as a direct level-two/higher-core overlap
program; 4--5 as a signed-shell and cap-inflation program; and 6--8 as a
far-negative completion-tail program. Ideas 9--10 support the dichotomies and
falsification audit.
