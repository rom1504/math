# Wave 50 independent idea ranking

Evidence basis: README, the Wave 49 Updated frontier, and STEERING at the
Wave 50 boundary.  This ranking was made before assigning agents.

1. **Annealed completion-CDF lower bound from exact optimality.**  Starting
   with (10.1232), exploit the fact that every full cut has energy at most the
   exact optimum to force a non-negligible average lower tail after averaging
   over `S` and the retained coordinates `a`.  Target
   `-log Z_t = O(n^(3/4-c))` for `t = O(n^(3/2-c))`, or a precise family that
   falsifies it.
2. **Deficit-stratified annealed transport.**  Decompose full grounds by the
   retained block deficit and norm of the linear completion field in
   (10.1231), then seek a scale on which both strata have enough mass.  This
   is testable and keeps the tolerance discarded by complement incidence.
3. **Adaptive-scale triple retention by spectral averaging.**  Average the
   exact cubic criterion (10.1234) over Johnson distances and prove that at
   least one balanced nonlocal scale with spectral gap ratio bounded below
   meets (10.1235), using exact-minimizer structure rather than a generic
   Boolean-set assertion.
4. **Partial-completion box mass from minimizer-specific column geometry.**
   Attack (10.1237) directly, seeking a signing of the omitted columns that
   cancels the child residual below `O(n^(9/4-c))`; quantify exactly what
   correlation supplied by full optimality is needed beyond uniform
   completion.
5. **Joint scale selection for box mass and coarea.**  Put the box-success
   mass and the cubic-retention surplus into one weighted average over
   nonlocal Johnson scales, so the same scale—not two unrelated existential
   choices—satisfies both clauses.
6. **Global all-maximal-bad anti-migration.**  Sum the all-state Bellman
   identity (10.1240) over every maximal selector with weights adapted to its
   state-dependent margin; test whether the assumption that all selectors
   are bad forces an impossible circulation or energy drift.
7. **Exact falsification search for tight decomposition.**  Enumerate the
   smallest exact minimizers beyond the known cases and solve the full
   all-selector state system, looking either for an all-bad configuration or
   for a dual certificate suggesting the missing global theorem.
8. **Conditional Paley--Zygmund for completion profiles.**  Use first and
   second moments of the completion CDF in (10.1232), conditioned on a
   favorable deficit/norm stratum, to obtain a lower rather than upper tail;
   explicitly test whether fourth moments remain below the required
   `exp(O(n^(3/4-c)))` scale.
9. **Replacement contradiction from successful nonlocal coarea.**  Assume
   every admissible nonlocal scale fails triple retention, combine the exact
   failures into a signed replacement flow, and try to construct a full
   signing above `q_n`, contradicting exact maximality.
10. **Minimizer-specialized partial-coloring theorem.**  Search for a vector
    balancing/partial-coloring statement whose hypotheses are forced by an
    exact quadratic minimizer and whose discrepancy exponent reaches the box
    target; generic law-free discrepancy bounds are known to be too weak.

Selected independent attacks: 1--2 as one annealed-tail program; 3--5 as one
box/coarea program; and 6--7 as one global tight-decomposition program.
