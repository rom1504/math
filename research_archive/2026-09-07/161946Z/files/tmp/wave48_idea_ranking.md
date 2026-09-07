# Wave 48 independent idea ranking

Evidence used after generation: Wave 47 `Updated frontier` (§10.100.4) and
`STEERING.md` at the Wave 47 cutoff.  This is a scratch ranking, not a claim
that any proposed lemma is true.

1. **Fixed-density omitted-block excess dichotomy.**  For complement columns
   `T=S^c` with `|T|=rho n`, combine the exact cap
   `Gamma_S <= 2 Q(A[T])` with the scalar-column condition to seek a dichotomy:
   either a project-mass fraction of active `T` has `Q(A[T]) >= c q_n`, or the
   aggregate low-slack responses yield an integral congestion contradiction.
   Success means proving the mass/high-slack pair in the current steering
   document; failure means an exact-minimizer family with project column mass
   but exponentially too little principal excess.  Highest relevance and a
   newly sharpened arithmetic input, although the integral response gap is a
   serious obstruction.

2. **Row-truncated incidence-size-biased coarea.**  For the hard common-core
   incidence graph, take the ground-lift law proportional to `1_C(z) a_z`,
   where `C` is the project-row truncation, and prove
   `E[1_C a B] <= (1-lambda)(1-eta) E[1_C a^2]` at a nonlocal core scale while
   retaining project mass.  This is an exact, falsifiable inequality whose
   conclusion gives hard extraction.  Finite nearest-core failures demand a
   genuinely nonlocal or truncated mechanism.

3. **Cross-block compatibility of tight-decomposition duals.**  Attach the
   one-block zero dual law to every maximal selector and use overlaps to prove
   that no jointly realizable family can satisfy all edge sign inequalities
   unless some selector already has a ground child.  A concrete sufficient
   theorem is a positive lower bound for the sum of fractional block margins
   over an overlap cover.  This targets exactly the gap left by the abstract
   antipodal obstruction and has unusually strong finite evidence.

4. **Fixed-density scalar-temperature Pareto interpolation.**  Study the
   log-mass/conditional-slack Pareto curve of complement columns as the scalar
   multiplier `lambda` varies near `n^{-3/2}`.  Prove that its integrated
   derivative cannot stay below both required thresholds throughout a compact
   density window.  This could produce one usable multiplier without a
   pointwise Harnack theorem.  It is distinct from geometric column Harnack but
   must be rejected if the derivative identity still averages over columns
   rather than locating one.

5. **Profile-conditioned direct bare-tail counting.**  Partition restrictions
   by a coarse vector of principal norms and response overlaps, then apply a
   second-moment or entropy argument inside one profile cell to obtain an
   arbitrary cut satisfying the bare mass/row objective.  The exact target is
   a cell with `-log mass + n^{-3/2} R_2 = O(n^{3/4-c})`.  This avoids the
   special complement arithmetic and remains the principal fallback.

6. **Multi-port regret transport.**  Generalize the exact adjacent-port
   decomposition from one added vertex to a block of `k` ports and seek a
   telescoping one-sided inequality in which inactive-port regrets cancel in
   expectation under a size-biased center law.  A bound of total positive
   regret `<= (1-lambda-o(1))` times total escape would prove hard extraction.
   It is more flexible than uniform coarea, but risks merely repackaging the
   currently uncontrolled one-sided growth.

7. **High Johnson-level correction under incidence bias.**  Expand both degree
   and escape functions on the relevant slice and prove that project-row
   truncation suppresses the high-level component enough for their covariance
   to be negative at the required scale.  A concrete advance is an
   `L2(nu)` bound on levels above `L_0` of size `e^{-c L_0}` relative to
   `E[a^2]`.  This is testable, but earlier local harmonic/Poincare routes warn
   that low-side smoothing alone cannot control high-side fan-in.

8. **All-size common-parent compactness without nested selectors.**  Use the
   finite property that one parent ground works separately at every size,
   forget selector nesting, and seek a Helly/compactness theorem for the family
   of size-indexed ground-child constraints.  An exact theorem giving one
   parent with a ground restriction at a prescribed geometric sequence of
   sizes would already imply the restriction iteration.  It is stronger than
   current evidence and must not rely on the falsified greedoid or nesting
   claims.

9. **Coupled discrepancy rounding of block dual laws.**  Instead of rounding
   each fractional replacement independently, round all overlapping block
   responses at once with a signing/discrepancy theorem.  A useful bound would
   make every active-state margin positive up to error `o(q_n)` whenever the
   aggregate fractional margin is `Omega(q_n)`.  The adaptive-witness gap may
   still be fatal; the only credible version must exploit cross-block
   correlations specific to one parent signing.

10. **Exact search for jointly realizable obstructions.**  Formulate a SAT/MILP
    or exhaustive switching-class search at order 9 for (a) failure of the
    all-size common-parent property, (b) a jointly consistent antipodal block
    dual, or (c) fixed-density scalar columns with mass but no high slack.
    This cannot prove the asymptotic theorem by itself, but a small genuine
    counterexample would retire a leading conjecture and identify the missing
    constraint.  Ranked last because computation is diagnostic rather than a
    convergence mechanism.

Selected independent attacks: 1, 2, and 3.  Idea 5 remains the leading
fallback if the structured complement attack fails; ideas 4, 6, 7, and 9 are
possible mechanisms within the selected routes rather than separate agents.
