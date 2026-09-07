# Wave 28 plan

Started: 2026-07-29T23:11:40Z

Read before ranking: README.md, STEERING.md (evidence cutoff Wave 25),
and the Wave 27 Updated frontier in ledger.md.

## Ranked ideas

1. **Row-good random-coset hit probability.**  For a balanced partition,
   study the law of
   `Q(P^T D H_S D P)` under random `D`, jointly with the row-good
   event from (10.833).  Seek the exact lower bound (10.838) with
   `-\log\delta=O(n^{3/4-c})`; if impossible from generic matrix data,
   isolate a sharp scale/structure obstruction rather than invoking a
   scalar sign tail.

2. **High-excess structural partition.**  Starting from an optimizer of
   `H_S` or an exact child ground, build an `A`-adapted partition and
   diagonal whose block quotient retains the excess while the complete
   coset is row-good.  Use child fields, column-vector balancing, or
   discrepancy; retain the global ratio-window quantifiers.

3. **Child-label entropy/collision dichotomy.**  Study the family of exact
   child grounds on selectors with `Y_A(S)>t`.  Prove either that it has
   an `exp(O(n^{3/4-c}))` representative/block-coset cover or that the
   restriction excess is already power-saving.  Candidate tools are
   Shearer entropy, VC/Sauer--Shelah, exchange identities, and exact
   signing minimality.

4. **Random hashing around a planted completion.**  Prove a CountSketch-like
   theorem: if a favorable completion `x` has low `R_2(x)`, a random
   partition with `D=\operatorname{diag}x` produces a whole row-good
   coset containing `x`, with a bound in terms of
   `R_2(x),\|A\|_{\rm op}`, and `n^2/k`.  Determine whether this can
   reduce shared completion to a finite-output map.

5. **Nonlinear row-penalized Laplace cover.**  Replace the circular linear
   mean dual (10.823) by a log-sum-exp/positive-part functional over
   compressed cosets.  Try to lower-bound a row-conditioned partition
   function whose layer cake gives (10.838), without assuming the mean
   restriction lemma.

6. **Multiscale hierarchical cosets.**  Use a tree of partitions with total
   free-sign dimension `O(n^{3/4-c})\).  Let each selector stop at a
   scale where its child-ground variation is compressible; prove the union
   remains row-good and quantify the accumulated alignment error.

7. **Weak-regularity/cut-decomposition at the minimizer scale.**  Test
   whether a Frieze--Kannan or Grothendieck decomposition of `A` can
   approximate every
   `H_S=\operatorname{diag}(\xi_S)A\operatorname{diag}(\xi_S)-p_2A`
   at error `O(n^{3/2-c})` using only the allowed entropy.  Record the
   exact exponent wall if the usual energy increment is too expensive.

8. **Shared low-row completion.**  Turn the per-selector completion formulas
   (10.826)--(10.830) into a stochastic/canonical completion map with
   mutual information or output entropy `O(n^{3/4-c})\).  Exploit that
   every completion of an exact child ground is loss-favorable, but do not
   assume positive variance on every fiber.

9. **Matched parent transport through a compressed prior.**  Use a
   row-regular block-coset law as the common parent prior in
   (10.799)--(10.800), and look for a data-processing comparison that
   controls both the barycenter entropy and adjacent finite-measure
   Hellinger term.  It must be minimizer-specific and avoid the generic
   LSI walls.

10. **Temperature-adaptive signed deck.**  Explore whether the finite
    low-temperature/high-temperature crossover in (10.844) can yield a
    weaker summable restriction comparison with `\beta=\beta_r`, or
    prove that the temperature error makes this incapable of convergence.
    Do not relabel it as fixed-temperature constant shortfall.

## Selected agent routes

- Route 1: row-good random-coset hit probability.
- Route 2: high-excess structural partition.
- Route 8: shared low-row completion and collision.

The main agent independently takes Route 4, with Route 5 as its immediate
fallback.  Default agent deadline is twenty minutes, with checks near five
minutes and extensions only for a concrete near-complete derivation.
