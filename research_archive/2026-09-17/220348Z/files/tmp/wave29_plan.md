# Wave 29 plan

Started: 2026-07-29T23:44:09Z

Read before ranking: README.md, STEERING.md (Wave 25 evidence cutoff), and
the Wave 28 Updated frontier in ledger.md.

## Ranked ideas

1. **Johnson-coherent low-signature witnesses.**  Choose favorable child
   completions jointly along adjacent selectors, using overlaps of
   near-ground sets and outside freedom so the relative coordinate
   signatures satisfy `K_b(F)<=k`.  Seek an exact Hall/coupling or
   Lipschitz-selection lemma; audit finite minimizers for a compatibility
   obstruction.

2. **Joint aligned-payoff/low-row density.**  For fixed `S`, analyze
   `E_{T,R}={x:|x^T H_Sx|>=T,\ x^TA^2x<=C}`.  Prove that a scalar aligned
   reverse tail at the exact boundary retains comparable mass after the row
   cap, using a joint exponential tilt, conditional Hanson--Wright, or an
   exact derivative identity.  A rigorous counterexample to automatic
   retention is equally useful.

3. **Matched parent Hellinger transport.**  Attack (10.799)--(10.800) for
   the actual likelihood `f=Ca`, not generic cube functions.  Use the
   barycenter identity, parent conditional fields, and minimizer-specific
   complete-signing constraints to compare parent entropy with adjacent
   finite-measure Hellinger transport.

4. **Adapted common-coset minimax.**  Optimize a distribution on row-good
   diagonals against the worst selector hit in (10.838), write its exact
   dual, and ask whether signing minimality separates every bad selector law.
   The functional must be nonlinear in coverage and must not collapse to the
   circular uniform mean dual (10.823).

5. **High-excess selector dichotomy.**  Split selectors by `Y_A(S)`.
   Cover the low-excess part automatically; if the high-excess set has at
   most `exp(O(n^(3/4-c)))` structured types, plant them.  If it is larger,
   seek an entropy/exchange contradiction from the associated child-ground
   family.

6. **Affine closure inside the low-row set.**  Find a block affine subspace
   of the abundant set `{x:R_2(x)<=C}` with dimension
   `Theta(n^(3/4-c)/log n)`, then maximize selector alignment over its
   translates.  Use additive-combinatorial energy or a PSD quadratic
   sublevel theorem; respect the A8 ternary-closure wall.

7. **Signature-field discrepancy.**  Construct signature classes directly
   so the necessary total variation in (10.858) is small, then control the
   full Boolean maximum of the refined quotient.  Candidate tools are
   vector balancing and multicolor discrepancy applied jointly to child
   fields.

8. **Multicoset fractional cover with adapted weights.**  Use the exact
   coset-saturation sets and solve the fractional set-cover dual before
   greedy rounding.  Seek a minimizer-specific lower bound for every
   selector weighting, rather than uniform random diagonals.

9. **Adjacent-ground sheaf obstruction.**  Regard child ground/near-ground
   sets as fibers over the Johnson graph and test whether local overlaps
   have nontrivial holonomy preventing a coherent global section.  A finite
   exact obstruction would sharply delimit idea 1; vanishing holonomy plus
   quantitative overlap could construct the section.

10. **Literature-guided correlated rounding.**  Search primary sources on
    list recovery, agreement theorems, affine extractors, and matrix/vector
    discrepancy for a theorem matching low-signature completion.  Import
    only a result whose hypotheses and exponents can be checked against the
    ledger.

## Selected agent routes

- Route 1: Johnson-coherent low-signature witness selection.
- Route 2: joint aligned-payoff/low-row density.
- Route 3: matched parent Hellinger transport.

The main agent independently takes Route 4, with Route 5 as fallback.  Default
agent deadline is twenty minutes, with checks near five-minute intervals and
extensions only for a concrete near-complete proof or falsifier.
