# Wave 14: ten ranked continuations after ledger 10.66

1. **Windowed sibling-certificate matching.**  Treat each positive residual
   bucket's witness `sigma u^T D u <= -r` as a signed resource and formulate a
   fractional hypergraph matching that prevents reuse of one sibling block or
   one incompatible sign across incomparable temporal histories.  Prove a
   window-level antichain bound or find a finite obstruction.

2. **Endpoint ground cross-Gram theorem.**  At a balanced parent, prove the
   exact identity `total residual(p,n)=|p^T A n|` for positive and negative
   grounds.  If every pair is neutral, their spans are `A`-orthogonal.  Combine
   this with rank/frame/Frobenius constraints to force a structural payment,
   and audit A8 and conference examples.

3. **Critical hybrid low-margin entropy via quadratic Littlewood--Offord.**
   Translate the low-slack cross profile in (10.549) into concentration of a
   Boolean quadratic form.  Use the robust-dependence theorem of Kwan--Sauermann
   (arXiv:2312.13826) or their algebraic inverse theorem (arXiv:1909.02089) to
   seek a dichotomy: subexponential near-top profile count, or low-rank/additive
   structure that can be handled separately.  Check whether the published
   hypotheses are actually strong enough; a polynomial point-mass bound alone
   is not enough.

4. **Chaining instead of a union bound for hybrid switching.**  Put the
   canonical metric on the internal randomized-switching process indexed by
   the low-margin cross states.  Bound its Gaussian/Rademacher complexity by
   overlap entropy rather than cardinality, with a target `o(n^(3/2))` at
   block size `sqrt(n)`.

5. **Pure-versus-mixed hybrid rounding.**  Compute the exact mixed minimax
   value for the A9 `3+6` wall and search for a general rounding theorem from a
   small-support mixed selector to one pure block tuple.  Quantify the
   integrality gap and relate it to the number of shared active states.

6. **Internal excess versus Hamming distance.**  Test and then prove or
   falsify a stability inequality relating `dist_H(D,M_m)` to
   `Q(D)-q_m` for induced blocks of global minimizers.  The unrestricted
   statement is likely false; search for a version averaged over an adaptive
   partition.

7. **Adaptive paving partition.**  Choose the critical blocks rather than
   fixing them, seeking simultaneously small cross-only overshoot and small
   total distance to local minimizers.  Audit whether rank-one matrix
   discrepancy/paving results such as Kyng--Luh--Song (arXiv:1901.06731)
   control the Boolean `Q` quantities or only operator norm.

8. **Amplify the A9 joint-selector wall.**  Attempt substitution, regular
   blow-up, or tensor-like compositions that preserve global minimality (or
   asymptotic minimality) while turning the fixed excess four into a critical
   `Theta(n^(3/2))` wall.  A successful construction would falsify the current
   hybrid route; a no-amplification theorem could become the needed positive
   input.

9. **Credit-bank potential and interval Hall cuts.**  Allow negative credits
   to be stored forward while residual resources may pay backward inside a
   contracted order window.  Derive the exact interval-cut dual on a chain and
   seek a scalar potential whose oscillation, not every prefix/suffix, is
   tail-summable.  Test it on the A7 and A5 timing walls.

10. **Endpoint-neutral temporal dichotomy.**  Characterize parents for which
    all optimized endpoint residuals vanish.  Prove that such a parent either
    has conference/orthogonality structure, forces a negative centered
    deletion increment in the same order window, or carries terminal excess.
    Exhaust small minimizers to identify the correct quantitative statement.

## Ranking and selected routes

- First: idea 1, because it directly targets the maximum-antichain error and
  uses the new exact negative sibling witness.
- Second: idea 2, because the balanced cross-Gram identity is exact and may
  turn the A8 wall into exploitable linear algebra.
- Third: idea 3, because it attacks the independent critical hybrid route and
  the literature provides a precise theorem to audit rather than a vague
  analogy.

Ideas 4, 5, 6, 8, and 9 are immediate successors if one selected route is
falsified.  Idea 7 is lower ranked because operator-norm paving may miss `Q`;
idea 10 is broader than the cross-Gram formulation in idea 2.
