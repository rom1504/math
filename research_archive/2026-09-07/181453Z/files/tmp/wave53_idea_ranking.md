# Wave 53 independent idea ranking

Evidence basis: README, STEERING at Wave 52, and ledger §10.105 Updated
frontier, reread after commit ec81d15. Scores are promise/testability/
independence on a 1--5 scale. None is a user directive.

1. **Bounded-coordinate K-functional amplification (5/5/5).** Derive a sharp
   deterministic lower bound for `K_{1,2}(beta,sqrt H)` from cross energy and
   coordinate concentration, then prove or falsify that exact-minimizer
   far-margin states satisfy the scalar profile on saved mass. This directly
   attacks (10.1277) without any independence assumption.
2. **Low-cross box witness through internal-field regularity (5/5/5).** Use the
   exact random-completion expectation to reduce the box row to cross energy
   plus `||A[S]y||_2^2`; seek a minimizer-specific averaging or selector choice
   that controls the latter at `O(n^(9/4-c))`. This is the complementary branch
   to (10.1278).
3. **Linear-core intersection entropy (5/4/5).** Rewrite large-core retention
   through the full intersection histogram `H_z(j)` and seek an entropy,
   compression, or forbidden-intersection inequality forcing saved positive
   excess at `-log rho_ell=omega(L)`. Explicitly distinguish a real global
   mechanism from a relabeled degree bound.
4. **Direct level-two pair-load incidence theorem (4/5/4).** Bypass spectral
   amplification and prove `P_2>=lambda_2(2)` at some controlled cap by double
   counting exact child-ground incidences. Wave 52 shows this alone supplies a
   polynomial center degree; the finite self-loop cannot be used.
5. **Cross-energy population dichotomy (4/4/4).** Average the exact identities
   for `g,V,D` over local states to force saved mass either below the box cross
   scale or above it with a usable K-profile. This must avoid the two-moment
   support wall by using structural conditioning, not another raw moment bound.
6. **Cancellation-stable large-radius exchange (3/4/4).** Iterate selector
   exchanges while charging actual row increments and cancellation rather than
   the one-flip triangle bound (10.1265). Prove a linear-radius low-row closure
   theorem or find an exact finite obstruction.
7. **Block-reveal completion tail (3/4/5).** Reveal outside spins in blocks,
   apply conditional Montgomery-Smith tails to fresh cross-linear pieces, and
   control the evolving quadratic part by a martingale cutoff. A successful
   theorem must retain `exp{-O(H)}` probability and avoid multiplying dependent
   marginal events.
8. **Weighted multiscale core transform (3/5/4).** Replace one fixed `ell` by a
   nonnegative generating function of `binom{|S cap T|}{ell}` and optimize its
   spectrum. Test whether a weighted kernel can obtain a saved useful mode while
   its direct-collision coefficient is `exp{-omega(L)}`.
9. **Outside-field truncation and exceptional-coordinate exchange (3/5/3).**
   Split K-profile failure into a regular outside field and a small set of huge
   coordinates; use maximal-selector exchange only on the exceptional set and
   the K-functional theorem on the regular part. This fairly uses, but does not
   overclaim, the Wave 52 port bounds.
10. **Projected vector balancing for the box row (2/4/5).** Project outside
    columns away from the span of large internal ground fields, sign the residual
    by a vector-balancing theorem, and pay the low-dimensional spike separately.
    It advances only if it beats the random-completion `O(n^2)` baseline plus the
    internal Gram term at the exact target scale.

Selected parallel attacks: 1, 2, and 3. The main agent will independently
derive scalar K-functional/profile inequalities and audit the agents' claimed
bridges; ideas 4 and 5 are immediate successors if a selected attack closes by
counterexample.
