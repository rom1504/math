# Wave 26 ranked plan

Evidence basis: `README.md`, §10.78.5 of `ledger.md`, and `STEERING.md` at the Wave 25 cutoff.

1. **Low-row-square reverse tail on the Johnson slice.** Choose a cut through a penalized `R_2`/selector-payoff principle and prove the fixed-cut event in (10.795) by a genuine lower-tail or level-set theorem for a degree-two slice polynomial.
2. **Noise interpolation from child or parent cuts.** Recompute the exact payoff, parent deficit, row square, and information under coordinate noise while retaining `B_{n,m}`; test whether a noise level now reaches (10.795), or prove the corrected capacity/cost obstruction.
3. **Penalized all-cut minimax.** Optimize `log U_m exp[-lambda(hat ell)_+] + eta R_2` over all cuts; use discrete signing minimality or an edge-flip local optimum to defeat every penalty without invoking the circular mean.
4. **Parent-cut entropy transport.** Express `Ent_{nu_beta}(Ca)` as a Gibbs likelihood-ratio entropy and seek a parent-cube log-Sobolev, competitor, or conditional-expectation bound at `O(n^(1/2-2c))`.
5. **Finite-measure Hellinger edge theorem.** Control both weight variation and lifted-channel Hellinger terms in (10.798) by a complete exact-minimizer comparison, independently of the parent-cut term.
6. **Signed-cycle deletion stability.** Use a multivariate cycle polynomial, vertex-support derivative, or cancellation-preserving transform to bound `P_{B[-i]}/P_B`; avoid absolute coefficients and linear plaquette aggregation.
7. **Low-cost completion of child grounds.** Given a ground on `S`, use vector discrepancy to choose outside signs with small full `R_2`, then quantify how many selectors share the resulting completions.
8. **Minimizer-cover/Delsarte correlation.** Combine the augmented-cut-code cover with the low-`R_2` subset of the cut cube to force one codeword with many centered-good selectors.
9. **External-surplus soft service.** Charge terminal Hall/replacement surplus only on selectors captured by one low-row-square cut, matching the new event rather than a worst-case schedule.
10. **Algebraic minimizer windows.** Test conference/Paley minimizers for an explicit low-`R_2` cut and fixed-density selector large deviation that could cover overlapping comparison windows with exact landing.

Selected independent attacks: 1 (direct leading lemma), 4 (parent-cut entropy half of the strongest alternative), and 6 (complete-signing local alternative). The main agent will pursue 2 and use it to decide whether noise smoothing should remain active.
