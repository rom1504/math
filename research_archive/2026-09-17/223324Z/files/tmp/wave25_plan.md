# Wave 25 ranked plan

Evidence basis: `README.md`, §10.77.5 of `ledger.md`, and `STEERING.md` at the Wave 24 cutoff.

1. **Priced exact-ground separation.** For every row-square price `theta`, minimize an explicitly chosen penalized functional over the exact parent ground face and prove a lower bound for `u_g[1+theta(C-c_g)]`; seek an exact-minimizer exchange or normal-cone certificate that yields the two-ground LP target.
2. **Loss/row-square structural dichotomy.** Tilt exact grounds jointly by selector loss and `R_2`; prove that either a low-`R_2` ground has adequate lower-tail mass or a high-`R_2` ground forces enough compensating deletion gain to construct a second favorable LP column.
3. **Full endpoint exponential transport.** Attack precisely `Ent_{nu_beta}(C a)+E_{Ca nu_beta} KL(rho_d||U)` in (10.774)--(10.775), using the exponential Johnson-edge Dirichlet form and exact-minimizer competitor inequalities, not base variance or mean truncation.
4. **Adaptive random-permutation martingale.** Reveal deleted vertices along a random permutation while changing the parent ground only on certified zero-cost directions; try optional stopping or a potential supermartingale that lands at every target order with a soft captured-row budget.
5. **Two-parameter pressure/Fenchel alternative.** Study `log E exp(-alpha ell-gamma R_2)` and its convex dual; determine whether failure of the desired joint event rigorously implies a leading restriction gap, thereby matching or improving converse (10.769).
6. **Replenishment coarea inequality.** Use `ell=R_T+b_T-d_T` together with level sets of nonnegative row sums and exact block-replacement inequalities to prove a lower-tail estimate specifically for active densities `rho>=1/2`.
7. **Signed-cycle deletion cancellation.** Seek an identity for the average/logarithmic derivative of `P_{B[-i]}/P_B` that preserves signed cancellation and can yield a `sqrt(r)` cavity reward; avoid coefficientwise absolute bounds and pairwise plaquette aggregation.
8. **Exact-ground degeneracy graph.** Classify zero-deficit sign/block moves on the ground-state face and prove either enough connectedness to redistribute row-square or a rigid structure strong enough to estimate selector coverage directly.
9. **External-surplus service with soft costs.** Fold the exact Hall/replacement formulation into the new captured-mass LP, charging surplus only on served selectors rather than requiring a worst-case service schedule.
10. **Explicit algebraic-minimizer comparison windows.** Use conference/Paley-type exact or near-exact families to construct target-specific restriction edges on enough overlapping fixed-ratio windows; audit exact landing and summability rather than extrapolating from a subsequence.

Selected independent attacks: 1 (exact-ground/LP geometry), 2 (probabilistic loss--row dichotomy), and 3 (parent-Gibbs endpoint transport). Route 7 remains the root-side independent fallback because Wave 24 found its exact target but also a strong generic obstruction.
