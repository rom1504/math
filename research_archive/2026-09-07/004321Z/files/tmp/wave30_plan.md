# Wave 30 plan

Started: 2026-07-30T00:08Z.

Read before ranking: `README.md`, refreshed `STEERING.md` (Wave 29 evidence
cutoff), and the Wave 29 `Updated frontier` in `ledger.md`.

## Ranked ideas

1. **Adversarial-law planted-coset overlap.**  Work directly with the dual
   in (10.873).  For arbitrary selector weights `w`, jointly choose a
   selector, a favorable completion, and a planted balanced hash; seek a
   collision/overlap theorem giving one row-good coset `w(H_a)>=exp(-O(L))`.
   The proof must improve on the `1/binom(n,m)` self-atom and charge every
   selector KL explicitly.

2. **Nonlinear soft fractional cover.**  Replace the hard hit indicator by
   an exponential or positive-part payoff over a row-good coset, prove a
   minimax lower bound, and recover hard mass by a one-sided moment argument.
   Optimize the threshold before linearizing so the uniform-selector mean
   wall (10.824) is not silently reintroduced.

3. **Higher-order Johnson agreement/list recovery.**  Treat low-row
   exact/near-child completions as lists on selectors and seek a robust
   agreement theorem producing few global coordinate profiles.  Use
   triangles or higher faces, not only edges, and include refinement row
   cost so the exact A9 cap-80 wall is a checked hypothesis rather than a
   contradiction.  Search primary literature if a theorem appears to fit.

4. **Random refinement of signature classes.**  Given a globally
   low-signature witness family, randomly split its signature classes into
   balanced blocks and control the full coset row maximum with a matrix
   concentration argument centered at the minimal signature span.  Determine
   the extra field quantity that exactly replaces the false automatic
   refinement step.

5. **Bad-dual-law perturbation of the minimizer.**  Assume a selector law
   `w` defeats every row-good coset.  Use a nonlinear hinge/Laplace envelope
   to construct a signing perturbation or full cut contradicting exact
   minimality.  A purely linear `H_w` separation is disallowed by (10.823).

6. **Multiplicity of planted hashes.**  Quantify how many balanced hashes
   plant each aligned low-row cut while keeping the whole coset row-good.
   Double-count selector--coset incidences to seek expansion beyond
   self-collision; falsify the route if incidence can concentrate on
   essentially disjoint candidate neighborhoods.

7. **Matched parent common-mode control.**  At
   `beta=gamma=Theta(n^(-1/2+c))`, bound the parent likelihood entropy using
   exact-minimizer constraints on the centered row square exposed by
   (10.886)--(10.887), while keeping the adjacent Hellinger term separate.
   Seek a uniform theorem, not the non-uniform Taylor series itself.

8. **Two-geometry endpoint decomposition.**  Decompose the endpoint
   likelihood into its selector-common component and Johnson-varying
   component.  Control the former by parent coordinate flips/minimality and
   the latter by (10.798), with constants audited against the A4 and J-I
   walls.

9. **Coset-averaged conditional row tilt.**  Test whether averaging over
   planted partitions enforces the integrated covariance bound
   (10.879)--(10.880) even though individual spins fail it.  This remains
   below ideas 1--8 because the scalar reverse-tail input is also missing.

10. **Signed-cycle deficit-profile stability.**  Revisit constant shortfall
    only through the exact child profile (10.847) or signed deletion ratio
    (10.780), looking for a minimizer-specific complex-analytic or
    log-concavity theorem.  Do not use unsigned coefficients or plaquette
    aggregation already known to yield only constant reward.

## Selected agent routes

- Route 1: adversarial selector-law planted-coset overlap and exact dual.
- Route 3: higher-order Johnson agreement/list recovery, including a primary
  literature search and a precise hypothesis match or obstruction.
- Route 7: matched parent common-mode control at the intended temperature.

The main agent independently takes Route 2, with Route 5 as fallback.  Default
agent deadline is twenty minutes, with checks near five-minute intervals and
extensions only for a concrete near-complete proof or falsifier.
