# Wave 12: same-window compatibility frontier

The verified weighted path theorem controls endpoint allocations at arbitrary
scales, but the temporal replenishment gaps in the deletion recurrence have not
been coupled to those allocations in the same order window.  The hybrid
replacement telescope is exact, but arbitrary optimal replacements can alter a
successor layer by order one even at order six.  Exact boundary recursion has
exponential state rank.

## Ten ranked routes

1. **Stopped-window temporal charging.**  Stop every endpoint chain when its
   order leaves `[rho n,n]`, keep only obligations born in that window, and seek
   a deterministic or averaged charging map from the centered gaps in (10.499)
   to the weighted buckets in (10.491), with no reuse outside the window.
2. **Dual certificate for the missing coupling.**  Formulate the strongest
   consequence of conservation, path-cover domination, global minimality, and
   the deletion identities as a finite laminar LP.  Dualize it: either exhibit
   a universal potential proving the needed coefficient or isolate a concrete
   realizable obstruction/extra hypothesis.
3. **Boundary-conditioned block minimizers.**  Replace a sibling block not by
   an arbitrary minimizer, but by a minimizer chosen lexicographically to fit
   the fixed cross mosaic and the already chosen tail.  Seek a Bellman choice
   for which the stability errors in (10.506) telescope or are one-sided.
4. **Orbit-averaged compatible replacement.**  Average all switchings and
   permutations of a fixed order-s minimizer against the fixed cross block.
   Use symmetrization or Khintchine bounds to prove that some orbit element has
   small positive deficit mismatch, then sum only across geometrically sized
   siblings.
5. **One-competitor randomized hybrid.**  Randomly switch every replacement
   block inside the single hybrid competitor of (10.504), keep original cross
   blocks fixed, and use global minimality only after averaging the complete
   variational expression rather than bounding layers separately.
6. **Mean terminal-excess jackknife.**  Double-count pairs of uniform
   deletions and use the overlap of their principal submatrices to control
   `E epsilon_H=E[Q(A[T])-q_m]` for at least one macroscopic child order.
7. **Adaptive multiplicative descent.**  Weaken the all-pairs tail condition
   (10.500): prove that from every large n one can choose one
   `m in [n/2,n)` with small normalized error, while ensuring the chosen orders
   cannot stall.  Determine the exact density/summability condition sufficient
   to bridge a liminf subsequence.
8. **Laplace representation of the reciprocal boundary partition.**  Insert
   `S^{-1}=integral_0^infinity exp(-tS) dt`, truncate its relevant t-range, and
   test whether the door-state sum regains a signed interpolation or Hölder
   inequality after this lift.
9. **Energy-resolution boundary compression.**  Quantize shore deficits and
   cross responses only to `o(n^{3/2})` accuracy.  Bound the number of relevant
   states for globally minimizing signings, rather than the full decoding
   table whose exact rank is `2^(n-2)`.
10. **Second-order graph-limit compactness.**  Search for a compact profile of
    minimizing signings at the `n^{3/2}` fluctuation scale (spectral measure,
    cut-field process, or exchangeable Gaussian limit) on which normalized Q
    is continuous and vertex deletion is sampling-consistent.

## Selected independent attempts

- Routes 1/2: derive or falsify a same-window temporal-to-endpoint charging
  lemma by an explicit primal/dual formulation.
- Routes 3/4: construct a compatible replacement rule, starting with exact
  orbit averaging and finite counterexample searches before claiming a bound.
- Routes 6/7: derive a macroscopic mean terminal-excess or adaptive descent
  lemma using two-level deletion identities and exact finite tests.
