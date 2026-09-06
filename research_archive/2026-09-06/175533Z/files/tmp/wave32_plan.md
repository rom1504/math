# Wave 32 idea ranking

Read before ranking: `README.md`, Wave 31 `Updated frontier`, and the Wave 31
version of `STEERING.md`.  Work at

```math
L_0=n^{3/4-c_0},\quad k_0=\Theta(L_0/\log n),\quad
r=\Theta(n^{1/4+c_0}\log n).
```

A batch cover by `T<=n^eta`, `eta<c_0`, eligible cosets suffices, with final
saving `c'=c_0-eta`.  Generic incidence, fixed-order agreement, KKT, and
uncentered harmonic accumulation are recorded walls.

## Ranked ideas

1. **Coordinatewise min-cut formula for completion forests (selected).** A
   full completion has free bits outside its selector.  For a fixed rooted
   tree of selectors, orientations, local favorable labels, and auxiliary
   root, minimize total Hamming edge length over all free bits.  The problem
   should decouple by coordinate into binary tree total-variation/min-cut
   problems.  Prove the exact formula, derive its dual/cut certificates, and
   identify the weakest local-label agreement statistic forcing total
   `O(k_0)` in fewer than `n^c0` groups.

2. **Low-row-center incidence power sum (selected).** For each low-row spin
   `x`, let `G_x` contain selectors admitting a favorable completion within a
   prescribed Hamming/tree radius of `x`.  Use the exact abundance of
   low-row centers, completion-cylinder freedom, Hamming-ball volume, and
   quadratic noise identities to lower-bound `sum_x w(G_x)^s` against every
   selector law—or produce a minimizer-scoped missing conditional moment.
   Normalize degrees by the total center count; raw abundance is insufficient.

3. **Centered harmonic screening at fixed density (selected).** Starting
   from (10.916)--(10.921), choose useful selector centers `alpha_S` and
   express the Doob increments or vertex-flip response through the
   renormalized screening residual.  Seek a real minimizer-specific variance
   bound, comparison theorem, or a rigorous wall showing which additional
   conditional susceptibility is necessary.  Do not revert to `E F_S`.

4. **Noise-to-low-row proximity around a favorable completion (main
   agent).** Compute the exact Hamming-sphere/noise expectation and variance
   of `R_2(x^F)` from (10.914).  Determine whether a favorable completion's
   sublinear ball must contain an `O(n^2)` center under exact minimality; if
   false from current data, construct a scoped PSD/complete-signing mechanism
   wall and isolate the needed local quadratic drift.

5. **Metric entropy of low-row cuts.** Use `E_U R_2=n(n-1)` plus the quadratic
   structure of `R_2` to estimate the covering/packing numbers of the low-row
   set at radius `k_0/r`.  Mere density one-half does not control worst-case
   fibers, so exploit Fourier degree two or signing constraints explicitly.

6. **Partial list-agreement with signature entropy.** Revisit primary
   agreement/list-recovery theorems now that only a `1/T` fraction is needed.
   Seek a weighted/adversarial-law theorem whose conclusion controls local
   label total variation or VC dimension, and verify every overlap premise;
   a generic citation without a derived premise is not progress.

7. **Discrete edge-flip susceptibility.** Apply exact-signing minimality to
   carefully correlated edge flips designed to expose
   `sum_i a_i(F)^2-r_i a_i(F)` in (10.914).  The goal is a true quadratic
   inequality, not another signed shore-total cycle identity.

8. **Adversarial-law structural dichotomy.** Separate laws concentrated on
   a few selector membership patterns from constant-distance code supports.
   Prove that any bad power-sum law must have a checkable irregular inclusion
   profile, then couple that profile to principal-ground constraints.  High
   min-entropy alone is explicitly insufficient.

9. **Finite forest-cover LPs.** Extend the exact `A_8,A_9` enumerations to
   compute minimum low-row-centered forest cover numbers for sampled/full
   slices and the associated partial power sums.  Use this only to discover
   identities or falsify proposed constants, not as asymptotic evidence.

10. **Multiscale forest clustering.** Allow coarse clusters with larger tree
    cost, refine only exceptional subclusters, and optimize the combined
    `T L_0` entropy and `D` row-cap costs.  Seek a hierarchical analogue of
    the one-level forest theorem that can spend part of `c_0` at each scale.

## Assignments

- Route 1: exact local-label/tree min-cut reduction and agreement criterion.
- Route 2: low-row-center incidence and power moments.
- Route 3: centered harmonic-screening variance.
- Main agent: exact Hamming-sphere drift of row square (Idea 4), with finite
  signing audits and a sharply scoped obstruction if the drift is wrong-way.
