# Strategic steering

Evidence cutoff: blank-slate direct campaign, ledger Section 10.145
(2026-08-21), started from commit `b5ec773`.

Status: **blank-slate campaign complete with a strike; no architecture is
authorized to continue automatically**.

## User-stated objective and workflow directives

Determine whether `M_n/n^(3/2)` converges.  Convergence to any constant is
success; `1/2` is conjectural but is not the objective.  A rigorous proof of
nonconvergence is also success.

The user explicitly authorized one four-to-six-hour blank-slate campaign on
the original problem.  It had to derive three to five genuinely different
architectures from first principles, use archived no-go results only after
freezing them, attack the best route with parallel falsification, and avoid
the frozen composition, pressure, restriction, transport, posterior-state,
and sparse-repair branches unless independently revived.  Success required a
new asymptotic statement, a strict simpler reduction, or a genuinely new
architecture surviving serious falsification.  The README verification,
Git, consolidation, and blank-slate rules remain in force.  These are workflow
directives only; all mathematical judgments below are agent-authored.

## Agent-authored rigorous frontier

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.
```

The interval is unchanged.  The campaign proved new radial-moment,
rectangular, and covering-multiplicity theorems, but none supplies a strict
reduction of convergence.

## Blank-slate abstraction audit

Four candidates were frozen before archive comparison.

1. **Radial moment geometry.**  Minimize the smallest central support cap
   matching the first `d` even Krawtchouk coordinates of the augmented
   cut-code coset.  A squared Chebyshev polynomial proves, uniformly in `A`,

   ```math
   Q(A)\operatorname{sech}\!\left(
   {\operatorname{arcosh}(2^{(n-1)/2})\over d}
   \right)\le L_d(A)\le Q(A).
   ```

   At `d=alpha n` the relative error is `O(alpha^(-2))`.  This is an exact
   finite LP interface, but the required large-`alpha` oscillation theorem is
   equivalent by this sandwich to convergence itself.  The state contains
   the archived linear-degree signed-Eulerian cancellation problem.

2. **Microcanonical cut-code coverage.**  Work directly with the multiplicity
   `Z_r(b)` of augmented cuts in a Hamming ball.  The exact partial
   inclusion--exclusion identity shows that below the elementary sphere
   constant every useful odd Bonferroni truncation needs exponentially many
   terms; an even coverage certificate needs `2^(n-O(sqrt(n)))` terms.  Its
   affine center rank is then linear, exposing essentially the full coset
   profile.  A pair saddle changes stability at normalized score `1/2`, but
   does not decide whether any hole exists.

3. **Rectangular Banach projection.**  The exact fixed-width minimum is
   `m E|sum_{j<=k} epsilon_j|+O(k2^k)`.  A new degree-two Fourier correction
   detects squared column correlations.  Eigenvector product rounding proves
   that this correction is `O(n^(-1/3))` on every signing with
   `Q(A)=O(n^(3/2))`; it therefore vanishes on the minimization-relevant
   branch and leaves a constant below the known frontier.

4. **Exact-minimizer stationarity.**  Edge flips force light/heavy cut
   witnesses, but a scalable signing with cap `n^2/4` is one-edge locally
   minimal.  An order-eight cap-12 signing is stable under every one- and
   two-edge flip although the exact optimum is 10.  Bounded local
   stationarity cannot characterize minimizers.

The archive comparison found genuine new statements in all four calculations,
but it maps their unresolved obligations to signed high-moment cancellation,
full coset multiplicity, flat-spectrum Boolean resonance, or global
optimization.  No candidate passes the campaign's strict-reduction test.

## Strongest surviving architecture and exact target

The only route not closed by its first quantitative test is **direct
augmented-cut-code coverage**, provided it is nonperturbative rather than a
raw finite-replica expansion.  No strictly smaller convergence lemma was
identified.  The exact remaining statement suggested by the pair transition
is:

> For every fixed `epsilon>0`, prove
> `F_{n,1/2-epsilon}(1)=Pr_b{Z_r(b)=0}=0` for every sufficiently large `n` by
> a nonperturbative covering/isoperimetric argument.

This would attack the absolute quadratic optimum directly and is independent
of the frozen bridge/pressure machinery.  With the known upper bound it is
equivalent to convergence to `1/2`, so it is a conditional research direction,
not a class-A route.  The smaller first milestone is uniform coverage for one
fixed `c>0.336493364431...`, which would improve the lower frontier.

## Strongest current obstructions

1. Linear-degree Krawtchouk moments approximate the maximum uniformly but
   retain exponential signed-Eulerian cancellation.
2. Direct Bonferroni needs exponentially many centers; full factorial moments
   invert to the full coset-multiplicity histogram.
3. Pair/finite-replica overlap data controls typical multiplicity, not the
   existence of a worst uncovered word.
4. Rectangular scalar projections top out below the rigorous lower frontier,
   and their spectral covariance gain vanishes on competitive signings.
5. One-edge and radius-two stationarity permit macroscopically or finitely
   suboptimal signings.
6. The earlier action-recovery, finite-temperature, bridge, local-profile,
   and bounded-state routes remain frozen under their recorded no-go results.

## Ranked alternatives and restart conditions

1. **Worst-coset cut-code isoperimetry.**  Resume only with a theorem or
   primary-source mechanism that bypasses factorial-moment truncation and
   stores less than a linear-rank center profile.
2. **External direct variational theorem.**  A Γ-limit, coding-radius, or
   extremal-geometric theorem must verify exact signs and all orders and give
   a strict reduction, not another compactness restatement.
3. **Genuine nonconvergence.**  Still requires two infinite order sequences
   with fixed normalized separation; arithmetic examples alone do not count.
4. **Frozen architectures.**  Reopen one only after a new ingredient proves a
   quantitative implication that escapes its specific archived obstruction.

## Decision and next refresh

The campaign is a **STRIKE**: it produced rigorous local mathematics and an
exponential obstruction to raw Bonferroni truncation, but no improved bound,
strict reduction, or surviving class-A architecture.  Do not continue
automatically.  Refresh steering
immediately after a new user-authorized campaign or a decisive imported
theorem/counterexample.  If ordinary numbered waves ever resume, Wave 61
remains the next scheduled blank-slate boundary.
