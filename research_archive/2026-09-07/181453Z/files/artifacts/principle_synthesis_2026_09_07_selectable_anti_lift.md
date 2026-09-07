# Selectable exact-child anti-lifts: finite audit and the actual convergence obligation

Date: 2026-09-07. Status: finite solver-certified results at orders 5--7;
bounded inconclusive diagnostics at order 8. No asymptotic transfer theorem.

## Exact scope and normalization

For a symmetric hollow full signing `A` of order `n` and a skew hollow full
sign matrix `C`, put

```
L(A,C) = [[A,C],[-C,-A]].
```

There are exactly `n` zero matching edges. Filling these costs at most `n`
in cap, regardless of the chosen signs. The exact cut formula is

```
Q(L) = 2 max_{z,w disjoint Boolean supports}
               (|z^T A w|+|z^T C w|) = 4 D(A,C).
```

The supports partition `[n]`, and both directions are included. In
particular a signed cancellation on one cut cannot prove the needed joint
magnitude bound: reversing the cut reverses the `C` response but preserves
the `A` response.

The present optimization is

```
D_n^sel = min {D(A,C): Q(A)=M_n, C skew full sign}.
```

It chooses **both** the exact child optimizer and its bridge. It does not
fix one width minimizer, demand success for all exact minimizers, or impose
a preferred bridge law. This is the relevant finite selectable quantifier.

## Computation and independently checked witnesses

The code is `computations/principle_synthesis_2026_09_07_selectable_anti_lift.py`.
It imposes all projective Boolean child constraints with the independently
archived certified value `M_n`, and all ordered disjoint-cut tests for the
lift. Both sets of signs are variables. Simultaneous switching of both
copies gauges the first row of `A` positive. Rooted graph complementation
and the usual remaining vertex permutations give valid symmetry reductions;
the unrestricted skew variables are transformed along with the child.
Returned matrices are separately tested on every projective Boolean state
of the full `2n`-vertex partial lift, confirming `max H=-min H=4D`.

| n | M_n | selectable D | solver conclusion | M_n/sqrt(2) |
|---|---:|---:|---|---:|
| 5 | 4 | 3 | OPTIMAL | 2.828427... |
| 6 | 5 | 5 | OPTIMAL | 3.535533... |
| 7 | 9 | 6 | OPTIMAL | 6.363961... |
| 8 | 10 | at most 8, at least 5 | FEASIBLE, not closed | 7.071067... |

The order-7 proof closed after the symmetry reduction in about 22 seconds.
The 45-second joint order-8 run remained inconclusive. A separate 90-second
order-8 joint run gave the same feasible value 8 but a weaker lower bound 4.
No missing solver proof was inferred from repeated agreement of incumbents.

As a bounded follow-up, the installed seven-vertex graph atlas has four
rooted graph-isomorphism representatives with original cap 10, paired by
complementation. The two representative types have spectra
`{+-sqrt(13), +-sqrt(5) each three times}` and
`{+-3 each three times, +-1}`. Fixed-child searches for atlas indices 580
and 722 each returned `6 <= D <= 8` after 90 seconds. These do **not**
settle whether the selectable target `D <= 7` is possible at order 8.
This atlas observation is a finite diagnostic, not an imported structural
classification theorem about asymptotic minimizers.

Full primary-run matrices and solver outputs are saved in
`computations/results/principle_synthesis_2026_09_07_selectable_anti_lift.json`;
the two atlas runs are saved separately in
`computations/results/principle_synthesis_2026_09_07_selectable_anti_lift_atlas.json`.
OR-Tools gives solver-certified computation, not a standalone formal proof
object. The returned upper-bound witnesses are independently exhaustive.

Thus the finite failures at `n=5,6` persist after allowing selection among
**all exact cap minimizers**. The earlier fixed-width-child audit did not by
itself establish that quantifier. Conversely, the exact selectable `n=7`
result satisfies the unfilled-lift target. None of these bounded-order
facts is a scalable obstruction: an `o(n^(3/2))` error can absorb all of them.

## What would actually imply an original limit

A sufficient doubling estimate would be, for some `eta>0`,

```
D_n^sel <= M_n/sqrt(2) + O(n^(3/2-eta)).
```

After filling the matching this gives

```
M_(2n)/(2n)^(3/2) <= M_n/n^(3/2) + O(n^(-min(eta,1/2))).
```

No compatibility of the selected minimizers across different orders is
needed for this value inequality: a fresh exact minimizer may be selected
at each order. Nevertheless, summable doubling loss alone does **not**
force a common limit across all orders. It controls dyadic descendants;
logarithmically periodic value profiles are not ruled out.

For example, the following genuinely sufficient two-multiplier statement
makes the extra obligation explicit. Let `a_n=M_n/n^(3/2)`. If there are
`K,eta>0` such that, for all sufficiently large `n`,

```
a_(2n) <= a_n + K n^(-eta),
a_(3n) <= a_n + K n^(-eta),
```

then `a_n` converges. Along any composition of the two multipliers starting
at `n_0`, the accumulated increase is at most
`K n_0^(-eta)/(1-2^(-eta))`. For any fixed `delta>0`, the semigroup
`{2^a 3^b: a,b>=0}` is eventually multiplicatively `(1+delta)`-dense.
Indeed irrationality of `log(3)/log(2)` supplies a finite `delta`-net of
residues modulo `log(2)`; for a sufficiently large target logarithm, the
needed nonnegative power of 2 can then be added. Choose a seed order `n_0`
with `a_(n_0)` close to the liminf. For all sufficiently large `N`, choose
`T=n_0 2^a3^b` in `[N,(1+delta)N]`. Principal cap monotonicity gives

```
a_N <= (T/N)^(3/2) a_T
    <= (1+delta)^(3/2)
       [a_(n_0)+K n_0^(-eta)/(1-2^(-eta))].
```

Let the seed order tend to infinity along the liminf subsequence and then
let `delta` tend to zero. This proves equality of limsup and liminf under
the stated hypotheses. An appropriately uniform unequal-child gluing
theorem could replace the second multiplier.

Neither required asymptotically lossless transfer estimate is established
here. The newly proved all-order anti-family upper below `.499` in
`principle_synthesis_2026_09_07_anti_invariant_weave.md` removes the old
putative Boolean half-floor for the family, but does not impose
`Q(A)=M_n` on its selected first-half child and therefore does not fill
this remaining seed-landing gap.
