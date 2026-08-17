# Near-minimizer absolute-overlap audit: frozen protocol

Status: **FROZEN BEFORE INSPECTING THE OUTCOMES OF THIS AUDIT**

Freeze date: 2026-08-17 UTC

This audit is deliberately narrow.  It asks whether the positive near-top
shells of available exact and one-parity-step near-minimizers look, at the
small certified orders, like one or two almost-projectively-parallel clusters.
It does not test a selector, a continuation, or an asymptotic statement.

## Objects and exact identities

For a hollow signing `A` of order `n`, write

```math
h_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_x|h_A(x)|,
```

and enumerate projective spins with `x_0=1`.  At additive deficit `d`, the
positive oriented shell is

```math
\mathcal S_d(A)=
\{z=\sigma c(x):\ \sigma h_A(x)\ge Q(A)-d\},
\qquad c(x)_{ij}=x_ix_j.
```

All registered deficits below keep `Q(A)-d>0`, so every projective `x` has at
most one admitted orientation.  Put `E=binom(n,2)`.  For two distinct shell
points, the frozen projective observables are

```math
R_E(z,z')={|\langle z,z'\rangle|\over E},
\qquad
R_V(x,x')={|\langle x,x'\rangle|\over n}.
```

They obey the exact identity

```math
R_E=
\left|{nR_V^2-1\over n-1}\right|.
```

Thus both are retained: the edge quantity is the proposed missing observable,
while the vertex quantity shows its projective origin and avoids a misleading
monotonicity claim near `R_V=0`.

## Frozen inputs and provenance strata

Orders `7,...,14` are audited.  Matrices are read first from the existing
machine-readable blind audit and are deduplicated by their actual upper
triangle.

1. **Certified exact minimizers.**  Through order 8 use every authoritative
   exhaustive switching/permutation/global-sign orbit representative.  At
   orders 9--14 use every available byte-distinct repository representative
   whose exhaustively recomputed cap equals the supplied exact `M_n`.  This is
   not claimed to be an orbit-uniform sample.
2. **Certified-cap one-step near-minimizers.**  Use available repository
   representatives with cap `M_n+2`.  Supplement them, in separately labelled
   strata, by deterministic hash-selected low-cap witnesses from the blind
   audit's greedy and cap-constrained searches.  Search provenance is
   heuristic; the saved matrix cap and its relation to the certified `M_n`
   are recomputed exactly.
3. **Random-signing controls.**  Generate 24 independent uniform upper
   triangles per order from fixed seed `20260817` and recompute every cap and
   shell exactly.  These controls are not cap-matched.
4. **Cardinality-matched geometry controls.**  For every audited physical
   shell, draw 32 fixed-seed uniform subsets of the projective cube having the
   same cardinality.  These controls distinguish geometry of the near-top
   shell from geometry caused merely by its size.

No finite population average is treated as orbit-weighted unless the source
is the authoritative exhaustive orbit list.

## Frozen deficits and statistics

For each matrix use the distinct values in

```text
d in {0, 2, 4, 2 floor(sqrt(n))}.
```

These are exact lattice deficits and are all `O(sqrt(n))`; at the audited
orders they are only finite thin-shell probes, not an asymptotic regime.

For every shell of size at least two, compute exactly:

- shell size and the two orientation counts;
- the minimum, quartiles, median, and maximum of `R_E` and `R_V` over all
  unordered distinct pairs;
- the fractions with `R_E>1/2`, `R_E>=3/4`, and `R_E>=9/10`;
- a deterministic greedy packing at the fixed separation

  ```math
  R_E(z,z')\le {1\over2}.
  ```

  When an exact maximum packing certificate is computationally inexpensive it
  is also recorded; otherwise the greedy number is explicitly only a lower
  bound.

For singleton shells the pair statistics are null and the packing is one.
For matched-subset controls, record empirical quantiles of the minimum overlap
and the same deterministic packing lower bound.

## Frozen interpretation and falsifiers

The phrase *two nearly antipodal clusters* is not assigned an unearned formal
model.  The following operational diagnostics are frozen instead.

- `packing_(1/2)>=3` is a finite certificate that the shell cannot be covered
  by two subsets each having all cross-projective edge lines mutually more
  than `1/2` aligned under the particular greedy witnesses found.  If the
  packing optimum is certified, the corresponding statement is exact.
- `packing_(1/2)<=2` **only when certified exact**, together with a large
  minimum `R_E`, is finite evidence compatible with a one/two-cap geometry.
- A shell is called *less diffuse than its cardinality-matched null* only via
  its registered empirical percentile; no visual judgment is used.

The proposed absolute-overlap direction is decisively unattractive at these
orders if generic random-signing shells and certified exact shells have the
same or stronger diffuse certificates throughout, or if the exact shells
remain one/two-cap-compatible even after the `O(sqrt(n))` deficit is admitted.
It receives finite encouragement only if exact/near-minimizer shells acquire
several fixed-separated projective directions more reliably than both random
signings and cardinality-matched subsets.  Neither outcome is an asymptotic
claim, and no trend over eight small orders is extrapolated.
