# Finite absolute-overlap audit of exact and near-minimizer shells

**Status:** finite, reproducible computation; **not an asymptotic claim** and
not a canonical theorem.

The observables and populations were frozen in
[`nearmin_absolute_overlap_observable_freeze.md`](../experiments/nearmin_absolute_overlap_observable_freeze.md)
before the outcomes were inspected.  The executable audit is
[`nearmin_absolute_overlap_audit.py`](../experiments/nearmin_absolute_overlap_audit.py),
and the machine-readable output is
[`nearmin_absolute_overlap_summary.json`](../experiments/nearmin_absolute_overlap_summary.json).
The larger per-witness record is
[`nearmin_absolute_overlap_results.json`](../experiments/nearmin_absolute_overlap_results.json).

## Bottom line

The fixed-scale *one or two orientation-antipodal caps* obstruction does not
occur in any available exact-minimizer active shell at orders 7--14.

For every one of the 42 audited certified exact witnesses, the exact active
shell (`d=0`) contains at least three points pairwise satisfying

```math
 { |\langle z,z'\rangle|\over \binom n2}\le {1\over2}.
```

All these `d=0` packing numbers are exact certificates.  At the even orders
8, 10, 12, and 14, **every pair** of active-shell points satisfies the
inequality, so the exact packing is the entire shell, of respective sizes
8, 40, 20, and 156.  Because the observable takes an absolute value, changing
an oriented word `z` to `-z` does not hide this separation: the test already
works in projective edge space.

The evidence is encouraging for the premise of the conditional physical
compiler in `nearmin_absolute_overlap_physical_compiler.md`, but it is not a
proof of a uniform or growing packing.  Cardinality-matched controls also show
that the geometry is not uniformly more diffuse than a generic projective
subset: the comparison alternates sharply with the order and is often reversed
at odd orders.

## 1. Exact identities and registered statistic

For projective vertex spins `x,x'` and oriented edge words
`z=sigma c(x)`, `z'=sigma' c(x')`, put

```math
 R_E={|\langle z,z'\rangle|\over E},
 \qquad R_V={|\langle x,x'\rangle|\over n},
 \qquad E={n\choose2}.
```

The audit checks the exact identity

```math
 R_E=\left|{nR_V^2-1\over n-1}\right|.
```

Its packing statistic is the largest certified subset, or an explicitly
labelled deterministic lower bound, with every pair satisfying `R_E<=1/2`.
The shell deficits are the distinct members of

```text
0, 2, 4, 2 floor(sqrt(n)).
```

At these orders this is merely a finite thin-shell probe.

### Authoritative inputs

The immediate matrix source is the previously frozen and machine-readable
[`nearmin_blind_structural_results.json`](../experiments/nearmin_blind_structural_results.json).
It records the original source file and key path for every repository matrix;
those paths are copied into every per-witness record of the detailed output.
The exact values used for classification are

```text
(M_7,...,M_14)=(9,10,12,13,17,18,20,21).
```

Their machine certificates are the order-specific `exact_m*.json` files
through order 10, `computations/results/certified_m11_m12.json`, and
`computations/results/certified_m13_m14.json`.  The complete orbit sources at
orders 7 and 8 are respectively
`computations/results/m7_minimizer_orbits.json` and
`computations/results/m8_minimizer_orbits.json`.  Beyond order 8 no complete
orbit claim is used.

## 2. Certified exact-minimizer results

The table reports

```text
minimum R_E ; packing size
```

over the available exact witnesses at each order.  A range means different
available witnesses gave different values.  `LB` means a certified explicit
packing lower bound rather than a maximum-packing certificate.  The last
column is the widest preregistered `O(sqrt(n))` shell.

| `n` | witness count | active shell `d=0` | shell `d=2` | shell `d=2 floor(sqrt(n))` |
|---:|---:|---:|---:|---:|
| 7 | 3 | `.0476--.1429 ; 3--7` | `.0476 ; 13--15` | `.0476 ; 22--30` |
| 8 | 2 | `.1429 ; 8` | `.07143 ; 28` | `.07143 ; 52` |
| 9 | 6 | `0 ; 11--14` | `0 ; 11--14` | `0 ; 40--49` |
| 10 | 9 | `.06667 ; 40` | `.06667 ; 80` | `.06667 ; 80` |
| 11 | 9 | `.01818 ; 5--15` | `.01818 ; 26--34` | `.01818 ; 96--106 LB` |
| 12 | 8 | `.03030 ; 20` | `.03030 ; 72` | `.03030 ; 272 LB` |
| 13 | 2 | `.02564 ; 36` | `.02564 ; 36` | `.02564 ; 234 LB` |
| 14 | 3 | `.010989 ; 156` | `.010989 ; 364 LB` | `.010989 ; 364 LB` |

The counts at orders 7 and 8 are the complete authoritative orbit lists under
switching, permutation, and global sign.  At orders 9--14 they are all
available byte-distinct repository witnesses whose cap was independently
enumerated to equal the certified `M_n`; they are **not** an orbit-uniform or
complete sample.  Repeated values in the table may therefore reflect
equivalent representatives and are not statistical multiplicities.

Two finite patterns are especially clear.

1. Every exact active shell has an exact three-point projective certificate.
   Thus the most severe fixed-radius two-cap model is already false before any
   deficit is allowed.
2. The packing fraction has a strong parity/residue effect.  At even orders
   8--14 it is `1` in the active shell.  Among the odd-order witnesses it can
   be much smaller: the order-13 active shell has size 78 but exact packing
   36.  A claim of generic equidistribution would therefore overstate the
   evidence.

The minimum `R_E` generally attains the smallest or nearly smallest value
allowed by the finite overlap lattice.  This is much stronger finite evidence
against two near-antipodal projective caps than the signed first moment alone.

## 3. One-step near-minimizers

The audit contains 69 byte-distinct matrices with exactly recomputed cap
`M_n+2`: 9 repository witnesses, 18 independent greedy-search witnesses, and
42 cap-constrained-walk witnesses.  The discovery procedures for the latter
two strata are heuristic, but the cap classification is exact.

The exact active shell can be spiky: 9 of the 69 witnesses have certified
packing at most two (five order-8 shells of size two, and four order-11 shells
of size one or three).  This is a useful warning against claiming projective
packing merely from being one parity step above optimal.

After admitting deficit `d=2`, however, every one of the 69 shells has a
certified explicit packing of at least three.  The minimum lower bounds by
order are

```text
n=8: 10,  n=9: 4,  n=10: 12,  n=11: 10,
n=12: 24, n=13: 16, n=14: 144.
```

This is finite evidence that the vanishing-width version of the absolute-
overlap premise may be more stable than its exact-active version.  It does not
establish a uniform theorem.

## 4. Random and cardinality-matched controls

The program also enumerates 24 independent uniform random signings at each
order.  Their active shells are usually much smaller: only 36 of 192 random
active shells have a packing certificate of at least three, versus all 42
exact-minimizer active shells.  At deficit two the random count rises to 104
of 192.  This comparison is **not cap-matched** and is confounded by shell
cardinality.

The stronger control samples 32 uniform projective subsets with the same
`(n, shell size)` as every physical shell.  On the exact active shells:

- 12 of 42 physical packing numbers exceed the maximum of their 32 matched
  greedy controls;
- 14 of 42 fall below the minimum matched control;
- the remainder overlap the matched range.

The order-level behavior is systematic but nonmonotone.  For example:

- order 10: physical `40`, matched range `28--34`;
- order 11, shell size 24: physical `15`, matched range `18--24`;
- order 13: physical `36`, matched range `65--74`;
- order 14: physical `156`, matched range `131--146`.

Thus exact minimization creates broad active shells at these orders, but does
not select one universal kind of within-shell projective geometry.  In
particular, the observed odd/even contrast should be treated as a finite
structural clue, not a trend toward an asymptotic conclusion.

## 5. Orientation-antipodal and signed-balance stress test

The absolute-overlap packing deliberately quotients the global orientation.
Consequently its three-point certificates cannot be artifacts of placing one
cluster near `z_0` and another near `-z_0`.  Some order-11 exact witnesses make
this especially transparent: all 24 active words have the same energy
orientation, yet their exact projective packing is 15.

After the frozen run, an auxiliary signed-barycentre statistic was added to
stress-test the dichotomy in the conditional compiler.  It is explicitly
labelled post hoc in the JSON.  For an oriented shell law uniform on its
atoms, it records

```math
 {1\over E}\sum_e|\mathbb E Z_e|
 \quad\hbox{and}\quad
 {1\over E}\sum_e(\mathbb E Z_e)^2.
```

The median first quantity in exact active shells is

| `n` | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| median bias | .524 | .429 | .339 | .289 | .327 | .273 | .256 | .231 |

By contrast, every observed one-step-near active shell with packing at most
two has bias at least `.571`.  Therefore this finite sample contains no
example that is simultaneously marginally diffuse and projectively confined
to two caps.  This does **not** invalidate the exact synthetic two-cap
countermodel: it only says that the countermodel was not realized by these
small exact/near-minimizer shells.

## 6. Consequence for the conditional physical compiler

At the finite orders audited, every exact minimizer satisfies the compiler's
premise (AO.1) with

```math
gamma=1/2,
\qquad d=0,
\qquad |U|\ge3.
```

The one-step-near witnesses satisfy it with `d=2`.  This is the intended first
falsification test, and the premise survives it decisively.  It upgrades the
conditional theorem from an entirely unobserved hypothesis to one realized
by every available certified exact witness.

What it does **not** provide is the needed asymptotic structural lemma.  The
finite packing sizes cannot be extrapolated; available representatives are
not complete beyond order 8; matched controls show large residue effects; and
the compiler's probabilistic `O(n^(5/4))` remainder is asymptotic rather than
an effective small-order certificate.  The exact open statement remains:

> Does every exact minimizer admit, in a shell of deficit `o(n^(3/2))`, a
> fixed-`gamma` projective packing whose size is uniformly at least three, or
> preferably tends to infinity?

The three-point version is already a strict finite structural target and
retains only witness energies and pairwise overlaps, not the full Boolean
landscape.  A growing version would feed directly into the conditional
physical contextual-packing theorem.

## 7. Reproduction and checks

Run from the repository root using the project virtual environment:

```bash
.venv/bin/python extremal_information/experiments/nearmin_absolute_overlap_audit.py
```

The saved checks pass:

- every selected cap independently reproduces and is classified as `M_n` or
  `M_n+2`;
- the edge/vertex overlap identity is checked numerically against its exact
  integer numerator;
- all registered shell statistics are invariant under random diagonal
  switching, vertex permutation, and global matrix sign;
- 345 of 422 physical shell packings are exact certificates; the rest retain
  only certified lower bounds and the trivial upper bound.

No result in this note is promoted beyond finite experimental evidence.
