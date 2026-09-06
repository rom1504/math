# Exact one-row insertion across ALL small minimizers

Status: exact exhaustive integer computation through order 9; independent
integer replay of one witness per insertion class. No asymptotic conclusion.

For `q_A=sum_(i<j) A_ij x_i x_j`, write

```math
E(A)=\min_{a\in\{\pm1\}^n}\max_x(|q_A(x)|+|a\cdot x|)-M(A).
```

All signings are enumerated after switching so that the first row off diagonal
is positive. There are `2^((n-1)(n-2)/2)` distinct gauged signings. Spins and
new rows are enumerated modulo global reversal. The C++ Gray-code scan uses
integer updates; every row is tested for every minimizing signing.

| n | M_n | Gauged minimizers | Histogram E: count | 3 M_n/(2n) |
|---|---:|---:|---|---:|
| 3 | 3 | 2 | 1: 2 | 1.5 |
| 4 | 4 | 6 | 0: 6 | 1.5 |
| 5 | 4 | 12 | 1: 12 | 1.2 |
| 6 | 5 | 12 | 4: 12 | 1.25 |
| 7 | 9 | 3240 | 1: 2520; 3: 720 | 1.928571... |
| 8 | 10 | 4200 | 2: 4200 | 1.875 |
| 9 | 12 | 1607760 | 1: 725760; 3: 882000 | 2 |

Thus choosing an arbitrary exact minimizer can fail to attain the next-order
minimum already at n=7 and n=9. A selectable-parent qualification is real,
not a technical afterthought. The uncorrected pointwise sharp insertion bound
also fails at n=6 and n=8. These are NOT counterexamples to an asymptotic
bound with an `O(n^(1/2-delta))` remainder.

At n=6 the best one-vertex extension raises the cap by 4, but at n=7 some
parents raise it only by 1. Finite insertion spikes therefore should not be
mistaken for normalized nonconvergence.

## Identical complete deficit spectrum, different insertion value

There are six absolute-energy histograms among all order-9 minimizers. One
histogram occurs with BOTH E=1 and E=3. The witnesses with edge codes 898008
and 6737136 have, on the 256 projective spins,

```text
absolute energy:   0    4    8   12
multiplicity:    60  111   60   25
```

The first extends optimally to cap 13, whereas the second cannot extend below
cap 15. Both have parent cap M_9=12. Therefore the entire absolute deficit
partition function, not merely the ground-state count, fails to determine
sharp insertion or select the good parent class. This is an exact finite
geometry obstruction, not a claim about asymptotic entropy criteria.

The independent replay includes both matrices, all row values, and assembled
optimal extensions. Their signed energy histograms differ, so no claim is
made that the signed (rather than absolute) histogram is also identical.

## Reproduction

```bash
g++ -O3 -std=c++17 computations/decisive_independent_all_small_extensions_2026_09_06.cpp -o computations/decisive_independent_all_small_extensions_2026_09_06
computations/decisive_independent_all_small_extensions_2026_09_06 9
.venv/bin/python computations/decisive_independent_extension_replay_2026_09_06.py
```

The independent NumPy replay reconstructs the gauged matrices from saved edge
bits, recomputes every energy and extension row using matrix multiplication,
then exhaustively checks the assembled parent. Its saved JSON includes matrices,
best rows, and full row-cap histograms:
`computations/decisive_independent_extension_replay_2026_09_06.json`.

One witness code per E class is: `(3,1):0`, `(4,0):1`, `(5,1):13`,
`(6,4):220`, `(7,1):828`, `(7,3):826`, `(8,2):53014`,
`(9,1):898008`, `(9,3):898023`. Free edges are lexicographic pairs
`1<=i<j<n`, and a set bit denotes coefficient -1.

## Live asymptotic target

The director's sufficient target is a selectable original-minimizer recurrence

```math
E(A_n)\le {3M_n\over2n}+O(n^{1/2-\delta}),\quad\delta>0.
```

Because `M_(n+1)<=M_n+E(A_n)`, this would bound the positive normalized drift
by a summable sequence (plus the elementary Taylor error), hence prove
convergence of the bounded normalized values. Selection separately at each n
is enough if these are exact minimizers: no nested parent chain is required.
For near-minimizers, their excess must also be at most the stated summable
recurrence error; merely `o(n^(3/2))` closeness is insufficient.

No scalable sub-half-cap example with `E(A)/sqrt(n)->infinity` has been found.
Simple constant block replication amplifies the normalized cap by sqrt(block
size) and therefore leaves the sub-half regime. Hadamard tensor lifts instead
introduce a vector-spin relaxation and do not preserve the seed cap constant.
