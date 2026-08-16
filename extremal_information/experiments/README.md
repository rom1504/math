# Exact finite laboratory

All programs in this directory are deterministic and use the repository-local
Python environment.  They are falsification tools; no finite pattern is
silently extrapolated.

Run from the repository root:

```bash
.venv/bin/python extremal_information/experiments/entropy_overlap_lab.py
.venv/bin/python extremal_information/experiments/pinned_query_rate_verify.py
.venv/bin/python extremal_information/experiments/build_quadratic_landscape_dataset.py
.venv/bin/python extremal_information/experiments/verify_code_replica_hierarchy.py
```

The scripts write only to this directory by default.

## `entropy_overlap_lab.py`

This verifies three independent claims exactly.

1. The codes
   `C={0000,0011,0101,0110}` and
   `D={0000,0011,0101,1001}` have the same ordered pair-distance enumerator
   and covering radii two and three.  Their Cartesian-power enumerators and
   additive radius separation are checked symbolically.
2. Two saved order-eight complete signings have identical energy histograms
   but caps `16` and `20` after the same all-negative one-vertex extension.
3. Every residual unlabeled graph in the NetworkX graph atlas is checked for
   orders four through eight.  No equal exact
   `(energy(x),energy(y),overlap(x,y))` signature has two different multisets
   of one-vertex extension caps in this range.

Output:
[`entropy_overlap_results.json`](entropy_overlap_results.json).

The third statement is a complete finite census only through order eight.

## `pinned_query_rate_verify.py`

This enumerates all `2^binom(n,2)` sign interactions through order five.  For
every landscape and every field direction it verifies the pinning identity,
recovers each edge coefficient by its degree-two Walsh coefficient, and
checks that all response vectors are distinct.

Output:
[`pinned_query_rate_results.json`](pinned_query_rate_results.json).

This script is a finite falsifier for the algebra in the information theorem;
the rate lower bound itself is proved analytically.

## `build_quadratic_landscape_dataset.py`

This builds the order-eight experimental dataset.  Switching makes the first
row positive; the remaining negative edges are represented by each unlabeled
graph on seven vertices from the NetworkX graph atlas.  The script groups all
`1044` residual graph representatives into `243` exact
energy--energy--overlap signature classes.

For one representative of every class it records:

- the exact matrix and energy histogram;
- positive, negative, and absolute extrema;
- exact near-cap counts;
- trace powers through degree eight;
- row sums; and
- the full histogram of one-vertex extension caps.

Output:
[`quadratic_landscape_order8.json`](quadratic_landscape_order8.json).

The dataset is intended for collision search and theorem discovery.  It is
not a machine-learning certificate and the chosen representative is not a
canonical switching/permutation normal form.

## `verify_code_replica_hierarchy.py`

This supplies exact finite checks for Theorem 3.3.  It verifies every proper
selected-codeword column profile for the parity constructions at
`r=3,5,7,9`, exhaustively computes covering radii at `r=3,5`, checks the
alternating-binomial identity, and enumerates the full ambient `t`-point
membership/distance census for the four-bit base pair.  The latter agrees
through `t=4` and first differs at `t=5`.

Output:
[`code_replica_hierarchy_results.json`](code_replica_hierarchy_results.json).
