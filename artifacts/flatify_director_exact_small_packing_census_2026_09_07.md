# Exact small-order support-packing census

Date: 2026-09-07. Status: **certified finite computation**, not asymptotic
evidence of a universal packing constant.

For a full signing A, let kappa(A) be the larger of the positive and
negative ternary fractional support-packing optima defined in
`flatify_construct_2026_09_07_affine_mub_partial_packing.md`.

| n | M_n | min_A kappa(A) | Number of labeled switched minimizers |
|---|---:|---:|---:|
| 3 | 3 | 3 | 2 |
| 4 | 4 | 4 | 6 |
| 5 | 4 | 5 | 12 |
| 6 | 5 | 6 | 12 |
| 7 | 9 | 68/7 | 3240 |
| 8 | 10 | 32/3 | 4200 |

At every displayed order a cap minimizer also minimizes kappa. At n=7,
different cap minimizers have kappa values 68/7, 21/2, and 56/5; at n=8
the values are 32/3 and 12. Thus selecting the child genuinely matters,
even among exact cap minimizers. It does not remove the packing gap at
these orders. No claim about selectable growing-order minimizers follows.

## Independent finite verification

`computations/flatify_director_small_packing_census_2026_09_07.py` uses
floating LP only to propose certificates. Every reported optimum has
nonnegative rational primal and dual variables with exactly equal values;
every ternary constraint is checked over integers and fractions.

To avoid relying on completeness of the graph atlas, the independent C++
program `computations/flatify_director_labeled_cap_census_2026_09_07.cpp`
enumerates every first-row-positive labeled signing, by integer Gray-code
updates of every spin energy. Switching makes this a complete class.
The separate checker
`computations/flatify_director_verify_labeled_census_2026_09_07.py`
then verifies that explicit residual vertex permutations of the stored
minimizing representatives cover **exactly** the complete C++ minimizer
list. It independently rechecks their rational packing certificates.

For every signing kappa(A)>=Q(A), by putting all packing mass on a full
spin vector. The next larger exact cap in the exhaustive census exceeds
or equals the displayed packing minimum, so no non-cap-minimizer can
improve that value. This last step avoids an assumption about atlas
completeness for the nonminimizing representatives too.

Results and source hashes are saved in
`computations/results/flatify_director_labeled_cap_census_2026_09_07.json`.
The generated executable is a reproducible build product; its C++ source
and exact checking program are the research-bearing artifacts.
