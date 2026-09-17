# Exact full-column response games on finite minimizer nearcodes

2026-09-17. These finite experiments include EVERY physical sign column,
not a Gaussian surrogate, a sampled dictionary, or a radial subclass.
The stored matrices' global-minimum labels are imported from archive
provenance. Their caps, complete nearcodes, and the primal/dual response
certificates are independently enumerated here; we do not give a new
global-minimality proof or infer asymptotic optimizer structure.

## 1. Game and exact certificates

For a stored full signing A of order n and integer deficit window T, put
`C_T={x: Q(A)-|H_A(x)|<=T}`. Define

```
beta(C_T)=min_nu max_(x in C_T) E_nu |h dot x|,
beta_iso(C_T)=the same minimum with E_nu hh^T=I.
```

All projective columns h=(1,h_2,...,h_n) are included. A global fair
sign makes any such law centered without changing correlations or
absolute responses. Both polarities of the energy and every nearcode
word are retained (projective representatives suffice).

A rational primal probability law gives an upper bound by evaluating
all responses, and in the isotropic case all pair correlations.
A lower certificate consists of nonnegative code weights pi_x summing
to at most one and a polynomial

```
sum_x pi_x |h dot x| >= b + sum_(i<j) d_ij h_i h_j
                         for EVERY physical sign column h.       (1)
```

For the unrestricted game all d_ij are zero. Averaging (1) under an
admissible column law proves its maximum response is at least b.
Exact matching bounds certify an optimum. Where exact reconstruction
fails, the numerical dual is rounded rationally, its code weights
normalized downward if needed, and b lowered by the exact worst
pointwise violation. This gives a rigorous rational interval; numerical
endpoint agreement is never promoted to equality.

An unrestricted certificate also certifies the isotropic game whenever
its lower bound matches an isotropic primal. This closes the n12 and
n14 ground-code values despite failed direct reconstruction of their
isotropic solver duals.

## 2. Complete ground-code results

Values below are UNNORMALIZED mean absolute overlaps; divide by sqrt(n)
when comparing asymptotic response constants. The final column retains
both maximal radial endpoint covariance constraints before mixing to I.

| Stored case | Q(A) | Projective ground words | beta | beta_iso | Balanced maximal-radial subclass |
|---|---:|---:|---:|---:|---:|
| n4 | 4 | 2 | 0 | 1 | 4/3 |
| n6 | 5 | 12 | 5/3 | 5/3 | 5/3 |
| n8 class0 | 10 | 8 | 3/2 | 7/4 | 7/4 |
| n8 class1 | 10 | 8 | 3/2 | 3/2 | 5/3 |
| n10 | 13 | 40 | 53/25 | 144/65 | 747/335 |
| n12 | 18 | 20 | 9/5 | 9/5 | 9/4 |
| n14 | 21 | 156 | 98/39 | 98/39 | 98/39 |

Every displayed value in this table is exact. The n12 radial endpoint,
initially confined to a narrow rational interval, was subsequently
closed by exact active-face row reduction and all-column verification.

Three independent effects are visible on actual finite examples:

- Isotropy has a genuine cost at n4, n8 class0, and n10.
- Maximal radial endpoint covariance has a further cost at n4, n8
  class1, n10, and especially n12: its best response is exactly
  9/4 versus exactly 9/5 for an unrestricted isotropic law (25% more).
- The same fixed covariance does not fix the absolute response:
  n8 class0 has two exactly certified balanced radial laws with
  responses 35/16 and 7/4, as recorded in the radial artifact.

These observations do not imply a monotone pattern in n or a limiting
value. In particular the maximal radial construction is neither
universally optimal nor universally suboptimal for these finite games.

## 3. Shell deterioration

| Case | T | Projective code size | beta | beta_iso |
|---|---:|---:|---:|---:|
| n4 | 2 | 6 | 4/3 | 3/2 |
| n4 | 4 | 8 | 3/2 | 3/2 |
| n6 | 2,4 | 32 | 15/8 | 15/8 |
| n8 class0 | 2 | 28 | 20/11 | 40/19 |
| n8 class0 | 4 | 52 | 52/25 | 40/19 |
| n8 class1 | 2 | 28 | 33/17 | 2 |
| n8 class1 | 4 | 52 | 2 | 40/19 |
| n10 | 2,4 | 120 | 103/45 | 79/33 |
| n12 | 2 | 92 | 452/185 | 1323/536 |
| n12 | 4 | 224 | 15223/5963 | 7086/2761 |
| n14 | 2,4 | 520 | 4837/1768 | 4837/1768 |

All 42 unrestricted/isotropic ground and shell games now have exact
matching rational primal and dual certificates. Direct independent
rational rounding of the n12 isotropic shell duals initially failed;
the first frozen intervals were narrower than 1e-9. Subsequent exact
row reduction on the entire active dual face recovered genuine endpoint
certificates. That intermediate failure remains recorded in the JSON
status history rather than being silently omitted.
At n14 the unrestricted primal law itself is an exact code-side dual:
its support lies in the nearcode, and its minimum response against ALL
physical columns is 4837/1768. Transferring that dual also closes the
isotropic game. This self-duality is explicitly replayed without a solver.

For example, the very cheap n12 ground response 9/5 rises above 2.46
already at deficit two under isotropy. Ground-only optimization cannot
silently stand in for a complete nearlevel or all-energy-shell theorem.
Likewise T=2 and T=4 coinciding in some rows is a property of the finite
energy spectrum, not an asserted continuity in the energy window.

## 4. Reproducibility and scope

The full-cube LP and rational replay code is
`computations/paper_discrepancy_2026_09_17_full_column_games.py`.
Run it with `--orders 4,6,8,10,12,14 --full`.
The additional exact affine-face recovery and self-dual replay are
`computations/paper_discrepancy_2026_09_17_column_dual_face.py` and
`computations/paper_discrepancy_2026_09_17_column_self_dual_probe.py`.
All 42 full rational primal laws and all-column lower certificates are
frozen at
`tmp/paper_portfolio_2026_09_17/discrepancy/full_column_game_certificates.json`.

The radial covariance optima and their separate response constraints
are detailed in `artifacts/paper_discrepancy_radial_correlation_2026_09_17.md`.
The finite games concern scalar L1 response only. They do not by
themselves give a sharp subGaussian law, independent-column realization,
control of off-code energies, a favorable Gaussian child value, or an
asymptotic extension of actual minimizers.
