# Exact minimizer optimized-bridge response protocol

## Question

For an exact-minimizer class representative `A` of order `n` and a hollow
sign query `C` of order `k`, compute exactly

```math
F_C(A)=\min_{B\in\{\pm1\}^{n\times k}}
Q\!\begin{pmatrix}A&B\\B^\mathsf T&C\end{pmatrix},
\qquad
Q(D)=\max_z\left|\sum_{i<j}d_{ij}z_i z_j\right|.
```

The finite question is whether the vector of these responses distinguishes
inequivalent exact-minimizer classes after switching gauge and the global-sign
identification in the authoritative class files are removed correctly.

## Authoritative inputs

The run consumes the exhaustive signed-permutation-and-global-sign orbit
classifications

```text
computations/results/m3_minimizer_orbits.json
...
computations/results/m8_minimizer_orbits.json
```

and the certified exact values `M_3,...,M_11` in
`computations/results/exact_m3.json` through `exact_m10.json` and
`computations/results/certified_m11_m12.json`.  SHA-256 hashes of every input
are embedded in the output.  In particular, byte-distinct matrices are never
treated as distinct classes; only the exhaustive class records are used.

## Query quotient

Up to diagonal switching and permutation there is one hollow sign query at
orders one and two.  At order three there are exactly two, classified by

```math
\tau(C)=c_{12}c_{13}c_{23}\in\{+1,-1\}.
```

The script exhaustively verifies these orbit counts.  Independent child
switches transport bridges by `B -> SBT`, giving a bijection of the complete
bridge fibre and preserving every parent energy.  The run checks the block
matrix identity for every projective pair `(S,T)` for every saved witness.

The source classes also quotient global matrix negation.  Since

```math
F_C(-A)=F_{-C}(A),
```

the two order-three coordinates swap under this residual action.  The class
signature is therefore

```text
(F_k=1, F_k=2, sorted(F_tau=+1, F_tau=-1)).
```

This is the comparison made in the result; an oriented pair is also retained.

## Exact algorithms

For `k=1,2`, all `2^(nk)` sign bridges are evaluated in integer NumPy chunks.
All projective parent spins are included.  The output records both the minimum
and the number of optimal labelled bridges.

For `k=3`, the bridge entries are Boolean CP-SAT variables.  For every parent
spin, the two exact linear inequalities imposing cap at most `T` are included.
The run starts at the certified global lower bound `M_(n+3)` and increases `T`
by the forced parity step two.  An `INFEASIBLE` status excludes a target; a
saved feasible bridge is exhaustively re-evaluated over parent spins.  Search
uses one worker and random seed zero.  A valid symmetry reduction restricts
the first bridge row to orbit representatives under signed-permutation
automorphisms of `C` and the exact symmetry `B -> -B`.

Thus the `k=3` classification is solver-certified exact finite computation,
with full solver statistics but no standalone proof object.  The result does
not promote those computations to a paper proof.

## Reproduction

From repository root, using the project virtual environment:

```bash
.venv/bin/python -u \
  extremal_information/experiments/exact_minimizer_optimized_bridge_response.py
```

The output is

```text
extremal_information/experiments/results/
  exact_minimizer_optimized_bridge_response.json
```

The expected response partitions are:

| `n` | canonical response signature | exact-minimizer classes |
|---:|---:|---:|
| 3 | `[4,4,5,7]` | `[0]` |
| 4 | `[4,5,9,9]` | `[0]` |
| 5 | `[5,9,10,10]` | `[0]` |
| 6 | `[9,10,12,12]` | `[0]` |
| 7 | `[10,12,13,15]` | `[0,2]` |
| 7 | `[12,12,13,15]` | `[1]` |
| 8 | `[12,15,17,17]` | `[0]` |
| 8 | `[12,15,19,19]` | `[1]` |

## Scope

This is a finite falsifier/ranking experiment.  It proves that optimized
future responses can retain information beyond the child cap and switching
class, including an order-eight separation first visible here at query order
three.  It gives no family of separated classes, no lower bound growing with
`n`, and no asymptotic response-gap theorem.
