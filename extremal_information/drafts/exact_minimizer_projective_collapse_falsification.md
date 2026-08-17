# Exact-minimizer projective-shell collapse: finite falsification study

Status: solver-certified finite computation plus independently enumerable
explicit witnesses; no standalone UNSAT proof objects.  This note is a finite
falsification report, not asymptotic evidence.

## Question and exact model

For an order-`n` signing `A`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
d_{\mathrm{proj}}(x,y)=\min\{d_H(x,y),n-d_H(x,y)\}.
```

After choosing an absolute ground state, switching it to the all-one word,
and globally negating `A` if necessary, the search asks whether an exact-cap
signing can satisfy

```math
|H_A(x)|\le M_n-d-2
\quad\text{whenever}\quad
d_{\mathrm{proj}}(x,\mathbf 1)>R.                 \tag{PC.1}
```

Thus its entire absolute deficit-`d` shell

```math
\{x:|H_A(x)|\ge M_n-d\}
```

lies in the projective vertex ball of radius `R` about the declared ground.
Only even `d` are used because every energy has the parity of
`binom(n,2)`.

The exact input values are

```text
(M_7,...,M_14)=(9,10,12,13,17,18,20,21).
```

They are recorded in `ACTIVE_STATE.md`.  Their finite certificate sources
are the order-specific `exact_m*.json` files through order 10,
`computations/results/certified_m11_m12.json`, and
`computations/results/certified_m13_m14.json`.

The implementation is
[`search_exact_minimizer_projective_collapse.py`](../experiments/search_exact_minimizer_projective_collapse.py).
It represents one spin from each antipodal line, constrains every such energy
to `[-M_n,M_n]`, and fixes the all-one energy to `M_n`.  Therefore every
feasible output has exact cap `M_n`, not merely cap at most `M_n`.

The gauge and symmetry reductions lose no candidates:

1. diagonal switching sends any selected absolute ground to all-one;
2. global coefficient negation makes its energy positive;
3. a vertex permutation sorts the negative degrees while preserving all-one
   and all projective distance constraints.

For every feasible result the saved matrix is exhaustively re-enumerated over
all `2^(n-1)` projective spins.  The output records the full absolute-energy
histogram and the radii of every requested shell state, and asserts that the
cap and radius are correct.  An `INFEASIBLE` status is an OR-Tools CP-SAT
certificate internal to the solver; no DRAT-like proof object is emitted.

## Decisive results

Let `R_min(n,d)` be the least `R` for which some exact-cap signing and some
absolute ground center have their entire deficit-`d` shell in that ball.
The maximum possible projective radius is `floor(n/2)`, so that value is
always a trivial upper bound once existence of an exact minimizer is known.

| `n` | `d=0` | `d=2` | `d=4` | basis |
|---:|---:|---:|---:|:---|
| 7 | **1** | **3** | **3** | SAT at `R=1` for `d=0`; UNSAT below the displayed values |
| 8 | **4** | **4** | **4** | active shell UNSAT for `R=1,2,3`; larger shells contain it |
| 9 | **4** | **4** | **4** | active shell UNSAT for `R=1,2,3`; larger shells contain it |
| 10 | unknown | **5** | **5** | `d=2,4` UNSAT for `R=1,2,3,4` |
| 11 | unknown | unknown | `2 <= R_min <= 5` | `d=4,R=1` UNSAT; `R=2` timed out |
| 12 | unknown | unknown | unknown | attempted instances timed out |
| 13 | unknown | unknown | unknown | attempted instances timed out |
| 14 | unknown | unknown | unknown | attempted instances timed out |

The order-seven counterexample is explicit.  Its entire absolute active shell
has four projective words at radii

```text
0, 1, 1, 1
```

from a gauged positive ground.  Its cap is exactly `M_7=9`.  Hence exact
minimization alone does **not** force even radius-two active-shell diffusion.
Admitting deficit two changes the same finite order sharply: every exact
minimizer shell then reaches the maximum radius three.

The exclusions at orders 8 and 9 say more than the absence of the particular
order-seven pattern: relative to every possible chosen absolute ground of
every exact minimizer, another absolute ground lies at the maximum projective
distance.  At order 10 this all-the-way-to-the-equator conclusion is certified
for the deficit-two and deficit-four shells, although the active-only question
remains unresolved by the allotted computation.

## Inconclusive runs and interpretation

All `UNKNOWN` statuses are solver timeouts and carry no mathematical
conclusion.  In particular, nothing is claimed for orders 12--14, and the
active-shell question at orders 10--11 is open in this experiment.  The
longest single run was the order-10 active-shell query at `R=4`, which remained
unknown after 180 seconds.

The result is discriminating in both directions:

- the `n=7,d=0` witness kills any finite theorem asserting automatic diffuse
  active shells for every exact minimizer;
- the exact exclusions for `n=8,9` and the `d=2,4` exclusions for `n=10`
  show that the collapse is not a universal small-order artifact and that one
  parity step of shell thickness can force genuinely global projective reach.

These statements do **not** provide a fixed normalized separation as
`n -> infinity`, nor a probability or counting bound.  They should therefore
be used only as finite falsifiers and as motivation for a theorem about thin
shells rather than exact active sets.

## Reproduction

All machine outputs are under
[`experiments/results/`](../experiments/results/) with names

```text
projective_collapse_nNN_dD_rR.json.
```

A representative invocation is

```bash
.venv/bin/python extremal_information/experiments/\
  search_exact_minimizer_projective_collapse.py \
  10 13 --vertex-radius 4 --shell-deficit 2 \
  --time-limit 90 --workers 4 \
  --output extremal_information/experiments/results/\
projective_collapse_n10_d2_r4.json
```

The JSON records Python and OR-Tools versions, worker and time limits, wall
time, branches, conflicts, status, and (when feasible) the matrix plus its
exhaustively verified profile.
