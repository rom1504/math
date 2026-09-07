# A concrete dense bridge between actual minimizing children

Status: exact finite positive witness, with bounded unsuccessful improvement
searches. No asymptotic recurrence or universal bridge theorem is claimed.

## Operation and exact check

For hollow sign children A,D and any sign bridge B, the parent C has the
exact cap identity

\[
 Q\begin{pmatrix}A&B\\B^T&D\end{pmatrix}
 =\max_{x,y}\left(|H_A(x)+H_D(y)|+|x^TBy|\right).
\]

Indeed, replacing every spin in one child by its negative leaves both
internal energies unchanged and reverses the bridge energy. This is the
absolute-cap identity, not a width replacement.

Search over relative row/column sign phases of a Sylvester Hadamard bridge,
both child polarities, all pairs of the graph-atlas actual minimizing
representatives, and a bounded additional permutation sample gives:

| Child orders | Actual child caps | Parent cap | Weighted equal-child target |
|---|---:|---:|---:|
| 4 + 4 | 4, 4 | 10 | 8 sqrt(7/3) = 12.2202... |
| 8 + 8 | 10, 10 | 30 | 20 sqrt(15/7) = 29.2770... |

The order-16 normalized cap is 30/64 = .46875. Thus a **single Hadamard
bridge between two actual optimal children is not subject to a universal
one-half floor**. This does not contradict floors for iterated self-weaves,
global Hadamard spectral classes, MUB amplification, or iid bridges. The
order-16 parent is not in the global Hadamard spectral class.

Full child, bridge, and parent matrices are saved in the linked result
files; no reconstruction from a random seed is needed. The independent
checker evaluates every projective parent spin (32,768 at order 16), also
checks the block identity, child caps, sign constraints, and BB^T=8I.

## The decisive shell is not the child ground-state shell

For the saved order-16 example, both children are the atlas-1005 order-eight
representative. The selected bridge has bilinear cap 20. The exact shell
envelope below is over all 16,384 projective child spin pairs.

| Absolute combined child energy | Maximum absolute bridge energy |
|---:|---:|
| 0 | 20 |
| 2 | 16 |
| 4 | 20 |
| 6 | 16 |
| 8 | 20 |
| 10 | 16 |
| 12 | 16 |
| 14 | 16 |
| 16 | 12 |
| 18 | 8 |
| 20 | 4 |

Exactly 50 projective pairs attain parent cap 30, all in the combined-child
energy-14 shell with bridge magnitude 16. At child energy 20 the bridge is
at most 4. This finite operation therefore pays attention to intermediate
child shells, not only paired ground states. An asymptotic construction
would need control across every shell.

## Exact local trade information and bounded failures

With both order-eight children frozen, every single bridge-edge flip raises
the saved cap to 32. Of all 2,016 double-edge flips, 2,006 give cap 34,
six give cap 32, and four retain cap 30. Exploring all cap-nonincreasing
double-edge moves from this bridge exhausts an eleven-state component;
no member has cap 28. This is an exact finite component computation, not
global infeasibility and not a statement about longer coordinated trades.

A separate Hadamard-preserving bridge search (row/column signs, swaps,
closed-quad trades) tried 48,000 proposals and retained cap 30. Unlike the
constructive track's whole-matrix trades, these keep both children fixed.
The unrestricted cap-28 MILP timed out after 100 seconds without a witness;
CP-SAT returned UNKNOWN after 120 seconds. Neither timeout certifies
infeasibility. Cap 28 remains unresolved for these frozen children.

This is distinct from, but compatible with, the archived Hadamard histogram
entropy discussion in `augmented_cut_code_second_order_amalgamation.md`:
that argument does not supply the required extreme-shell alignment. The
present finite witness supplies one explicit alignment, not its scalable
replacement.

## Reproduction

- `computations/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.py`
  and `computations/results/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json`.
- Independent integer checker:
  `computations/flatify_adversary_2026_09_07_hadamard_bridge_witness_check.py`
  and its same-stem result JSON.
- `computations/flatify_adversary_2026_09_07_frozen_hadamard_trades.py`
  and its same-stem result JSON.
- `computations/flatify_adversary_2026_09_07_neutral_bridge_trade_graph.py`
  and its same-stem result JSON (all eleven states and predecessors saved).
- `computations/flatify_adversary_2026_09_07_frozen_bridge_milp.py` and
  `computations/flatify_adversary_2026_09_07_frozen_bridge_cpsat.py`, with
  same-stem result JSONs preserving the unsuccessful solver statuses.
