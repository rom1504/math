# Bounded joint Hadamard-bridge probe with actual optimal order12 children

The existing order8 probe found a cap30 parent with actual cap10 children.
We tested the next available Hadamard child order,12, rather than repeating
that finite search. The stored solver-certified value is M12=18, with the
explicit witness in `computations/results/extension_nested_m11_to_12.json`.
The lower-bound provenance is preserved in `certified_m11_m12.json`; it is
a solver certificate, not a new standalone proof in this note.

## Setup and bounded search

Both internal children were fixed to that actual order12 witness, allowing
both relative global polarities. The cross block was a Paley Hadamard12,
with sampled row/column permutations and optimized row/column sign phases.
Every tested bridge therefore had exact sign entries and orthogonal rows.

For projective child spin representatives X (first coordinate1), the exact
parent cap is

 max_{u,v} |H_A(u)+sigma H_A(v)|+|u^T C v|.

The formula accounts for independent whole-child reversal exactly. There
is no assumption that favorable center cancellation pays for the bridge.

The bounded run used100 permutation restarts and8012 phase-neighbor
evaluations. Accepted candidates were evaluated over all2048² projective
child pairs; rejected candidates could stop as soon as they exceeded the
current cap cutoff, and are not described as exact full profiles.

Result: best cap64. Both the matrix and search history are saved in
`computations/results/flatify_construct_2026_09_07_bridge12_probe.json`.
An independent full-parent Gray-code evaluator checked all8388608
projective order24 spins and found minimum-62, maximum64, cap64. Its
replay is `flatify_construct_2026_09_07_bridge12_verify.json`.

The finite row-normalized child target is36sqrt(23/11)=52.05591399...;
the asymptotic leading proxy is36sqrt2=50.91168825... . Cap52 would meet
the former, but was NOT found. The existing repository conference26
restriction calculation supplies order24 cap56, already better than this
probe. No improved finite upper bound results.

This was not an exhaustive search over Hadamard equivalences, optimal
order12 representatives, or bridges. Its failure does not prove such a
favorable parent is impossible. It does show that the favorable order8
example did not recur in this bounded joint order12 test. No asymptotic
recurrence or convergence claim is made.

## Exact fixed-Paley phase-cube follow-up

A sharper bounded calculation then solved the phase problem EXACTLY for
one fixed Paley12 row/column ordering and the same fixed optimal child.
For EACH relative child polarity, the minimum over all2048² projective
row/column phase pairs is62.

This is certified without evaluating all phase/spin quadruples. If K_ij
is the absolute Hadamard interaction of projective spin representatives,
a violating spin pair(i,j) imposes the necessary phase constraint

 |e(i xor p)+sigma e(j xor q)|+K_ij<=T.

At target T=60,1092 such constraints for sigma=+1 and1087 for sigma=-1
exclude the entire phase cube. At T=62, explicit feasible phase pairs
were found for both polarities. Since parent energies are even, exclusion
at60 and feasibility at62 prove the exact minimum.

Search data and all exclusion constraints are in
`flatify_construct_2026_09_07_bridge12_phase_certificate.json`. A separate
verifier recomputed every constraint's Hadamard entry and replayed all
phase exclusions, then independently enumerated all8388608 projective
full-parent spins of BOTH62 witnesses. Both have minimum-62 and maximum62.
The replay is `flatify_construct_2026_09_07_bridge12_phase_verify.json`.

This improves the heuristic64 result but remains non-favorable. The exact
minimum statement is ONLY for the declared fixed row/column ordering and
fixed optimal child, not for arbitrary Hadamard permutations or arbitrary
selectable optimal12 children.

## Why splitting a rounding into many steps is not itself a repair

The accompanying global-operation question asked whether correlated
rounding increments could cancel their leading discrepancy. For the
original block-diagonal weighted target B, ANY final full-sign matrix C
has cross-block discrepancy entries of magnitude1. Thus each row of
E=C-B has Euclidean norm at least the square root of the opposite block
size, regardless of how all old edges were changed or how the increments
were correlated. The already audited Khintchine/polarization inequality

 Q(E)>=(1/(4sqrt2))sum_i ||E_i||2

therefore forces Q(C-B)=Omega(N^(3/2)) on comparable splits. A multistep
construction cannot evade this endpoint statement by cancelling stochastic
variances along the path.

This ONLY excludes a uniform-discrepancy justification against B. The
requested inequality compares Q(C) with Q(B), not Q(C-B), so a genuinely
one-sided favorable correlation with the whole energy landscape could
still succeed. Orthogonally moving the target also changes that landscape
and needs its own paid cap estimate. We did not silently convert this
endpoint discrepancy observation into an impossibility theorem for the
authorized global old-edge operation.
