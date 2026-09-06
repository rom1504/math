# Stored optimizer nesting and exact extension replay, orders 10--14

Status: exact exhaustive spin/row replay of the specified stored matrices.
Previously stored solver lower certificates are referenced, not rerun here.

| Stored matrix | n | cap | Optimal one-vertex predecessors | Best extension | E |
|---|---:|---:|---:|---:|---:|
| exact_m10 | 10 | 13 | all 10 | 19 | 6 |
| first ten vertices of nested_10_in_11_cap17 | 10 | 13 | all 10 | 17 | 4 |
| nested_10_in_11_cap17 | 11 | 17 | exactly 1 | 18 | 1 |
| heuristic_m11 (now known optimal) | 11 | 17 | none | 20 | 3 |
| extension_nested_m11_to_12 | 12 | 18 | all 12 | 24 | 6 |
| heuristic_m12 (NOT optimal) | 12 | 20 | none | 22 | 2 |
| bridge_6_7_sign1_cap20 | 13 | 20 | none | 21 | 1 |
| conference_completion_m13 | 14 | 21 | all 14 | 27 | 6 |

The relevant known optimum values are M_9=12, M_10=13, M_11=17, M_12=18,
M_13=20, M_14=21. Exactness of values 11--14 depends on the repository's
stored solver-certified lower conclusions, with their documented lack of a
standalone lower proof object; the matrix caps/deletions/extensions here are
direct exhaustive integer checks.

In particular there DOES exist an optimal order-11 matrix with an optimal
order-10 child. The stored nested witness supplies it. There also exists an
optimal order-11 matrix without ANY optimal predecessor, so universal
heredity of exact minimizers is false even at this order. That fact does not
falsify selectable-parent heredity.

The first two order-10 matrices have identical optimal predecessor counts but
different insertion costs. Predecessor optimality alone cannot select a good
row-extension parent. The stored order-13 optimal witness has no optimal
order-12 child, but the calculation does not quantify over all order-13
optimizers and therefore does not prove universal nonnesting.

Reproduction and source hashes are in
`computations/decisive_independent_nesting_audit_2026_09_06.py` and its JSON
output. The script exhaustively computes each parent cap, every one-vertex
deletion cap, and all 2^(n-1) inequivalent extension rows.
