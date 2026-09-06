"""Exact finite checks for causal_signed_flow_r13.md."""

from itertools import product
from math import sqrt
from random import Random

import numpy as np


def q_norm(matrix: np.ndarray) -> int:
    spins = np.asarray(list(product((-1, 1), repeat=len(matrix))), dtype=int)
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins)
    return int(np.max(np.abs(energies)))


A7 = np.asarray(
    [
        [0, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, 1, -1, -1, 1],
        [1, 1, 0, 1, -1, 1, -1],
        [1, 1, 1, 0, 1, -1, -1],
        [1, -1, -1, 1, 0, -1, -1],
        [1, -1, 1, -1, -1, 0, -1],
        [1, 1, -1, -1, -1, -1, 0],
    ],
    dtype=int,
)

A5 = np.asarray(
    [
        [0, -1, 1, -1, 1],
        [-1, 0, -1, 1, 1],
        [1, -1, 0, 1, 1],
        [-1, 1, 1, 0, 1],
        [1, 1, 1, 1, 0],
    ],
    dtype=int,
)

assert q_norm(A7) == 18
B = A7[1:, 1:]
assert q_norm(B) == 10
assert all(q_norm(np.delete(np.delete(B, i, 0), i, 1)) == 8 for i in range(6))
assert q_norm(A5) == 8
assert all(q_norm(np.delete(np.delete(A5, i, 0), i, 1)) == 8 for i in range(5))

alpha7 = 18 / (7 * sqrt(7))
g7 = 8 - alpha7 * 5 ** 1.5
v6 = 10 - alpha7 * 6 ** 1.5
assert g7 < 0 < g7 - v6

x7 = np.asarray((1, -1, -1, -1, -1, -1, -1))
y6 = np.asarray((1, -1, -1, -1, -1, -1))
assert x7 @ A7 @ x7 == -18
assert tuple((np.diag(x7) @ (-A7) @ np.diag(x7)).sum(axis=1)) == (6, 0, 0, 0, 4, 4, 4)
assert y6 @ B @ y6 == -10
assert tuple((np.diag(y6) @ (-B) @ np.diag(y6)).sum(axis=1)) == (1, 1, 1, 1, 1, 5)

alpha5 = 8 / (5 * sqrt(5))
prefix5 = 8 - alpha5 * 4 ** 1.5
g5 = 2 - alpha5 * 2 ** 1.5
assert g5 < 0 < prefix5


# Exhaustively compare general Hall deficiency with the antichain formula on
# 1,000 random capacity assignments on a fixed branching tree.
parents = (-1, 0, 0, 1, 1, 2, 2)
n_tree = len(parents)
desc = [set() for _ in range(n_tree)]

# Build descendant sets without any recursion-dependent ordering assumptions.
for v in range(n_tree):
    desc[v] = {v}
    changed = True
    while changed:
        changed = False
        for w in range(n_tree):
            if parents[w] in desc[v] and w not in desc[v]:
                desc[v].add(w)
                changed = True

comparable = [
    [w in desc[v] or v in desc[w] for w in range(n_tree)]
    for v in range(n_tree)
]
rng = Random(13036925)
for _ in range(1_000):
    demand = [rng.randrange(5) for _ in range(n_tree)]
    supply = [rng.randrange(5) for _ in range(n_tree)]
    hall = 0
    for mask in range(1 << n_tree):
        neighborhood = set()
        demand_sum = 0
        for v in range(n_tree):
            if mask >> v & 1:
                demand_sum += demand[v]
                neighborhood |= desc[v]
        hall = max(hall, demand_sum - sum(supply[w] for w in neighborhood))

    h = [
        sum(demand[w] - supply[w] for w in desc[v])
        for v in range(n_tree)
    ]
    antichain = 0
    for mask in range(1 << n_tree):
        chosen = [v for v in range(n_tree) if mask >> v & 1]
        if all(not comparable[v][w] for i, v in enumerate(chosen) for w in chosen[i + 1 :]):
            antichain = max(antichain, sum(h[v] for v in chosen))
    assert hall == antichain

print(
    "PASS (including 1,000 exact Hall/antichain trials):",
    f"A7 two-step G={g7:.12f}, V6={v6:.12f}, suffix x2={g7-v6:.12f};",
    f"A5 prefix={prefix5:.12f}, G={g5:.12f}",
)
