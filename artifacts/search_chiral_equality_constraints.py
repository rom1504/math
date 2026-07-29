#!/usr/bin/env python3
"""Necessary-condition MILP for equality compressed four-lifts."""

from __future__ import annotations

from dataclasses import dataclass
import itertools

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

from audit_dependent_4lift import B_MINUS, B_TRIANGLE, exact_q


@dataclass(frozen=True)
class EdgeIndex:
    order: int
    edges: list[tuple[int, int]]
    index: dict[tuple[int, int], int]

    @classmethod
    def create(cls, order: int) -> "EdgeIndex":
        edges = list(itertools.combinations(range(order), 2))
        return cls(order, edges, {edge: k for k, edge in enumerate(edges)})

    def at(self, u: int, v: int) -> int:
        if u > v:
            u, v = v, u
        return self.index[(u, v)]


def seed_extremizers(seed: np.ndarray) -> tuple[int, list[np.ndarray], list[np.ndarray]]:
    q, _, _ = exact_q(seed)
    n = len(seed)
    positives = []
    negatives = []
    for bits in itertools.product((-1, 1), repeat=n - 1):
        x = np.asarray((1,) + bits, dtype=np.int8)
        e = int(x @ seed @ x)
        if e == q:
            positives.append(x)
        elif e == -q:
            negatives.append(x)
    return q, positives, negatives


def necessary_system(seed: np.ndarray, fibre: int = 4):
    n = len(seed)
    order = n * fibre
    edge_index = EdgeIndex.create(order)
    q, positives, negatives = seed_extremizers(seed)

    rows: list[dict[int, int]] = []
    lower: list[int] = []
    upper: list[int] = []
    labels: list[str] = []

    def add(coeff: dict[int, int], lo: int, hi: int, label: str) -> None:
        rows.append(coeff)
        lower.append(lo)
        upper.append(hi)
        labels.append(label)

    # Cross-block compression: sum b_uv = 8 a_ij. With b=2z-1,
    # sum z_uv = 12 for a=+1 and 4 for a=-1.
    for i, j in itertools.combinations(range(n), 2):
        coeff = {}
        for a in range(fibre):
            for b in range(fibre):
                k = edge_index.at(fibre * i + a, fibre * j + b)
                coeff[k] = 1
        target = 12 if seed[i, j] == 1 else 4
        add(coeff, target, target, f"compress({i},{j})")

    # Equality on centered/chiral seed forces total diagonal-fibre
    # energy d_B=0, equivalently sum of all within-fibre edge signs 0.
    # In z variables this is sum z = (# within edges)/2.
    coeff = {}
    for i in range(n):
        for a, b in itertools.combinations(range(fibre), 2):
            k = edge_index.at(fibre * i + a, fibre * i + b)
            coeff[k] = 1
    add(coeff, len(coeff) // 2, len(coeff) // 2, "zero diagonal-fibre shift")

    # A repeated positive maximizer must be a global maximum, hence
    # every switched local field is >= 1 (degree is odd). A repeated
    # negative maximizer must have every switched local field <= -1.
    for energy_sign, states in ((1, positives), (-1, negatives)):
        for state_id, x in enumerate(states):
            X = np.repeat(x, fibre)
            for u in range(order):
                coeff = {}
                constant = 0
                for v in range(order):
                    if u == v:
                        continue
                    sign = int(X[u] * X[v])
                    k = edge_index.at(u, v)
                    # b_k=2z_k-1.
                    coeff[k] = 2 * sign
                    constant -= sign
                if energy_sign == 1:
                    add(
                        coeff,
                        1 - constant,
                        np.inf,
                        f"positive local field state={state_id} vertex={u}",
                    )
                else:
                    add(
                        coeff,
                        -np.inf,
                        -1 - constant,
                        f"negative local field state={state_id} vertex={u}",
                    )

    matrix = lil_matrix((len(rows), len(edge_index.edges)), dtype=np.float64)
    for i, row in enumerate(rows):
        for j, value in row.items():
            matrix[i, j] = value
    return (
        q,
        positives,
        negatives,
        edge_index,
        LinearConstraint(matrix.tocsr(), np.asarray(lower), np.asarray(upper)),
        labels,
    )


def solve(name: str, seed: np.ndarray) -> None:
    q, positives, negatives, edge_index, constraints, labels = necessary_system(seed)
    print(
        name,
        "seed order",
        len(seed),
        "Q",
        q,
        "positive/negative extremizers",
        len(positives),
        len(negatives),
        "variables",
        len(edge_index.edges),
        "constraints",
        len(labels),
    )
    result = milp(
        c=np.zeros(len(edge_index.edges)),
        integrality=np.ones(len(edge_index.edges)),
        bounds=Bounds(0, 1),
        constraints=constraints,
        options={"time_limit": 300, "presolve": True},
    )
    print("status", result.status, result.message)
    if result.x is not None:
        print("feasible objective", result.fun)


def main() -> None:
    solve("B_MINUS", B_MINUS)
    solve("B_TRIANGLE", B_TRIANGLE)


if __name__ == "__main__":
    main()
