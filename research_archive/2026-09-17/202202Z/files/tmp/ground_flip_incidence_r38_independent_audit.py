#!/usr/bin/env python3
"""Independent exact audit of Wave 38 ground-flip incidence claims."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9
from ground_flip_incidence_r38_check import fractional_cover


def spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def cap(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in spins(len(a)))


def flipped(a: np.ndarray, selector: tuple[int, ...]) -> np.ndarray:
    inside = np.zeros(len(a), dtype=bool)
    inside[list(selector)] = True
    keep = np.outer(inside, inside)
    out = np.where(keep, a, -a)
    np.fill_diagonal(out, 0)
    return out


def cut_weight(edges: np.ndarray, mask: np.ndarray) -> int:
    crossing = np.logical_xor.outer(mask, mask)
    return int(np.sum(edges[crossing]) // 2)


def audit(a: np.ndarray, m: int, name: str) -> None:
    n = len(a)
    q = cap(a)
    selectors = list(itertools.combinations(range(n), m))
    states = [(x, sigma) for x in spins(n) for sigma in (-1, 1)]
    rows = np.asarray([int((x * (a @ x)) @ (x * (a @ x))) for x, _ in states])
    incidence = np.zeros((len(selectors), len(states)), dtype=np.int8)
    ground = np.zeros_like(incidence)
    max_factor_error = 0
    edge_i, edge_j = np.triu_indices(n, 1)
    subset_masks = np.asarray(
        [[(mask_bits >> i) & 1 for i in range(n)] for mask_bits in range(1 << n)],
        dtype=bool,
    )
    crossing = np.logical_xor(
        subset_masks[:, edge_i], subset_masks[:, edge_j]
    ).astype(np.int64)

    for si, selector in enumerate(selectors):
        selector_set = set(selector)
        b = flipped(a, selector)
        qb = cap(b)
        delta = qb - q
        assert delta >= 0
        for di, (x, sigma) in enumerate(states):
            s = sigma * np.outer(x, x) * a
            t = sigma * np.outer(x, x) * b
            parent_energy = int(s.sum())
            switched_energy = int(t.sum())
            parent_deficit = q - parent_energy
            switched_deficit = qb - switched_energy
            assert parent_deficit >= 0 and switched_deficit >= 0
            incidence[si, di] = switched_energy >= q
            ground[si, di] = switched_energy == qb

            parent_degrees = np.sum(s, axis=1)
            if incidence[si, di]:
                assert switched_deficit <= delta
                assert set(np.flatnonzero(parent_degrees > delta / 4)).issubset(selector_set)

            inside = np.zeros(n, dtype=bool)
            inside[list(selector)] = True
            svec = s[edge_i, edge_j]
            tvec = t[edge_i, edge_j]
            internal_edges = inside[edge_i] & inside[edge_j]
            cs = crossing @ svec
            ct = crossing @ tvec
            cs_internal = crossing @ (svec * internal_edges)
            max_factor_error = max(
                max_factor_error, int(np.max(np.abs(ct - (2 * cs_internal - cs))))
            )
            assert np.all(cs >= -parent_deficit / 4)
            assert np.all(ct >= -switched_deficit / 4)
            omitted_subsets = ~np.any(subset_masks & inside[None, :], axis=1)
            assert np.all(cs[omitted_subsets] >= -parent_deficit / 4)
            assert np.all(cs[omitted_subsets] <= switched_deficit / 4)
            if incidence[si, di]:
                assert np.all(cs[omitted_subsets] <= delta / 4)
            if ground[si, di]:
                assert np.all(ct >= 0)

    assert max_factor_error == 0
    assert np.all(ground <= incidence)
    assert np.all(ground.sum(axis=1) > 0)

    degree = ground.sum(axis=0)
    for di, (x, sigma) in enumerate(states):
        parent_degrees = sigma * x * (a @ x)
        p = int(np.sum(parent_degrees > 0))
        upper = math.comb(n - p, m - p) if p <= m else 0
        assert int(degree[di]) <= upper

    huge = 10**18
    min_g = np.min(np.where(ground, rows[None, :], huge), axis=1)
    min_i = np.min(np.where(incidence, rows[None, :], huge), axis=1)
    tau_g = fractional_cover(ground)
    tau_i = fractional_cover(incidence)
    print(
        name,
        {
            "q": q,
            "ground_min_hist": {int(v): int(np.sum(min_g == v)) for v in np.unique(min_g)},
            "incidence_min_hist": {int(v): int(np.sum(min_i == v)) for v in np.unique(min_i)},
            "mean_cap_ground_misses": int(np.sum(min_g > n * (n - 1))),
            "mean_cap_incidence_misses": int(np.sum(min_i > n * (n - 1))),
            "tau_ground_numeric": tau_g,
            "tau_incidence_numeric": tau_i,
        },
    )


def main() -> None:
    audit(A6, 5, "A6")
    audit(A8, 6, "A8")
    audit(A9, 7, "A9")
    print("PASS ground_flip_incidence_r38_independent_audit")


if __name__ == "__main__":
    main()
