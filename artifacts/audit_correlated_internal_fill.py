#!/usr/bin/env python3
"""Exact small-order audit of cap-kernel and thick-band internal fills.

The script enumerates switching classes through n=7, selects one
representative of every endpoint profile of every centered-width
minimizer, and counts:

* sign vectors satisfying all exact-cap equations;
* sign vectors satisfying every thick-cap band;
* the minimum absolute midpoint among both sets.

All Hamiltonians and constraints use exact integer arithmetic.
"""

from __future__ import annotations

import argparse
import itertools
from collections import OrderedDict

import numpy as np


def spin_table(n: int) -> np.ndarray:
    states = 1 << (n - 1)
    x = np.ones((states, n), dtype=np.int8)
    for state in range(states):
        for i in range(n - 1):
            if (state >> i) & 1:
                x[state, i] = -1
    return x


def edge_table(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def signing_from_mask(
    n: int, free_edges: list[tuple[int, int]], mask: int
) -> np.ndarray:
    a = np.zeros((n, n), dtype=np.int8)
    for i in range(n - 1):
        a[i, n - 1] = a[n - 1, i] = 1
    for e, (i, j) in enumerate(free_edges):
        a[i, j] = a[j, i] = -1 if ((mask >> e) & 1) else 1
    return a


def find_minimizers(n: int, spins: np.ndarray) -> tuple[int, list[int]]:
    free_edges = [
        (i, j) for i in range(n - 1) for j in range(i + 1, n - 1)
    ]
    all_edges = edge_table(n)
    chi = np.stack(
        [spins[:, i] * spins[:, j] for i, j in all_edges], axis=1
    ).astype(np.int16)
    best = 10**9
    masks: list[int] = []
    for mask in range(1 << len(free_edges)):
        a = signing_from_mask(n, free_edges, mask)
        coeff = np.array([a[i, j] for i, j in all_edges], dtype=np.int16)
        energy = chi @ coeff
        width = int((energy.max() - energy.min()) // 2)
        if width < best:
            best = width
            masks = [mask]
        elif width == best:
            masks.append(mask)
    return best, masks


def enumerate_sign_vectors(m: int) -> np.ndarray:
    count = 1 << m
    g = np.ones((count, m), dtype=np.int8)
    for mask in range(count):
        for e in range(m):
            if (mask >> e) & 1:
                g[mask, e] = -1
    return g


def audit_n(n: int) -> list[dict[str, int]]:
    spins = spin_table(n)
    all_edges = edge_table(n)
    free_edges = [
        (i, j) for i in range(n - 1) for j in range(i + 1, n - 1)
    ]
    chi = np.stack(
        [spins[:, i] * spins[:, j] for i, j in all_edges], axis=1
    ).astype(np.int16)
    width, minimizers = find_minimizers(n, spins)

    representatives: OrderedDict[tuple[int, ...], dict[str, int]] = (
        OrderedDict()
    )
    for mask in minimizers:
        a = signing_from_mask(n, free_edges, mask)
        coeff = np.array([a[i, j] for i, j in all_edges], dtype=np.int16)
        energy = chi @ coeff
        top = int(energy.max())
        bottom = int(energy.min())
        tops = np.flatnonzero(energy == top)
        bottoms = np.flatnonzero(energy == bottom)

        for ts, bs in itertools.product(tops, bottoms):
            xt = spins[ts]
            switched = a * np.outer(xt, xt).astype(np.int8)
            relative_bottom = xt * spins[bs]
            side = relative_bottom == -1
            usize = int(side.sum())
            vsize = n - usize

            internal = [
                (i, j) for i, j in all_edges if side[i] == side[j]
            ]
            cross = [(i, j) for i, j in all_edges if side[i] != side[j]]
            m = len(internal)

            delta = np.stack(
                [
                    ((1 - spins[:, i] * spins[:, j]) // 2)
                    for i, j in internal
                ],
                axis=1,
            ).astype(np.int16)
            cross_energy = np.zeros(len(spins), dtype=np.int16)
            for i, j in cross:
                cross_energy += (
                    switched[i, j] * spins[:, i] * spins[:, j]
                )
            cap = np.abs(cross_energy) == width
            cap_rows = delta[cap]
            rank = int(np.linalg.matrix_rank(cap_rows.astype(float)))
            midpoint = (top + bottom) // 2
            profile = (
                abs(midpoint),
                min(usize, vsize),
                max(usize, vsize),
                m,
                rank,
                int(cap.sum()),
            )
            if profile in representatives:
                continue

            g = enumerate_sign_vectors(m).astype(np.int16)
            cap_ok = np.all(cap_rows @ g.T == 0, axis=0)
            bands = width - np.abs(cross_energy)
            band_ok = np.all(
                2 * np.abs(delta @ g.T) <= bands[:, None], axis=0
            )

            cap_indices = np.flatnonzero(cap_ok)
            band_indices = np.flatnonzero(band_ok)
            cap_mid = (
                int(np.min(np.abs(g[cap_indices].sum(axis=1))))
                if len(cap_indices)
                else -1
            )
            band_mid = (
                int(np.min(np.abs(g[band_indices].sum(axis=1))))
                if len(band_indices)
                else -1
            )
            original_h = np.array(
                [switched[i, j] for i, j in internal], dtype=np.int16
            )
            original_ok = bool(
                np.all(2 * np.abs(delta @ original_h) <= bands)
            )
            representatives[profile] = {
                "n": n,
                "W": width,
                "midpoint": midpoint,
                "u": min(usize, vsize),
                "v": max(usize, vsize),
                "m": m,
                "cap_rank": rank,
                "cap_rows": int(cap.sum()),
                "cap_sign_count": int(cap_ok.sum()),
                "cap_min_abs_midpoint": cap_mid,
                "band_sign_count": int(band_ok.sum()),
                "band_min_abs_midpoint": band_mid,
                "original_h_band_ok": int(original_ok),
            }
    return list(representatives.values())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=7)
    args = parser.parse_args()
    for n in range(3, args.max_n + 1):
        for row in audit_n(n):
            print(" ".join(f"{key}={value}" for key, value in row.items()))


if __name__ == "__main__":
    main()
