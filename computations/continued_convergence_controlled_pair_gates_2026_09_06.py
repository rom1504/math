"""Exact simultaneous controlled-translation gates on arbitrary tensor seeds.

Default audit exhausts all order-four hollow sign seeds and all compatible
translation assignments on an order-four Walsh outer factor. Integer
quadratic caps are exact, with antipodal configurations quotiented.
"""

import argparse
import itertools
import json

import numpy as np

from continued_convergence_symmetric_switching_2026_09_06 import walsh


def sign_vectors(n):
    codes = np.arange(1 << (n - 1), dtype=np.uint64)
    x = np.ones((len(codes), n), dtype=np.int64)
    x[:, 1:] = 1 - 2 * ((codes[:, None] >> np.arange(n - 1, dtype=np.uint64)) & 1).astype(np.int64)
    return x


def controlled_gate_twice(n, outer_order, directions):
    """Twice an orthogonal controlled-permutation matrix, no floating point."""
    assert n % 2 == 0 and len(directions) == n // 2
    twice_r = np.zeros((n * outer_order, n * outer_order), dtype=np.int64)
    eye = np.eye(outer_order, dtype=np.int64)
    for pair, direction in enumerate(directions):
        perm = eye[np.arange(outer_order) ^ direction]
        block = np.block([[eye + perm, eye - perm], [eye - perm, eye + perm]])
        indices = np.arange(2 * pair * outer_order, 2 * (pair + 1) * outer_order)
        twice_r[np.ix_(indices, indices)] = block
    assert np.array_equal(twice_r @ twice_r.T, 4 * np.eye(n * outer_order, dtype=np.int64))
    return twice_r


def compatible_directions(n, outer_order):
    for directions in itertools.product(range(outer_order), repeat=n // 2):
        if all(bin(a & b).count("1") % 2 == 0 for a, b in itertools.combinations(directions, 2)):
            yield directions


def run_exhaustive(n, outer_order):
    assert n == 4 and outer_order <= 4
    outer = walsh(outer_order)
    edges = list(itertools.combinations(range(n), 2))
    seed_codes = np.arange(1 << len(edges), dtype=np.uint64)
    edge_signs = 1 - 2 * ((seed_codes[:, None] >> np.arange(len(edges), dtype=np.uint64)) & 1).astype(np.int64)
    seeds = np.repeat(np.eye(n, dtype=np.int64)[None, :, :], len(seed_codes), axis=0)
    for column, (i, j) in enumerate(edges):
        seeds[:, i, j] = edge_signs[:, column]
        seeds[:, j, i] = edge_signs[:, column]
    seed_x = sign_vectors(n)
    seed_energies = np.einsum("xi,sij,xj->sx", seed_x, seeds - np.eye(n, dtype=np.int64), seed_x)
    seed_caps = np.max(np.abs(seed_energies), axis=1) // 2
    x = sign_vectors(n * outer_order)
    records = []
    for directions in compatible_directions(n, outer_order):
        twice_r = controlled_gate_twice(n, outer_order, directions)
        # For every seed, certify sign-valuedness and exact similarity.
        for seed in seeds:
            parent = np.kron(seed, outer)
            four_new = twice_r @ parent @ twice_r.T
            assert np.all(np.isin(four_new, (-4, 4)))
            assert np.trace(four_new) == 4 * np.trace(parent)
        twice_y = (x @ twice_r).reshape(len(x), n, outer_order)
        diagonal = np.einsum("xiu,uv,xiv->x", twice_y, outer, twice_y)
        coefficients = np.stack([
            2 * np.einsum("xu,uv,xv->x", twice_y[:, i], outer, twice_y[:, j])
            for i, j in edges
        ])
        four_values = edge_signs @ coefficients + diagonal[None, :]
        assert np.all(four_values % 4 == 0)
        # trace(outer)=0 for the tested orders>1, so these are hollow caps too.
        caps = np.max(np.abs(four_values), axis=1) // 8
        records.append({"directions": list(directions), "caps_by_seed_code": caps.tolist()})
    initial = np.array(records[0]["caps_by_seed_code"])
    all_caps = np.array([record["caps_by_seed_code"] for record in records])
    best = all_caps.min(axis=0)
    worst = all_caps.max(axis=0)
    return {
        "seed_order": n, "outer_order": outer_order,
        "seed_caps": seed_caps.tolist(),
        "initial_parent_caps": initial.tolist(),
        "best_parent_caps": best.tolist(),
        "worst_parent_caps": worst.tolist(),
        "improved_seed_codes": np.flatnonzero(best < initial).tolist(),
        "worsened_seed_codes": np.flatnonzero(worst > initial).tolist(),
        "records": records,
        "scope": "Exact finite caps; neither growing-seed landing nor asymptotic convergence is certified.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-order", type=int, default=4, choices=[2, 4])
    parser.add_argument("--output")
    args = parser.parse_args()
    result = run_exhaustive(4, args.outer_order)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps({key: result[key] for key in (
        "seed_order", "outer_order", "improved_seed_codes", "worsened_seed_codes", "scope"
    )}, indent=2))


if __name__ == "__main__":
    main()
