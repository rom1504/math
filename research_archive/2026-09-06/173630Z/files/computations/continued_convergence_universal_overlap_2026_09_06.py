"""Universal overlapping reflector search with exact common-quotient audit.

Every selected move preserves signs for every constant-diagonal full seed
of the requested order. Common quotient witnesses are allowed to pull back
to non-Boolean outer eigenvectors: only testing the original Boolean outer
eigenvectors would miss births.
"""

import argparse
import itertools
import json
import math

import numpy as np
import sympy as sp

from continued_convergence_controlled_pair_gates_2026_09_06 import sign_vectors
from continued_convergence_symmetric_switching_2026_09_06 import walsh


def all_seed_parents(n, s):
    edges = list(itertools.combinations(range(n), 2))
    seeds = []
    for signs in itertools.product((1, -1), repeat=len(edges)):
        b = np.eye(n, dtype=np.int64)
        for (i, j), sign in zip(edges, signs):
            b[i, j] = b[j, i] = sign
        seeds.append(b)
    return np.array([np.kron(b, walsh(s)) for b in seeds])


def universal_moves(family):
    n = family.shape[-1]
    # Pair-product hashing finds all closed quadruples without an O(n^4)
    # scan. A row descriptor includes every seed and every column.
    descriptors = family.transpose(1, 0, 2).reshape(n, -1)
    groups = {}
    for a, b in itertools.combinations(range(n), 2):
        descriptor = descriptors[a] * descriptors[b]
        descriptor *= descriptor[0]
        key = np.packbits(descriptor > 0).tobytes()
        groups.setdefault(key, []).append((a, b))
    quads = set()
    for pairs in groups.values():
        for first, second in itertools.combinations(pairs, 2):
            joined = first + second
            if len(set(joined)) == 4:
                quads.add(tuple(sorted(joined)))
    for quad in sorted(quads):
        products = np.prod(family[:, quad, :], axis=1)
        if not np.all(products == products[0, 0]):
            continue
        parity = int(products[0, 0])
        for tail in itertools.product((1, -1), repeat=3):
            v = np.array((1,) + tail, dtype=np.int64)
            if np.prod(v) != parity:
                continue
            q = np.array(quad)
            inside = family[:, q[:, None], q]
            local_two_r = 2 * np.eye(4, dtype=np.int64) - np.outer(v, v)
            four_inside = local_two_r @ inside @ local_two_r
            if not np.all(np.isin(four_inside, (-4, 4))):
                continue
            prods = np.einsum("i,sij->sj", v, family[:, q, :])
            outside = np.ones(n, dtype=bool)
            outside[q] = False
            if not np.any(prods[:, outside]) and np.array_equal(four_inside, 4 * inside):
                continue
            yield quad, v


def apply_move(family, quad, v):
    q = np.array(quad)
    local_two_r = 2 * np.eye(4, dtype=np.int64) - np.outer(v, v)
    four_inside = local_two_r @ family[:, q[:, None], q] @ local_two_r
    prods = np.einsum("i,sij->sj", v, family[:, q, :])
    new = family.copy()
    new[:, q, :] -= v[None, :, None] * prods[:, None, :] // 2
    new[:, :, q] = new[:, q, :].transpose(0, 2, 1)
    new[:, q[:, None], q] = four_inside // 4
    assert np.all(np.abs(new) == 1)
    assert np.array_equal(new, new.transpose(0, 2, 1))
    return new


def boolean_eigenvectors(h, eigenvalue):
    """Exhaust the exact eigenspace, not the ambient cube."""
    n = len(h)
    basis = (sp.Matrix(h) - eigenvalue * sp.eye(n)).nullspace()
    if not basis:
        return np.empty((0, n), dtype=np.int64)
    b = sp.Matrix.hstack(*basis)
    d = b.shape[1]
    assert d <= 16
    pivots = list(b.T.rref()[1])
    reconstruction = b * b[pivots, :].inv()
    denominator = sp.ilcm(*[entry.q for entry in reconstruction])
    numerator = np.array(reconstruction * denominator, dtype=np.int64)
    signs = sign_vectors(d)
    candidates = signs @ numerator.T
    valid = candidates[np.all(np.abs(candidates) == int(denominator), axis=1)] // int(denominator)
    assert np.all(valid @ h == eigenvalue * valid)
    return valid


def common_quotients(family, rotation_numerator, denominator, seed_n, outer_s):
    root = math.isqrt(outer_s)
    assert root * root == outer_s
    # family[0] is the all-one full seed, so its nonzero eigenvalues are ±n√s.
    seed_spins = sign_vectors(seed_n)
    records = []
    for orientation in (1, -1):
        ys = boolean_eigenvectors(family[0], orientation * seed_n * root)
        for y in ys:
            pulled = rotation_numerator.T @ y
            # The pulled eigenspace is 1_seed ⊗ E_{±√s}(H_s).
            shaped = pulled.reshape(seed_n, outer_s)
            assert np.all(shaped == shaped[0])
            v_num = shaped[0]
            assert np.array_equal(walsh(outer_s) @ v_num, orientation * root * v_num)
            pure_inputs = np.array([np.kron(u, v_num) for u in seed_spins])
            outputs = pure_inputs @ rotation_numerator.T
            if np.all(np.abs(outputs) == denominator * denominator):
                records.append({"orientation": orientation,
                                "outer_numerator": v_num.tolist(), "outer_denominator": denominator})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-order", type=int, default=3, choices=[3, 4])
    parser.add_argument("--outer-order", type=int, default=4, choices=[4, 16])
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=9060900)
    parser.add_argument("--output")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    family = all_seed_parents(args.seed_order, args.outer_order)
    original = family.copy()
    dim = family.shape[-1]
    r_num = np.eye(dim, dtype=np.int64)
    denominator = 1
    records = []
    for step in range(args.steps):
        moves = list(universal_moves(family))
        if not moves:
            break
        quad, v = moves[int(rng.integers(len(moves)))]
        family = apply_move(family, quad, v)
        q = np.array(quad)
        two_r = 2 * np.eye(dim, dtype=np.int64)
        two_r[np.ix_(q, q)] -= np.outer(v, v)
        r_num = two_r @ r_num
        denominator *= 2
        # Cancel common factors, preventing avoidable integer growth.
        while np.all(r_num % 2 == 0):
            r_num //= 2
            denominator //= 2
        assert np.array_equal(r_num @ r_num.T, denominator**2 * np.eye(dim, dtype=np.int64))
        assert np.array_equal(r_num @ original @ r_num.T, denominator**2 * family)
        quotients = common_quotients(family, r_num, denominator, args.seed_order, args.outer_order)
        record = {"step": step, "universal_moves": len(moves), "quad": list(quad), "v": v.tolist(),
                  "rotation_denominator": denominator,
                  "maximum_row_support": int(np.count_nonzero(r_num, axis=1).max()),
                  "common_quotients": quotients}
        records.append(record)
        print(json.dumps({"step": step, "moves": len(moves), "row_support": record["maximum_row_support"],
                          "common_quotients": len(quotients)}), flush=True)
        if not quotients:
            break
    result = {"parameters": vars(args), "path": records, "rotation_numerator": r_num.tolist(),
              "rotation_denominator": denominator,
              "scope": "Universal finite sign-preserving family and exact common-quotient audit; no asymptotic cap assertion."}
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")


if __name__ == "__main__":
    main()
