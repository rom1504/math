"""Exact symmetric quadruple-reflection moves and small Boolean-cap tests.

R=I-vv^T/2, with v supported on four coordinates. Valid moves return
R H R as a full symmetric sign matrix, preserving its entire spectrum.
This is not ordinary diagonal Seidel switching.
"""

import argparse
import itertools
import json

import numpy as np


def walsh(order):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < order:
        h = np.block([[h, h], [h, -h]])
    assert len(h) == order
    return h


def cap(h):
    n = len(h)
    assert n <= 22
    result = 0
    low, high = 10**20, -(10**20)
    for offset in range(0, 1 << (n - 1), 1 << 14):
        codes = np.arange(offset, min(offset + (1 << 14), 1 << (n - 1)), dtype=np.uint64)
        x = np.ones((len(codes), n), dtype=np.int64)
        x[:, 1:] = 1 - 2 * ((codes[:, None] >> np.arange(n - 1, dtype=np.uint64)) & 1).astype(np.int64)
        values = np.sum((x @ h) * x, axis=1)
        low = min(low, int(values.min()))
        high = max(high, int(values.max()))
        result = max(result, int(np.abs(values).max()))
    return {"full_twice_cap": result, "twice_range": [low, high], "half_cap": result / 2}


def quadruple_move(h, quad, v):
    q = np.asarray(quad, dtype=np.int64)
    v = np.asarray(v, dtype=np.int64)
    products = v @ h[q, :]
    if not np.all(np.isin(products, (-4, 0, 4))):
        return None
    inside = h[np.ix_(q, q)]
    vv = np.outer(v, v)
    four_new_inside = 4 * inside - 2 * vv @ inside - 2 * inside @ vv + vv @ inside @ vv
    if not np.all(np.isin(four_new_inside, (-4, 4))):
        return None
    new = h.copy()
    new[q, :] -= np.outer(v, products) // 2
    new[:, q] = new[q, :].T
    new[np.ix_(q, q)] = four_new_inside // 4
    assert np.array_equal(new, new.T)
    assert np.all(np.abs(new) == 1)
    # For this finite verifier, check the full similarity identity exactly.
    two_r = 2 * np.eye(len(h), dtype=np.int64)
    two_r[np.ix_(q, q)] -= vv
    assert np.array_equal(two_r @ h @ two_r, 4 * new)
    assert np.trace(new) == np.trace(h)
    return new


def all_moves(h):
    for quad in itertools.combinations(range(len(h)), 4):
        q = np.asarray(quad)
        prod = np.prod(h[q, :], axis=0)
        if not np.all(prod == prod[0]):
            continue
        patterns = h[q, :].T.copy()
        patterns *= patterns[:, 0, None]
        patterns = np.unique(patterns, axis=0)
        for v in patterns:
            new = quadruple_move(h, quad, v)
            if new is not None and not np.array_equal(new, h):
                yield quad, v, new


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=16)
    parser.add_argument("--steps", type=int, default=12)
    parser.add_argument("--seed", type=int, default=9060812)
    parser.add_argument("--output")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    h = walsh(args.order)
    start = cap(h) if args.order <= 22 else None
    records = []
    for step in range(args.steps):
        moves = list(all_moves(h))
        if not moves:
            records.append({"step": step, "valid_nontrivial_moves": 0})
            break
        quad, v, new = moves[int(rng.integers(len(moves)))]
        h = new
        assert np.array_equal(h @ h, args.order * np.eye(args.order, dtype=np.int64))
        records.append({"step": step, "valid_nontrivial_moves": len(moves), "quad": list(quad),
                        "v": v.tolist(), "cap": cap(h) if args.order <= 22 else None})
    result = {"parameters": vars(args), "initial": start, "path": records, "final_matrix": h.tolist(),
              "scope": "Exact finite symmetric flat-sign orbit; no asymptotic transfer claim."}
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            stream.write(rendered + "\n")


if __name__ == "__main__":
    main()
