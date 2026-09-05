"""Independent exact small-order tests of the signed Sylvester lift.

Run with .venv/bin/python -B computations/fresh_limit_algebra_lift.py --n 6.
"""
import argparse
import itertools
import json
import numpy as np


def spins(n):
    # Global reversal leaves every quadratic energy unchanged.
    masks = np.arange(1 << (n - 1), dtype=np.uint64)
    tail = 1 - 2 * ((masks[:, None] >> np.arange(n - 1, dtype=np.uint64)) & 1)
    return np.concatenate((np.ones((len(masks), 1), dtype=np.int16), tail.astype(np.int16)), axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=6)
    args = parser.parse_args()
    n = args.n
    X = spins(n)
    Y = spins(2 * n)
    upper = list(itertools.combinations(range(n), 2))
    free = list(itertools.combinations(range(1, n), 2))
    F = np.array([X[:, i] * X[:, j] for i, j in upper]).T
    best = 10 ** 10
    records = []
    for code in range(1 << len(free)):
        A = np.ones((n, n), dtype=np.int16) - np.eye(n, dtype=np.int16)
        for bit, (i, j) in enumerate(free):
            A[i, j] = A[j, i] = 1 - 2 * ((code >> bit) & 1)
        coeff = np.array([A[i, j] for i, j in upper])
        vals = F @ coeff
        q = int(np.max(np.abs(vals)))
        if q < best:
            best = q
            records = []
        if q == best:
            records.append((code, A.copy()))
    lift_values = []
    for code, A in records:
        L = np.block([[A, A + np.eye(n, dtype=np.int16)], [A + np.eye(n, dtype=np.int16), -A]])
        qlift = int(np.max(np.abs(np.sum((Y @ L) * Y, axis=1) // 2)))
        lift_values.append((qlift, code))
    lift_values.sort()
    print(json.dumps({'n': n, 'M_n': best, 'number_switch_normalized_minimizers': len(records),
                      'min_lift_Q': lift_values[0][0], 'max_lift_Q': lift_values[-1][0],
                      'target': (2 ** 1.5) * best,
                      'best_lift_codes': lift_values[:5], 'worst_lift_codes': lift_values[-5:]}))


if __name__ == '__main__':
    main()
