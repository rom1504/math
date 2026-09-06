"""Exact-integer finite checks for the rank-one weave transfer audit."""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def sylvester(m):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < m:
        h = np.block([[h, h], [h, -h]])
    assert len(h) == m
    return h


def cap(a):
    n = len(a)
    assert n <= 20
    values = []
    for tail in itertools.product((-1, 1), repeat=n - 1):
        x = np.array((1,) + tail, dtype=np.int64)
        values.append(int(x @ a @ x) // 2)
    return max(map(abs, values))


def weave(s, hs, selectors):
    m = len(s)
    return np.block([
        [s[i, j] * np.outer(hs[i][selectors[i], j], hs[j][selectors[j], i])
         for j in range(m)]
        for i in range(m)
    ])


def main():
    source = ROOT / "computations/results/exact_m8.json"
    a = np.array(json.loads(source.read_text())["matrix"], dtype=np.int64)
    m = len(a)
    assert cap(a) == 10
    diagonal = np.array([1] * (m // 2) + [-1] * (m // 2), dtype=np.int64)
    s = a + np.diag(diagonal)
    rng = np.random.default_rng(20260906)
    h = sylvester(m)
    hs = [h[:, rng.permutation(m)] * rng.choice((-1, 1), size=m)
          for _ in range(m)]
    full_t = [np.arange(m) for _ in range(m)]
    full = weave(s, hs, full_t)
    assert np.array_equal(full @ full, m * m * np.eye(m * m, dtype=np.int64))
    assert int(np.trace(full)) == 0
    pairs = [(i, i + 1) for i in range(0, m, 2)]
    pi = dict(pairs + [(j, i) for i, j in pairs])

    for sigma in (-1, 1):
        u = np.ones(m, dtype=np.int64)
        for i, j in pairs:
            u[j] = sigma * s[i, j]
        x = np.concatenate([u[i] * hs[i][:, pi[i]] for i in range(m)])
        assert np.array_equal(full @ x, sigma * m * x)

    records = []
    for k in (3, 5, 7, 8):
        ts = [np.sort(rng.choice(m, k, replace=False)) for _ in range(m)]
        restricted = weave(s, hs, ts)
        restricted -= np.diag(np.diag(restricted))
        means = []
        max_witness = 0
        for sigma in (-1, 1):
            total = 0
            for phases in itertools.product((-1, 1), repeat=m // 2):
                u = np.ones(m, dtype=np.int64)
                for z, (i, j) in zip(phases, pairs):
                    u[i] = z
                    u[j] = sigma * s[i, j] * z
                x = np.concatenate([u[i] * hs[i][ts[i], pi[i]] for i in range(m)])
                value = int(x @ restricted @ x) // 2
                total += value
                max_witness = max(max_witness, abs(value))
            assert total % (1 << (m // 2)) == 0
            means.append(total // (1 << (m // 2)))
        assert means[1] - means[0] == m * k * k
        floor = m * k * k // 2
        assert max(map(abs, means)) >= floor
        records.append({"k": k, "conditional_means": means,
                        "exact_floor": floor, "found_witness_cap": max_witness})

    alternate = s.copy()
    ds = np.ones((m, m), dtype=np.int64)
    for i in range(m):
        for j in range(i + 1, m):
            alternate[i, j] = alternate[j, i] = int(rng.choice((-1, 1)))
            ds[i, j] = s[i, j] * alternate[i, j]
    absorbed_hs = [hs[i] * ds[i] for i in range(m)]
    assert np.array_equal(full, weave(alternate, absorbed_hs, full_t))

    square = np.kron(s, s)
    hollow_square = square - np.diag(np.diag(square))
    x = s.ravel()
    witness = int(x @ hollow_square @ x) // 2
    fourth = int(np.trace(s @ s @ s @ s))
    assert 2 * witness == fourth - int(np.trace(s)) ** 2
    assert witness >= (m ** 3 - m ** 2) // 2

    result = {
        "seed_order": m, "enumerated_seed_cap": cap(a),
        "seed_normalized_cap": cap(a) / m ** 1.5,
        "full_weave_exact_hollow_cap": m ** 3 // 2,
        "full_weave_normalized_cap": 0.5,
        "restricted_checks": records,
        "seed_absorption_entrywise_verified": True,
        "self_square_witness_cap": witness,
        "self_square_normalized_witness": witness / m ** 3,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
