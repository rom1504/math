"""Exact MUB block construction and seeded numerical Boolean-cap search.

Any returned witness is rigorously evaluable using integer arithmetic;
failure to exceed the seed cap is only a failed numerical search.
"""
import json
from pathlib import Path
import numpy as np


def mul4(a, b):
    z = 0
    while b:
        if b & 1:
            z ^= a
        a <<= 1
        if a & 4:
            a ^= 7
        b >>= 1
    return z


def bases4():
    q = 4
    h = np.array([[(-1) ** bin(a & b).count("1") for b in range(q)]
                  for a in range(q)], dtype=np.int64)
    bases = []
    for slope in [None, 0, 1, 2, 3]:
        rows = []
        for intercept in range(q):
            points = ([(intercept, z) for z in range(q)] if slope is None
                      else [(z, mul4(slope, z) ^ intercept) for z in range(q)])
            for character in range(q):
                row = np.zeros(q * q, dtype=np.int64)
                for z, (a, b) in enumerate(points):
                    row[q * a + b] = h[character, z]
                rows.append(row)
        bases.append(np.array(rows, dtype=np.int64))
    for i, u in enumerate(bases):
        assert np.array_equal(u @ u.T, q * np.eye(q * q, dtype=np.int64))
        for v in bases[:i]:
            assert np.all(np.abs(u @ v.T) == 1)
    return bases


def best_witness(c, rng, trials=300):
    n = len(c)
    best = -1
    bestx = None
    for polarity in [1, -1]:
        d = polarity * c
        for _ in range(trials):
            x = rng.choice([-1, 1], n)
            fields = d @ x
            for step in range(20 * n):
                gains = -2 * x * fields
                gain = gains.max()
                if gain <= 0:
                    break
                j = int(rng.choice(np.flatnonzero(gains == gain)))
                x[j] *= -1
                fields += 2 * x[j] * d[:, j]
            energy = int(x @ c @ x // 2)
            if abs(energy) > best:
                best, bestx = abs(energy), x.copy()
    return best, bestx


def main():
    bases = bases4()
    rng = np.random.default_rng(202609070541)
    outputs = []
    for k in [4, 5]:
        source = f"computations/results/exact_m{k}.json"
        a = np.array(json.loads(Path(source).read_text())["matrix"], dtype=np.int64)
        c = np.zeros((16 * k, 16 * k), dtype=np.int64)
        for i in range(k):
            for j in range(k):
                if i != j:
                    c[16*i:16*(i+1), 16*j:16*(j+1)] = a[i, j] * bases[i] @ bases[j].T
        value, x = best_witness(c, rng)
        seed_q = 4
        outputs.append(dict(k=k, fibre_dimension=16, seed_source=source,
                            zero_within_fibre=True, seed_transfer_target=64*seed_q,
                            numerical_witness_cap=value, witness=x.tolist(),
                            exceeds_exact_transfer=value > 64*seed_q))
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
