"""Exact pentagon witness for quadratic-phase MUBs, including Heisenberg audit."""
import itertools
import json
from pathlib import Path
import numpy as np
from flatify_construct_2026_09_07_kerdock_diagnostic import bases16


def parity(x):
    return bin(int(x)).count("1") % 2


def quadratic(bits):
    n = int(np.log2(len(bits)))
    constant = int(bits[0])
    linear = [int(bits[1 << i]) ^ constant for i in range(n)]
    pair = {(i, j): int(bits[(1 << i) | (1 << j)]) ^ constant ^ linear[i] ^ linear[j]
            for i in range(n) for j in range(i+1, n)}
    for v in range(len(bits)):
        value = constant
        for i in range(n):
            if (v >> i) & 1:
                value ^= linear[i]
        for (i, j), c in pair.items():
            if (v >> i) & 1 and (v >> j) & 1:
                value ^= c
        assert value == int(bits[v])


def witness_for_subset(bases, ids, orbit=False):
    t = len(bases[0])
    ref = next(i for i in range(len(bases)) if i not in ids)
    relative = [bases[i] @ bases[ref].T // 4 for i in ids]
    phases = [(v[0] < 0).astype(np.int64) for v in relative]
    for f in phases:
        quadratic(f)
    def bilinear(i, u, v):
        f = phases[i]
        return int(f[0] ^ f[u] ^ f[v] ^ f[u ^ v])
    u = 1
    v = next(v for v in range(t)
             if (bilinear(0, u, v) ^ bilinear(1, u, v)) == 1
             and (bilinear(0, u, v) ^ bilinear(2, u, v)) == 1)
    c = bilinear(1, u, v)
    support = [0, u, v, u ^ v]
    patterns = np.array([[-1, 1, 1, -1], [-1, 1, -1, -1],
                         [-1, 1, -1, -1], [0, 0, 2, 0],
                         [0, 0, 2, 0]], dtype=np.int64)
    patterns[:, 3] *= (-1) ** c
    w = np.zeros((5, t), dtype=np.int64)
    w[:, support] = patterns
    xs = np.array([relative[i] @ w[i] // 2 for i in range(5)])
    assert np.all(np.abs(xs) == 1)
    a = np.array(json.loads(Path("computations/results/exact_m5.json").read_text())["matrix"], dtype=np.int64)
    energy = sum(int(a[i, j] * xs[i] @ (bases[ids[i]] @ bases[ids[j]].T // 4) @ xs[j])
                 for i in range(5) for j in range(i+1, 5))
    assert energy == 320
    if orbit:
        samples = []
        for translation in range(t):
            for character in range(t):
                transformed = np.array([[(-1) ** parity(character & z) * w[i, z ^ translation]
                                         for z in range(t)] for i in range(5)])
                assert np.array_equal(transformed @ transformed.T, w @ w.T)
                x = np.array([relative[i] @ transformed[i] // 2 for i in range(5)])
                assert np.all(np.abs(x) == 1)
                samples.append(x)
        samples = np.array(samples)
        for i in range(5):
            assert np.array_equal(samples[:, i].T @ samples[:, i],
                                  t*t*np.eye(t, dtype=np.int64))
    return dict(bases=list(ids), unused_reference=ref, u=u, v=v,
                common_quadratic_bit=c, energy=energy, witness=xs.ravel().tolist(),
                heisenberg_orbit_checked=orbit)


def main():
    bases = bases16()
    results = [witness_for_subset(bases, ids, orbit=(ids == (0, 2, 4, 6, 8)))
               for ids in itertools.combinations(range(9), 5)]
    payload = dict(exact_subsets=len(results), exact_energy=320,
                   target_from_Q_A5=256, arbitrary_fibre_completion_average=320,
                   results=results)
    output = Path("computations/results/flatify_construct_2026_09_07_kerdock_exact_witness.json")
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
