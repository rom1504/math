"""Exact real Kerdock MUB construction; numerical seed-lift witnesses only."""
import itertools
import json
from pathlib import Path
import numpy as np
from flatify_construct_2026_09_07_mub_diagnostic import best_witness


def mul8(a, b):
    z = 0
    while b:
        if b & 1:
            z ^= a
        a <<= 1
        if a & 8:
            a ^= 11
        b >>= 1
    return z


def trace8(a):
    z, v = 0, a
    for _ in range(3):
        z ^= v
        v = mul8(v, v)
    assert z in (0, 1)
    return z


def polar(a, u, v):
    x, s = u & 7, u >> 3
    y, r = v & 7, v >> 3
    ax, ay = trace8(mul8(a, x)), trace8(mul8(a, y))
    return (trace8(mul8(mul8(a, a), mul8(x, y)))
            ^ (ax & ay) ^ (s & ay) ^ (r & ax))


def bases16():
    n = 16
    h = np.array([[(-1) ** bin(a & b).count("1") for b in range(n)]
                  for a in range(n)], dtype=np.int64)
    bases = [4 * np.eye(n, dtype=np.int64)]
    for a in range(8):
        form = np.array([[polar(a, 1 << i, 1 << j) for j in range(4)]
                         for i in range(4)], dtype=np.int64)
        assert np.all(np.diag(form) == 0)
        phases = []
        for v in range(16):
            q = sum(int(form[i, j]) for i in range(4) for j in range(i+1, 4)
                    if (v >> i) & 1 and (v >> j) & 1) % 2
            phases.append((-1) ** q)
        bases.append(h * phases)
    for i, u in enumerate(bases):
        assert np.array_equal(u @ u.T, 16 * np.eye(n, dtype=np.int64))
        for v in bases[:i]:
            assert np.all(np.abs(u @ v.T) == 4)
    return bases


def main():
    bases = bases16()
    a = np.array(json.loads(Path("computations/results/exact_m5.json").read_text())["matrix"], dtype=np.int64)
    rng = np.random.default_rng(202609070610)
    results = []
    for ids in itertools.combinations(range(9), 5):
        c = np.block([[np.zeros((16, 16), dtype=np.int64) if i == j
                       else a[i, j] * (bases[ids[i]] @ bases[ids[j]].T // 4)
                       for j in range(5)] for i in range(5)])
        value, witness = best_witness(c, rng, trials=40)
        results.append(dict(bases=list(ids), witness_cap=value,
                            witness=witness.tolist()))
    results.sort(key=lambda r: r["witness_cap"])
    payload = dict(dimension=16, exact_mub_count=9,
                          seed_cap=4, naive_transfer_target=256,
                          minimum_found_witness=results[0]["witness_cap"],
                          maximum_found_witness=results[-1]["witness_cap"],
                          results=results)
    output = Path("computations/results/flatify_construct_2026_09_07_kerdock_diagnostic.json")
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
