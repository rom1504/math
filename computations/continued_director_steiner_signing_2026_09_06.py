"""Exact asymmetric-spectrum hollow signings from binary Steiner triples.

This is a deterministic stress-test family, not a claim of minimizing Q.
The checked identity is S^2=(r+3)S+(N-1)I, with r=2^(m-1)-1.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def make_signing(m):
    v = 2**m - 1
    r = (v - 1) // 2
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < r + 1:
        h = np.block([[h, h], [h, -h]])
    triples = sorted({tuple(sorted((a, b, a ^ b)))
                      for a in range(1, v + 1) for b in range(a + 1, v + 1)})
    assert len(triples) == v * r // 3
    incidence = [[] for _ in range(v)]
    pair_count = {}
    for row, triple in enumerate(triples):
        for point in triple:
            incidence[point - 1].append(row)
        for i in range(3):
            for j in range(i + 1, 3):
                pair = (triple[i], triple[j])
                pair_count[pair] = pair_count.get(pair, 0) + 1
    assert len(pair_count) == v * (v - 1) // 2
    assert set(pair_count.values()) == {1}
    n = v * (r + 1)
    e = np.zeros((len(triples), n), dtype=np.int64)
    for point, rows in enumerate(incidence):
        assert len(rows) == r
        e[np.array(rows), point * (r + 1):(point + 1) * (r + 1)] = h[1:]
    assert np.array_equal(e @ e.T, 3 * (r + 1) * np.eye(len(triples), dtype=np.int64))
    s = e.T @ e - r * np.eye(n, dtype=np.int64)
    assert np.array_equal(s, s.T)
    assert np.all(np.diag(s) == 0)
    assert np.all(np.abs(s + np.eye(n, dtype=np.int64)) == 1)
    assert np.array_equal(s @ s, (r + 3) * s + (n - 1) * np.eye(n, dtype=np.int64))
    return s, {"m": m, "v": v, "r": r, "N": n, "positive_eigenvalue": 2*r+3,
               "positive_multiplicity": len(triples), "negative_eigenvalue": -r,
               "negative_multiplicity": n-len(triples),
               "normalized_cubic_trace_squared": {"numerator": (r+3)**2, "denominator": n-1}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for m in range(2, args.max_m + 1):
        _, record = make_signing(m)
        records.append(record)
        print("PASS", json.dumps(record), flush=True)
    if args.output:
        args.output.write_text(json.dumps({"exact": True, "records": records}, indent=2) + "\n")


if __name__ == "__main__":
    main()
