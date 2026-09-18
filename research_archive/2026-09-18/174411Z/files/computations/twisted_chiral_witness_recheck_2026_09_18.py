#!/usr/bin/env python3
"""Independent all-spin evaluation of submitted twisted-double witnesses.

Enumerates x_0=+1 and every y. Integer energies are evaluated as float64
matrix products of +/-1 arrays; all sums remain exactly representable.
This independently avoids the cut-profile implementation used to search.
"""

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def spins(n, projective=False):
    k = n - int(projective)
    ids = np.arange(1 << k, dtype=np.uint64)
    x = 2 * ((ids[:, None] >> np.arange(k, dtype=np.uint64)) & 1).astype(np.int64) - 1
    return np.column_stack((np.ones(len(x), dtype=np.int64), x)) if projective else x


def exhaustive_double(a, b, d):
    n = len(a)
    xs = spins(n, True)
    ys = spins(n)
    hx = np.einsum("bi,ij,bj->b", xs, a, xs) // 2
    hy = np.einsum("bi,ij,bj->b", ys, a, ys) // 2
    c = b + np.diag(d)
    yfloat = ys.astype(float).T
    best = -1
    witness = None
    histogram = np.zeros(4 * n * n + 1, dtype=np.int64)
    shift = 2 * n * n
    for start in range(0, len(xs), 256):
        xx = xs[start:start + 256]
        values = (xx @ c).astype(float) @ yfloat + hx[start:start + len(xx), None] - hy[None, :]
        assert np.all(values == np.rint(values))
        integer = values.astype(np.int64)
        histogram += np.bincount(integer.ravel() + shift, minlength=len(histogram))
        index = int(np.argmax(abs(integer)))
        cap = int(abs(integer.flat[index]))
        if cap > best:
            i, j = np.unravel_index(index, integer.shape)
            best = cap
            witness = {"x": xx[i].tolist(), "y": ys[j].tolist(),
                       "energy": int(integer[i, j]), "H_A_x": int(hx[start+i]),
                       "H_A_y": int(hy[j]), "cross": int(xx[i] @ c @ ys[j])}
    assert int(histogram.sum()) == (1 << (2*n-1))
    assert np.array_equal(histogram, histogram[::-1])
    return {"cap": best, "maximizer": witness, "projective_states": int(histogram.sum()),
            "energy_histogram": {str(i - shift): int(count) for i, count in enumerate(histogram) if count}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-n", type=int, default=12)
    args = parser.parse_args()
    reports = []
    seen = set()
    for source in args.inputs:
        data = json.loads(source.read_text())
        for index, record in enumerate(data["records"]):
            a = np.asarray(record["child_matrix"], dtype=np.int64)
            if len(a) > args.max_n:
                continue
            b = np.asarray(record["bridge_hollow_matrix"], dtype=np.int64)
            p = np.asarray(record["p"], dtype=np.int64)
            s = np.asarray(record["s"], dtype=np.int64)
            d = np.asarray(record["d"], dtype=np.int64)
            key = (a.tobytes(), b.tobytes(), d.tobytes())
            if key in seen:
                continue
            seen.add(key)
            assert sorted(p.tolist()) == list(range(len(a)))
            assert np.all(abs(s) == 1) and np.all(abs(d) == 1)
            assert np.array_equal(b, a[np.ix_(p, p)] * s[:, None] * s[None, :])
            parent = np.block([[a, b + np.diag(d)], [b + np.diag(d), -a]])
            assert np.array_equal(parent, np.asarray(record["parent_matrix"]))
            report = exhaustive_double(a, b, d)
            assert report["cap"] == record["cap"]
            report.update({"source": str(source), "record_index": index, "seed_order": len(a),
                           "claimed_cap": record["cap"], "reconstruction_verified": True})
            reports.append(report)
            print(json.dumps({k: report[k] for k in ["seed_order", "cap", "projective_states"]}), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"method": "all x0=+1, all y; direct block energy", "records": reports}, indent=2) + "\n")


if __name__ == "__main__":
    main()
