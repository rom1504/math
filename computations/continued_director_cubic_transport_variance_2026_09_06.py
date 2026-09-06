"""Exact seed tests for cubic transported-noise variance.

For a symmetric full sign seed S of order k, D=S/sqrt(k), C=D^2,
and T=D(C entrywise-cubed)D.  Each diagonal of T has denominator k^4.
Such seeds give actual hollow-signing sequences after tensoring a growing
symmetric Hadamard and deleting the diagonal.  This tests pointwise
variance >=1, not a claim about actual minimizing signings.
"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-exact", type=int, default=5)
    parser.add_argument("--samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=730615)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    records = []
    for k in range(2, 17):
        pairs = list(zip(*np.triu_indices(k)))
        if k <= args.max_exact:
            words = itertools.product((-1, 1), repeat=len(pairs))
            count = 1 << len(pairs)
            kind = "exhaustive"
        else:
            words = rng.choice((-1, 1), size=(args.samples, len(pairs)))
            count = args.samples
            kind = "sampled"
        best = None
        best_s = None
        best_avg_sqrt = float("inf")
        for word in words:
            s = np.empty((k, k), dtype=np.int64)
            for (i, j), value in zip(pairs, word):
                s[i, j] = s[j, i] = value
            c = s @ s
            tnum = s @ (c ** 3) @ s
            assert np.trace(tnum) == np.sum(c ** 4)
            assert np.trace(tnum) >= k ** 5
            diagonal = np.diag(tnum)
            value = int(np.min(diagonal))
            if best is None or value < best:
                best, best_s = value, s.tolist()
            best_avg_sqrt = min(best_avg_sqrt, float(np.mean(np.sqrt(diagonal / k**4))))
        item = {"k": k, "kind": kind, "count": count,
                "minimum_diagonal_numerator": best, "denominator": k**4,
                "pointwise_ge_one_survives": best >= k**4,
                "witness": best_s,
                "minimum_average_sqrt_exploratory": best_avg_sqrt}
        records.append(item)
        print(json.dumps(item), flush=True)
    result = {"purpose": "cubic transport pointwise variance falsifier",
              "seed": args.seed, "records": records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
