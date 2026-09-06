"""Exact small Walsh spectrum-type counts; no asymptotic inference."""

import argparse
import collections
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import hadamard


def record(m):
    h = hadamard(m).astype(np.int16)
    counts = collections.Counter()
    count = 1 << m
    for start in range(0, count, 4096):
        ids = np.arange(start, min(start + 4096, count), dtype=np.uint64)
        x = (1 - 2 * ((ids[:, None] >> np.arange(m, dtype=np.uint64)) & 1).astype(np.int16))
        z = np.abs(x @ h)
        for row in z:
            counts[tuple(np.bincount(row // 2, minlength=m // 2 + 1).tolist())] += 1
    rows = []
    for key, size in counts.items():
        orbit = math.factorial(m) // math.prod(math.factorial(c) for c in key)
        support = m - key[0]
        score = math.log2(size) - (math.log2(orbit) + support) / 2
        rows.append({"absolute_value_counts": {str(2 * j): c for j, c in enumerate(key) if c},
                     "boolean_row_count": size, "magnitude_orbit_size": orbit,
                     "support": support, "log2_count_minus_half_type_entropy": score})
    rows.sort(key=lambda r: r["log2_count_minus_half_type_entropy"], reverse=True)
    assert sum(r["boolean_row_count"] for r in rows) == count
    return {"m": m, "enumerated_rows": count, "distinct_types": len(rows), "types": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", type=int, default=[4, 8, 16])
    parser.add_argument("--output", default="computations/results/continued_director_walsh_type_count_2026_09_06.json")
    args = parser.parse_args()
    result = {"status": "exact complete finite enumeration; logarithms floating; no asymptotic claim",
              "records": [record(m) for m in args.orders]}
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summaries": [
        {"m": r["m"], "types": r["distinct_types"], "largest_scores": r["types"][:4]}
        for r in result["records"]]}, indent=2))
