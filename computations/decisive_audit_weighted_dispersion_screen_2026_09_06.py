"""Exact small weighted quadratic minima under switching gauge.

Screen the unproved assertion that the sharp normalized unweighted lower
constant also controls sum(w)*sqrt(sum(w*w)). No asymptotic inference.
"""

import itertools
import json
import math
import numpy as np


def screen(n):
    edges = list(itertools.combinations(range(n), 2))
    free = [(i, j) for i, j in edges if i > 0]
    spins = np.array([(1,) + s for s in itertools.product((-1, 1), repeat=n-1)], dtype=np.int64)
    products = np.array([spins[:, i] * spins[:, j] for i, j in edges], dtype=np.int64)
    weight_patterns = [np.array((2,) * k + (1,) * (n-k), dtype=np.int64) for k in range(n+1)]
    weighted_products = [products * np.array([w[i]*w[j] for i,j in edges])[:, None] for w in weight_patterns]
    best = [10**10] * len(weight_patterns)
    codes = [None] * len(weight_patterns)
    free_edge_indices = [edges.index(e) for e in free]
    for start in range(0, 1 << len(free), 1024):
        ids = np.arange(start, min(start + 1024, 1 << len(free)), dtype=np.int64)
        signs = np.ones((len(ids), len(edges)), dtype=np.int64)
        signs[:, free_edge_indices] = 1 - 2 * ((ids[:, None] >> np.arange(len(free))) & 1)
        for k, wp in enumerate(weighted_products):
            caps = np.abs(signs @ wp).max(axis=1)
            at = int(caps.argmin())
            if int(caps[at]) < best[k]:
                best[k], codes[k] = int(caps[at]), int(ids[at])
    baseline = best[0] / n**1.5
    rows = []
    for k,w in enumerate(weight_patterns):
        ratio = best[k] / (int(w.sum()) * math.sqrt(int(w @ w)))
        rows.append({"twos": k, "weighted_cap": best[k], "gauge_code": codes[k],
                     "weighted_normalized": ratio, "unweighted_normalized": baseline,
                     "dispersion_lower_test": ratio >= baseline - 1e-14})
    return {"n": n, "gauge_signings": 1 << len(free), "cases": rows}


if __name__ == "__main__":
    print(json.dumps([screen(n) for n in range(3, 8)], indent=2))
