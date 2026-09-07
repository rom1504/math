"""Bounded actual-cap diagnostic for correlated skew Gaussian orientations.

Every reported cap is exhaustive at the chosen finite order. No scaling
claim or favorable Gaussian parent bound is inferred from the experiment.
"""
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


def queries(a):
    n = len(a)
    edges = list(itertools.combinations(range(n), 2))
    offsets, coefficients = [], []
    for subset in range(1, (1 << n) - 1, 2):
        inside = [i for i in range(n) if subset & (1 << i)]
        outside = [i for i in range(n) if not subset & (1 << i)]
        free = [i for i in range(n) if i not in (inside[0], outside[0])]
        for choices in itertools.product((-1, 1), repeat=len(free)):
            spin = np.ones(n, dtype=int)
            spin[free] = choices
            u, v = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
            u[inside], v[outside] = spin[inside], spin[outside]
            offsets.append(abs(u @ a @ v))
            coefficients.append([u[i] * v[j] - u[j] * v[i] for i, j in edges])
    return np.array(offsets), np.array(coefficients, dtype=float), edges


def one_seed(a, samples, rng):
    n = len(a)
    offsets, tests, edges = queries(a)
    row, col = np.array(edges).T
    upper = rng.normal(size=(samples, len(edges)))
    g = np.zeros((samples, n, n))
    g[:, row, col], g[:, col, row] = upper, -upper
    anti = a @ g + g @ a
    # A symmetric Gaussian driver is needed for a skew commutator.
    symmetric = np.zeros_like(g)
    symmetric[:, row, col] = upper
    symmetric[:, col, row] = upper
    diag = rng.normal(scale=np.sqrt(2), size=(samples, n))
    symmetric[:, np.arange(n), np.arange(n)] = diag
    comm = a @ symmetric - symmetric @ a
    modes = {"iid_skew": g, "anticommutator": anti, "commutator_control": comm}
    result = {}
    for name, matrices in modes.items():
        signs = np.sign(matrices[:, row, col])
        assert not np.any(signs == 0)
        caps = 2 * np.max(offsets[:, None] + np.abs(tests @ signs.T), axis=0)
        result[name] = {
            "samples": samples,
            "minimum_cap": int(caps.min()),
            "mean_cap": float(caps.mean()),
            "histogram": dict(sorted(Counter(map(int, caps)).items())),
        }
    # The exact variance formula is checked on a complete Gaussian basis.
    variance_checks = 0
    for unused in range(30):
        bits = rng.integers(0, 2, n)
        spin = rng.choice((-1, 1), n)
        u, v = spin * bits, spin * (1 - bits)
        actual = 0
        for i, j in edges:
            basis = np.zeros((n, n), dtype=int)
            basis[i, j], basis[j, i] = 1, -1
            value = u @ (a @ basis + basis @ a) @ v
            actual += int(value) ** 2
        predicted = int(v @ v) * int((a @ u) @ (a @ u))
        predicted += int(u @ u) * int((a @ v) @ (a @ v))
        predicted -= 2 * int(u @ a @ v) ** 2
        predicted += 2 * int(u @ a @ u) * int(v @ a @ v)
        assert actual == predicted
        variance_checks += 1
    for i, j in edges:
        actual = 0
        for k, ell in edges:
            basis = np.zeros((n, n), dtype=int)
            basis[k, ell], basis[ell, k] = 1, -1
            actual += int((a @ basis + basis @ a)[i, j]) ** 2
        assert actual == 2 * (n - 2)
    spins = np.array([(1,) + tail for tail in itertools.product((-1, 1), repeat=n - 1)])
    child = int(np.max(np.abs(np.sum(spins * (spins @ a), axis=1))) / 2)
    return {"n": n, "child_cap": child, "cut_tests": len(tests),
            "variance_checks": variance_checks, "laws": result}


if __name__ == "__main__":
    rng = np.random.default_rng(202609072011)
    results = []
    for n in (6, 7, 8):
        path = Path(__file__).parent / "results" / ("exact_m%d.json" % n)
        matrix = np.asarray(json.loads(path.read_text())["matrix"], dtype=int)
        results.append(one_seed(matrix, 500, rng))
    print(json.dumps({"status": "BOUNDED_DIAGNOSTIC", "results": results,
                      "qualification": "Finite exact cap tests; no asymptotic claim"}, indent=2))
