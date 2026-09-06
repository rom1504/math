"""Exact rational finite replay of isotropic-response extension bounds."""
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("exact_isotropy_helpers", ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)


def normalized_key(matrix):
    a = np.asarray(matrix, dtype=np.int16)
    n = len(a)
    signs = np.ones(n, dtype=np.int16)
    if n > 1:
        signs[1:] = a[0, 1:]
    b = signs[:, None] * a * signs[None, :]
    return n, tuple(int(b[i, j]) for i in range(n) for j in range(i+1, n))


def exact_response(key):
    n, coefficients = key
    edges = list(itertools.combinations(range(n), 2))
    spins = np.asarray([(1,) + tail for tail in itertools.product((-1, 1), repeat=n-1)], dtype=np.int16)
    features = np.asarray([[x[i] * x[j] for i, j in edges] for x in spins], dtype=np.int16)
    values = abs(features @ np.asarray(coefficients, dtype=np.int16))
    eq = np.vstack((np.ones(len(spins)), features.T))
    rhs = np.zeros(1 + len(edges)); rhs[0] = 1
    lp = linprog(-values.astype(float), A_eq=eq, b_eq=rhs, bounds=(0, None), method="highs")
    assert lp.success
    law = helpers.exact_weights(lp, features)
    assert law is not None
    response = sum(weight * int(values[index]) for index, weight in law)
    dual = [F(float(z)).limit_denominator(1000000) for z in lp.eqlin.marginals]
    valid = lambda y: y is not None and y[0] == -response and all(
        y[0] + sum(z * int(f) for z, f in zip(y[1:], row)) <= -int(value)
        for row, value in zip(features, values))
    if not valid(dual):
        tight = np.flatnonzero(abs(eq.T @ lp.eqlin.marginals + values) < 1e-7)
        dual = helpers.exact_linear_solution(
            [[1] + [int(z) for z in features[i]] for i in tight],
            [-int(values[i]) for i in tight], lp.eqlin.marginals)
    assert valid(dual)
    assert response * response <= F(n*n*(n-1), 4)
    return response, {
        "order": n, "root_normalized_edge_coefficients": coefficients,
        "response": str(response),
        "primal_law": [{"spin": spins[index].tolist(), "weight": str(weight)}
                       for index, weight in law],
        "dual_upper_constant": str(-dual[0]),
        "dual_upper_edge_coefficients": [str(-z) for z in dual[1:]],
        "certificate": "exact nonnegative isotropic primal law and pointwise quadratic majorant agree",
    }


def below_sum_of_radicals(t, a, b):
    """Exact comparison t <= sqrt(a)+sqrt(b) for nonnegative rationals."""
    assert t >= 0 and a >= 0 and b >= 0
    if t*t <= a+b:
        return True
    return (t*t-a-b)**2 <= 4*a*b


def main():
    cache = {}
    certificates = []
    def response(a):
        key = normalized_key(a)
        if key not in cache:
            value, cert = exact_response(key)
            cache[key] = value
            certificates.append(cert)
        return cache[key]

    comparisons = 0
    order_counts = {}
    for total in range(2, 6):
        residual = list(itertools.combinations(range(1, total), 2))
        for edge_signs in itertools.product((-1, 1), repeat=len(residual)):
            b = np.ones((total, total), dtype=np.int16) - np.eye(total, dtype=np.int16)
            for (i, j), sign in zip(residual, edge_signs):
                b[i, j] = b[j, i] = sign
            large = response(b)
            for old in range(1, total):
                for subset in itertools.combinations(range(total), old):
                    small = response(b[np.ix_(subset, subset)])
                    increment = large - small
                    assert increment >= 0
                    added = total-old
                    assert below_sum_of_radicals(increment / added, F(old), F(added-1, 4))
                    comparisons += 1
                    order_counts[total] = order_counts.get(total, 0) + 1
    result = {
        "scope": "All root-normalized full signing extensions up to total order five, every nonempty proper old subset",
        "bound": "0 <= I(B)-I(A) <= r sqrt(n) + r sqrt(r-1)/2",
        "exact_extension_comparisons": comparisons,
        "comparison_counts_by_total_order": order_counts,
        "distinct_exact_response_LPs": len(certificates),
        "certificates": certificates,
    }
    output = ROOT / "computations/results/transfer_adversary_isotropic_response_extension_2026_09_06.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "certificates"}))
    print("PASS: exact rational primal/dual response values; exact radical comparisons")


if __name__ == "__main__":
    main()
