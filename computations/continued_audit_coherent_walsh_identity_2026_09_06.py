"""Finite exact-seed audit of the coherent degree-three Walsh decomposition.

The algebraic identity is exact; float tolerances only check its evaluation.
No stored output is overwritten.
"""

import itertools
import json

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def analytic_decomposition(q):
    gram = q @ q.T
    power2_gram = (q**2) @ (q**2).T
    power3_gram = (q**3) @ (q**3).T
    first = 3 * np.diag(gram)[:, None] * q - 2 * q**3
    third_covariance = 6 * (gram**3 - 3 * gram * power2_gram + 2 * power3_gram)
    return first, third_covariance


def exact_seed_check():
    a, _ = make_signing(2)
    n = len(a)
    q = (a @ a) / (n - 1)
    s = np.array(list(itertools.product((-1.0, 1.0), repeat=n)))
    d = (s @ q.T)**3
    first, third_covariance = analytic_decomposition(q)
    actual_first = d.T @ s / len(s)
    residual = d - s @ first.T
    actual_third = residual.T @ residual / len(s)
    assert np.max(np.abs(first - actual_first)) < 1e-10
    assert np.max(np.abs(third_covariance - actual_third)) < 1e-9
    return {"order": n, "all_sign_seeds": len(s), "first_error": float(np.max(np.abs(first - actual_first))), "third_covariance_error": float(np.max(np.abs(third_covariance - actual_third)))}


def scale_checks():
    records = []
    for exponent in [2, 3, 4, 5]:
        a, _ = make_signing(exponent)
        n = len(a)
        b = a / np.sqrt(n - 1)
        q = (a @ a) / (n - 1)
        first, third_covariance = analytic_decomposition(q)
        records.append({
            "order": n,
            "first_walsh_energy": float(np.sum(b * (first @ first.T)) / (2 * n)),
            "third_walsh_energy": float(np.sum(b * third_covariance) / (2 * n)),
            "third_walsh_mean_variance": float(np.trace(third_covariance) / n),
            "third_covariance_max_absolute_row_sum": float(np.max(np.sum(np.abs(third_covariance), axis=1))),
        })
    return records


if __name__ == "__main__":
    identity_first, identity_third = analytic_decomposition(np.eye(6))
    assert np.array_equal(identity_first, np.eye(6))
    assert np.array_equal(identity_third, np.zeros((6, 6)))
    print(json.dumps({"exact_seed_check": exact_seed_check(), "identity_cube_has_zero_third_walsh_part": True, "steiner_scaling_diagnostics": scale_checks()}, indent=2))
