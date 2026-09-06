"""Replay odd reciprocal-variance bounds and exact even-power failure."""

import json
from fractions import Fraction

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing
from continued_feedback_variance_falsifier_2026_09_06 import sylvester


def odd_check(a, name):
    n = len(a)
    b = a / np.sqrt(n - 1)
    q = b @ b
    variances = []
    records = []
    for p in (1, 3, 5, 7):
        variance = np.diag(b @ (q**p) @ b)
        assert min(variance) > 0
        reciprocal_mean = np.mean(1 / variance)
        assert reciprocal_mean < 1 + 1e-10
        variances.append(variance)
        records.append({"p": p, "minimum": float(min(variance)),
                        "reciprocal_mean": float(reciprocal_mean),
                        "mean_standard_deviation": float(np.mean(np.sqrt(variance)))})
    mixture = np.array((0.0, 0.2, 0.3, 0.5)) @ variances
    assert np.mean(1 / mixture) < 1 + 1e-10
    return {"name": name, "order": n, "odd_powers": records,
            "mixture_reciprocal_mean": float(np.mean(1 / mixture))}


def even_check(m):
    twin = np.array([[1, -1], [-1, 1]], dtype=np.int64)
    a = np.kron(twin, sylvester(m))
    np.fill_diagonal(a, 0)
    n, d = len(a), len(a) - 1
    qnumerator = a @ a
    numerator = np.diag(a @ (qnumerator**2) @ a)
    # d^3 is the exact denominator of B (Q circ Q) B.
    claimed_numerator = 3 * d * d - 3 * d + 1
    assert np.all(numerator == claimed_numerator)
    q = qnumerator / d
    positive_twins = np.kron(np.ones((2, 2)), np.eye(m))
    three_covariance_source = np.exp(-3) * np.sinh(q)**2 * np.sinh(positive_twins)
    three_covariance_variance = np.diag(a @ three_covariance_source @ a) / d
    tau2 = (np.exp(-1) * np.sinh(1))**3
    beta = np.exp(-3) * np.sinh((d - 1) / d)**2 * np.sinh(1)
    predicted_three_covariance_variance = beta / d + tau2 - beta
    assert np.max(np.abs(three_covariance_variance - predicted_three_covariance_variance)) < 1e-12
    return a, {"hadamard_order": m, "signing_order": n,
               "every_even_variance": str(Fraction(claimed_numerator, d**3)),
               "every_even_variance_decimal": claimed_numerator / d**3,
               "three_distinct_covariance_globally_odd_variance_over_tau2": predicted_three_covariance_variance / tau2}


if __name__ == "__main__":
    rng = np.random.default_rng(9060902)
    upper = np.triu(rng.choice((-1, 1), size=(25, 25)), 1)
    examples = [(upper + upper.T, "fixed_random")]
    steiner, _ = make_signing(4)
    examples.append((steiner, "Steiner"))
    even = []
    for m in (2, 4, 8, 16, 32, 64):
        opposite, record = even_check(m)
        even.append(record)
        if m in (8, 64):
            examples.append((opposite, "opposite_twins"))
            apex = np.ones((2 * m + 1, 2 * m + 1), dtype=np.int64)
            apex[1:, 1:] = opposite
            np.fill_diagonal(apex, 0)
            examples.append((apex, "apex_opposite_twins"))
    print(json.dumps({"odd_floating_checks": [odd_check(a, name) for a, name in examples],
                      "even_exact_integer_checks": even}, indent=2))
