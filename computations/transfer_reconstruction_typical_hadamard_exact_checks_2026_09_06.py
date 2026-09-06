"""Bounded checks of exact retained-Hadamard identities and smooth cosine tests.

Gram trace and every row-spin energy are checked in integer arithmetic.
Cosine expectations/covariances are numerical regressions, explicitly labeled.
"""

from itertools import product
import json
import math

from transfer_reconstruction_exact_checks_2026_09_06 import paley_twelve


def sylvester(order):
    h = [[1]]
    while len(h) < order:
        h = [row + row for row in h] + [row + [-value for value in row] for row in h]
    assert len(h) == order
    return h


def case(h, selector):
    m, k = len(h), len(selector)
    rows = [h[i] for i in selector]
    gram = [[sum(row[i] * row[j] for row in rows) for j in range(m)] for i in range(m)]
    assert all(gram[i][i] == k for i in range(m))
    assert sum(value * value for row in gram for value in row) == m * m * k
    spectra = []
    for x in product((-1, 1), repeat=k):
        w = [sum(row[j] * sign for row, sign in zip(rows, x)) for j in range(m)]
        assert sum(value * value for value in w) == k * m
        spectra.append(w)
    results = []
    for frequency in (0.3, 1.0, 1.7):
        test_norm = max(1.0, frequency**4)
        values = [sum(math.cos(frequency * value / math.sqrt(k)) for value in w) / m for w in spectra]
        mean = sum(values) / len(values)
        variance = sum((value - mean)**2 for value in values) / len(values)
        exact_cosine_mean = math.cos(frequency / math.sqrt(k))**k
        gaussian_mean = math.exp(-frequency * frequency / 2)
        assert abs(mean - exact_cosine_mean) < 2e-13
        mean_bound, variance_bound = test_norm / (6 * k), 4 * test_norm**2 / k
        assert abs(mean - gaussian_mean) <= mean_bound + 2e-13
        assert variance <= variance_bound + 2e-13
        covariance_sum = 0.0
        worst_covariance_error = -math.inf
        for i in range(m):
            for j in range(m):
                agreement = (k + gram[i][j]) // 2
                disagreement = k - agreement
                moment = 0.5 * (
                    math.cos(2 * frequency / math.sqrt(k))**agreement
                    + math.cos(2 * frequency / math.sqrt(k))**disagreement
                )
                covariance = moment - exact_cosine_mean**2
                rho = gram[i][j] / k
                gaussian_covariance = 0.5 * (
                    math.exp(-frequency**2 * (1 + rho))
                    + math.exp(-frequency**2 * (1 - rho))
                ) - gaussian_mean**2
                assert gaussian_covariance >= -2e-13
                assert gaussian_covariance <= test_norm**2 * rho**2 + 2e-13
                assert covariance - gaussian_covariance <= 3 * test_norm**2 / k + 2e-13
                worst_covariance_error = max(worst_covariance_error, covariance - gaussian_covariance)
                covariance_sum += covariance
        assert abs(variance - covariance_sum / m**2) < 3e-13
        results.append({"frequency": frequency, "mean_error_over_bound": abs(mean-gaussian_mean)/mean_bound,
                        "variance_over_bound": variance/variance_bound,
                        "maximum_pair_covariance_error": worst_covariance_error})
    return {"m": m, "k": k, "selector": selector, "all_spin_energies_exact": len(spectra),
            "raw_gram_trace_square_exact": m*m*k, "cosine_checks_numerical": results}


def main():
    cases = [
        case(sylvester(4), [0, 2]),
        case(sylvester(8), [0, 1, 3]),
        case(sylvester(8), [0, 2, 3, 5, 7]),
        case(paley_twelve(), [0, 2, 5, 8]),
        case(paley_twelve(), [0, 1, 4, 6, 8, 11]),
    ]
    print(json.dumps({"cases": cases,
                      "scope": "exact finite matrix identities; cosine checks numerical; no uniform limit inferred"}, indent=2))


if __name__ == "__main__":
    main()
