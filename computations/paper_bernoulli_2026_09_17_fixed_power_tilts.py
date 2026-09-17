"""Exact Walsh square-tilt identities; finite fixed-power diagnostics."""

import json
import math
import numpy as np


def fwht(values):
    out = np.array(values, copy=True)
    step = 1
    while step < len(out):
        view = out.reshape(-1, 2 * step)
        left = view[:, :step].copy()
        right = view[:, step:].copy()
        view[:, :step] = left + right
        view[:, step:] = left - right
        step *= 2
    return out


def cube(n):
    indices = np.arange(2 ** n, dtype=np.int64)
    spins = 1 - 2 * ((indices[:, None] >> np.arange(n)) & 1)
    kernel = np.abs(spins.sum(axis=1))
    return spins, kernel


def sylvester(order):
    result = np.ones((1, 1), dtype=np.int64)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    return result


def main():
    rng = np.random.default_rng(20260917)
    checked = 0
    matrices = 0
    fourier_checks = 0
    exact_rows = []
    for n in range(5, 13):
        spins, kernel = cube(n)
        transform_kernel = fwht(kernel)
        mass_kernel = int(kernel.sum())
        d = n * (n - 1) // 2
        m = n - 1 if n % 2 == 0 else n
        assert int(transform_kernel[3]) * m == mass_kernel
        assert int(transform_kernel[15]) * m * (m - 2) == -mass_kernel
        fourier_checks += 2
        for _ in range(10):
            raw = rng.choice((-1, 1), size=(n, n))
            matrix = np.triu(raw, 1)
            matrix += matrix.T.copy()
            fields = spins @ matrix
            energies = (fields * spins).sum(axis=1) // 2
            cap = int(np.abs(energies).max())
            squared = energies ** 2
            assert int(squared.sum()) == 2 ** n * d
            raw_convolution = fwht(fwht(squared) * transform_kernel)
            assert np.all(raw_convolution % (2 ** n) == 0)
            convolution = raw_convolution // (2 ** n)
            row_square = (fields ** 2).sum(axis=1) - n * (n - 1)
            right = mass_kernel * (d * m * (m - 2) + (m - 1) * row_square - squared + d)
            assert np.array_equal(convolution * m * (m - 2), right)
            assert np.max(np.abs(fields).sum(axis=1)) <= 4 * cap
            assert np.max((fields ** 2).sum(axis=1)) <= 4 * (n - 1) * cap
            checked += 2 ** n
            matrices += 1
        exact_rows.append({"order": n, "matrices": 10, "all_queries_per_matrix": 2 ** n})
    diagnostics = []
    for n in (8, 12, 16):
        ambient = 1
        while ambient < n:
            ambient *= 2
        matrix = sylvester(ambient)[:n, :n].copy()
        np.fill_diagonal(matrix, 0)
        spins, kernel = cube(n)
        fields = spins @ matrix
        energies = (fields * spins).sum(axis=1) // 2
        cap = int(np.abs(energies).max())
        transformed_kernel = fwht(kernel.astype(float))
        power_results = []
        for power in (1, 2, 3, 4):
            weights = energies.astype(float) ** (2 * power)
            convolution = fwht(fwht(weights) * transformed_kernel) / (2 ** n)
            response = convolution / weights.sum() / math.sqrt(n)
            assert np.min(response) >= -1e-10
            power_results.append({"power": power, "minimum_normalized_response": float(response.min()),
                                  "maximum_normalized_response": float(response.max()),
                                  "maximum_deviation_from_kappa": float(np.max(np.abs(response - math.sqrt(2 / math.pi))))})
        diagnostics.append({"order": n, "normalized_cap": cap / n ** 1.5,
                            "power_responses_float_diagnostic_only": power_results})
    print(json.dumps({"status": "PASS", "exact_matrices": matrices,
                      "exact_all_query_walsh_and_cap_checks": checked,
                      "exact_fourier_coefficient_checks": fourier_checks,
                      "exact_orders": exact_rows, "finite_diagnostics": diagnostics},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
