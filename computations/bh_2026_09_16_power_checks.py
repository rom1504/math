"""Exact Walsh-power checks for the 2026-09-16 polynomial-BH audit.

Run with the repository venv. Small-order checks are evidence for the
identities, not substitutes for the proofs in the associated artifact.
"""

import argparse
from collections import Counter
import hashlib
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def bit_count(value):
    return bin(value).count("1")


def fwht(values):
    out = list(values)
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2 * width):
            for j in range(start, start + width):
                left, right = out[j], out[j + width]
                out[j], out[j + width] = left + right, left - right
        width *= 2
    return out


def sylvester_signing(n):
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < n:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    matrix = matrix[:n, :n].copy()
    np.fill_diagonal(matrix, 0)
    return matrix


def random_signing(n, rng):
    matrix = np.triu(rng.choice([-1, 1], size=(n, n)), 1)
    return matrix + matrix.T


def paley_conference(prime):
    assert prime % 4 == 1
    matrix = np.ones((prime + 1, prime + 1), dtype=np.int64)
    np.fill_diagonal(matrix, 0)
    residues = {x * x % prime for x in range(1, prime)}
    for i in range(prime):
        for j in range(i + 1, prime):
            matrix[i + 1, j + 1] = matrix[j + 1, i + 1] = (
                1 if (i - j) % prime in residues else -1
            )
    assert np.array_equal(matrix @ matrix, prime * np.eye(prime + 1, dtype=np.int64))
    return matrix


def stored_witness(n):
    if n <= 10:
        source, key = f"computations/results/exact_m{n}.json", "matrix"
    elif n == 11:
        source, key = "computations/results/nested_10_in_11_cap17.json", "matrix"
    elif n == 12:
        source, key = "computations/results/extension_nested_m11_to_12.json", "parent_matrix"
    elif n == 13:
        source, key = "computations/results/conference_completion_m13.json", "conference_matrix"
    elif n == 14:
        source, key = "computations/results/heuristic_m14_from_conference.json", "matrix"
    else:
        raise ValueError(n)
    payload = json.loads((ROOT / source).read_text())
    matrix = np.asarray(payload[key], dtype=np.int64)[:n, :n]
    assert matrix.shape == (n, n)
    assert np.array_equal(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0)
    assert set(matrix[np.triu_indices(n, 1)]) <= {-1, 1}
    return matrix, source, key


def hafnian(matrix, vertices):
    if not vertices:
        return 1
    first, rest = vertices[0], vertices[1:]
    return sum(
        int(matrix[first, last])
        * hafnian(matrix, tuple(v for v in rest if v != last))
        for last in rest
    )


def energies(matrix):
    n = len(matrix)
    return [
        sum(
            int(matrix[i, j])
            * (-1 if (mask >> i) & 1 else 1)
            * (-1 if (mask >> j) & 1 else 1)
            for i in range(n)
            for j in range(i + 1, n)
        )
        for mask in range(1 << n)
    ]


def mask_of(vertices):
    return sum(1 << v for v in vertices)


def exact_check(matrix):
    n = len(matrix)
    count = 1 << n
    values = energies(matrix)
    coeffs = {}
    for k in (1, 2, 3):
        transformed = fwht([int(value) ** k for value in values])
        assert all(value % count == 0 for value in transformed)
        coeffs[k] = [value // count for value in transformed]
        assert sum(value * value for value in coeffs[k]) * count == sum(
            int(value) ** (2 * k) for value in values
        )
        for mask, value in enumerate(coeffs[k]):
            size = bit_count(mask)
            assert value == 0 or (size % 2 == 0 and size <= 2 * k)
        if 2 * k <= n:
            for vertices in itertools.combinations(range(n), 2 * k):
                haf = hafnian(matrix, vertices)
                assert abs(haf) % 2 == 1
                assert coeffs[k][mask_of(vertices)] == math.factorial(k) * haf

    a2, a3 = matrix @ matrix, matrix @ matrix @ matrix
    edges = n * (n - 1) // 2
    assert coeffs[2][0] == edges
    assert coeffs[3][0] == int(np.trace(a3))
    for i, j in itertools.combinations(range(n), 2):
        mask = (1 << i) | (1 << j)
        assert coeffs[2][mask] == 2 * int(a2[i, j])
        expected = 6 * int(a3[i, j]) + (3 * edges - 12 * n + 16) * int(
            matrix[i, j]
        )
        assert coeffs[3][mask] == expected
    for vertices in itertools.combinations(range(n), 4):
        correlation = 0
        for i, j in itertools.combinations(vertices, 2):
            left, right = (v for v in vertices if v != i and v != j)
            correlation += int(matrix[left, right]) * int(a2[i, j])
        stars = sum(
            math.prod(int(matrix[v, w]) for w in vertices if w != v)
            for v in vertices
        )
        assert coeffs[3][mask_of(vertices)] == 6 * (correlation - 2 * stars)
    return {"n": n, "cap": max(abs(value) for value in values), "identities": "PASS"}


def pnorm(values, exponent):
    return sum(abs(value) ** exponent for value in values) ** (1 / exponent)


def numerical_norm_check(matrix, powers):
    n = len(matrix)
    values = energies(matrix)
    cap = max(abs(value) for value in values)
    normalized = [value / cap for value in values]
    rows = []
    for k in powers:
        c = [value / len(values) for value in fwht([value ** k for value in normalized])]
        p = 4 * k / (2 * k + 1)
        # Off-support roundoff is removed using the exact support identity.
        c = [value if bit_count(mask) % 2 == 0 and bit_count(mask) <= 2 * k else 0.0
             for mask, value in enumerate(c)]
        moment = (sum(abs(value) ** (2 * k) for value in normalized) / len(values)) ** (1 / (2 * k))
        functional = pnorm(c, p) ** (1 / k)
        d = sum(math.comb(n, r) for r in range(0, min(n, 2 * k) + 1, 2))
        upper_factor = d ** (1 / (4 * k * k))
        squared_sum = sum(value * value for value in c)
        alpha = p / 2
        probabilities = [value * value / squared_sum for value in c if value]
        renyi = math.log(sum(value ** alpha for value in probabilities)) / (1 - alpha)
        entropy_ratio = math.exp(renyi / (4 * k * k))
        weighted_square = 0.0
        for r in range(1, min(n, 2 * k) + 1):
            level = [value for mask, value in enumerate(c) if bit_count(mask) == r]
            q = 2 * r / (r + 1)
            weighted_square += pnorm(level, q) ** 2 / r ** 10
        weighted_root = weighted_square ** (1 / (2 * k))
        ratio = functional / moment
        assert ratio >= 1 - 1e-10
        assert ratio <= upper_factor + 1e-10
        assert abs(ratio - entropy_ratio) <= 1e-9
        assert weighted_root <= (math.e * n) ** (1 / (2 * k)) * moment + 1e-10
        rows.append({
            "k": k,
            "moment_over_cap": moment,
            "power_functional_over_cap": functional,
            "functional_over_moment": ratio,
            "entropy_formula_ratio": entropy_ratio,
            "support_upper_factor": upper_factor,
            "weighted_root_over_cap": weighted_root,
        })
    return {"n": n, "cap": cap, "rows": rows}


def conference_cancellation_check(prime):
    matrix = paley_conference(prime)
    n = len(matrix)
    values = energies(matrix)
    edges = n * (n - 1) // 2
    hafnian_rows = []
    for size in sorted(set([4, 6, min(8, n)])):
        histogram = Counter(hafnian(matrix, vertices)
                            for vertices in itertools.combinations(range(n), size))
        matchings = math.prod(range(1, size, 2))
        hafnian_rows.append({
            "boundary_size": size,
            "number_of_matching_terms_per_minor": matchings,
            "hafnian_histogram": dict(sorted(histogram.items())),
            "sum_squared_hafnians": sum(haf * haf * count for haf, count in histogram.items()),
            "random_sign_expected_sum_squared_hafnians": math.comb(n, size) * matchings,
            "max_hafnian_abs_over_total_matching_terms": max(abs(haf) for haf in histogram) / matchings,
        })
    pairing_rows = []
    for k in sorted(set([max(2, n // 2), n, min(2 * n, edges)])):
        exact_sum = sum(int(value) ** (2 * k) for value in values)
        moment = exact_sum // len(values)
        assert exact_sum % len(values) == 0
        pairings = math.comb(edges, k) * math.factorial(2 * k) // (2 ** k)
        pairing_rows.append({
            "k": k,
            "moment": str(moment),
            "positive_exact_pairing_contribution": str(pairings),
            "pairing_to_full_moment_ratio": pairings / moment,
            "remaining_eulerian_sum": str(moment - pairings),
        })
    return {
        "family": "Paley symmetric conference",
        "prime": prime,
        "n": n,
        "conference_identity": "PASS",
        "cap": max(abs(value) for value in values),
        "hafnian_rows": hafnian_rows,
        "pairing_rows": pairing_rows,
    }


def sign_averaging_check(n, max_k=8):
    edge_list = list(itertools.combinations(range(n), 2))
    edges = len(edge_list)
    total_moments = {k: 0 for k in range(2, max_k + 1)}
    total_z = Fraction(0)
    best_z, best_matrix = None, None
    for sign_mask in range(1 << edges):
        matrix = np.zeros((n, n), dtype=np.int64)
        for index, (i, j) in enumerate(edge_list):
            matrix[i, j] = matrix[j, i] = 1 if (sign_mask >> index) & 1 else -1
        values = energies(matrix)
        z = Fraction(0)
        for k in range(2, max_k + 1):
            numerator = sum(int(value) ** (2 * k) for value in values)
            assert numerator % (1 << n) == 0
            moment = numerator // (1 << n)
            total_moments[k] += moment
            gaussian_moment = math.prod(range(1, 2 * k, 2)) * edges ** k
            z += Fraction(moment, (2 ** k) * gaussian_moment)
        total_z += z
        if best_z is None or z < best_z:
            best_z, best_matrix = z, matrix.tolist()
    averaged_rows = []
    for k in range(2, max_k + 1):
        averaged = Fraction(total_moments[k], 1 << edges)
        independent = Fraction(sum(math.comb(edges, count) * (edges - 2 * count) ** (2 * k)
                                   for count in range(edges + 1)), 1 << edges)
        assert averaged == independent
        gaussian_moment = math.prod(range(1, 2 * k, 2)) * edges ** k
        assert averaged <= gaussian_moment
        averaged_rows.append({"k": k, "averaged_moment": str(averaged),
                              "independent_sign_moment": str(independent),
                              "gaussian_upper_moment": str(gaussian_moment)})
    average_z = total_z / (1 << edges)
    assert best_z <= average_z <= Fraction(1, 2)
    return {"n": n, "all_signings_enumerated": 1 << edges, "max_k": max_k,
            "average_z": str(average_z), "best_z": str(best_z),
            "best_control_matrix": best_matrix, "averaged_moments": averaged_rows,
            "scope": "Finite check of the coefficient-sign averaging and simultaneous-moment argument."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20260916)
    exact = []
    for n in range(2, 11):
        exact.append({"family": "Sylvester", **exact_check(sylvester_signing(n))})
        for sample in range(3):
            exact.append({"family": "random", "sample": sample, **exact_check(random_signing(n, rng))})
    numerical = [
        {"family": "Sylvester", **numerical_norm_check(sylvester_signing(n), [1, 2, 3, 4, 8, 16, 32, 64])}
        for n in [6, 8, 10, 12]
    ]
    stored = []
    expected_caps = {6: 5, 7: 9, 8: 10, 9: 12, 10: 13, 11: 17, 12: 18, 13: 20, 14: 21}
    for n in range(6, 15):
        matrix, source, key = stored_witness(n)
        check = exact_check(matrix)
        assert check["cap"] == expected_caps[n]
        stored.append({
            "source": source,
            "source_key": key,
            "matrix_sha256_int8_row_major": hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest(),
            "scope": "Stored minimum-cap witness, cap rechecked exhaustively; global optimality not reproved.",
            **check,
            "norm_checks": numerical_norm_check(matrix, [2, 3, n // 2, n, 2 * n]),
        })
    contraction = (1 - 2 * math.exp(-8 / 3)) ** (-13 / 24) * 13 ** (1 / 24) * 3 * math.sqrt(2) * (2 / 3) ** 5
    result = {
        "status": "PASS",
        "scope": "Finite verification only; asymptotic theorem is proved in artifact.",
        "seed": 20260916,
        "primary_contraction_constant": contraction,
        "exact_checks": exact,
        "numerical_norm_checks": numerical,
        "stored_witness_checks": stored,
        "conference_cancellation_checks": [conference_cancellation_check(prime) for prime in [5, 13]],
        "coefficient_sign_averaging_checks": [sign_averaging_check(n) for n in [4, 5]],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "random_sylvester_exact_matrices": len(exact), "stored_witnesses": len(stored), "norm_rows": sum(len(row["rows"]) for row in numerical) + sum(len(row["norm_checks"]["rows"]) for row in stored), "conference_orders": [6, 14], "contraction": contraction, "output": str(args.output)}))


if __name__ == "__main__":
    main()
