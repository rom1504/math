"""Integer/rational sanity checks for global monotone polarity balancing.

No small-order optimality or asymptotic concentration claim is inferred from
these checks. Run with the repository .venv/bin/python.
"""

from fractions import Fraction
import json
import math
import numpy as np


def spins(n):
    words = np.arange(1 << (n - 1), dtype=np.int64)
    result = np.ones((len(words), n), dtype=np.int64)
    for j in range(1, n):
        result[:, j] = 1 - 2 * ((words >> (j - 1)) & 1)
    return result


def check_case(matrix, rng, label):
    n = len(matrix)
    x = spins(n)
    energies = np.sum((x @ matrix) * x, axis=1) // 2
    p, r = int(energies.max()), -int(energies.min())
    if p < r:
        matrix = -matrix
        energies = -energies
        p, r = r, p
    original = energies.copy()
    y = x[int(np.argmin(energies))].copy()
    y_energy = ((x @ y) ** 2 - n) // 2
    eligible = [(i, j) for i in range(n) for j in range(i + 1, n)
                if matrix[i, j] == y[i] * y[j]]
    d = n * (n - 1) // 2
    length = len(eligible)
    assert 2 * length == d - r
    assert 2 * p >= n
    rng.shuffle(eligible)
    first_crossing = 0 if p <= r else None
    last_gap = p - r
    rows = []
    for t in range(length + 1):
        if t:
            i, j = eligible[t - 1]
            energies -= 2 * matrix[i, j] * x[:, i] * x[:, j]
        actual_p = int(energies.max())
        actual_r = -int(energies.min())
        assert actual_r == r + 2 * t
        gap = actual_p - actual_r
        if t:
            assert gap - last_gap in (-4, -2, 0)
        last_gap = gap
        if first_crossing is None and gap <= 0:
            first_crossing = t
            assert abs(gap) <= 2
        if length:
            mean_numerators = (length - t) * original - t * y_energy
            assert int(mean_numerators.max()) <= length * p
            assert -int(mean_numerators.min()) == length * (r + 2 * t)
            deviation = int(np.max(np.abs(length * energies - mean_numerators)))
            rows.append((t, Fraction(deviation, length)))
    assert first_crossing is not None
    crossing_error = max((dev for t, dev in rows if t <= first_crossing), default=Fraction(0))
    crossing_r = r + 2 * first_crossing
    assert crossing_r <= p + crossing_error + 2
    return {"label": label, "n": n, "P": p, "R": r,
            "eligible_edges": length, "first_crossing": first_crossing,
            "crossing_cap": crossing_r,
            "prefix_error_to_crossing": str(crossing_error)}


def asymptotic_parameter_check():
    # Check only the explicit algebraic sufficient inequalities for idealized
    # real upper choices Delta,R<=C n^1.5; these are NOT signing realizations.
    records = []
    for c in (0.5, 1.0, 2.0):
        for n in (10 ** 8, 10 ** 10, 10 ** 12):
            delta = c * n ** 1.5
            r = c * n ** 1.5
            d = n * (n - 1) / 2
            length = (d - r) / 2
            width = math.sqrt(n * (delta + n))
            count = math.ceil(delta / 2 + 64 * width)
            u = (n + 1) * math.log(2) + math.log(8 * (count + 1) * (length + 1))
            error = math.sqrt(8 * count * u) + 4 * u / 3
            assert length >= n ** 2 / 8
            assert count <= length
            assert u <= 2 * n
            assert error < 40 * width
            assert 2 * count > delta + error
            records.append({"C": c, "n": n, "error_over_w": error / width,
                            "endpoint_margin_over_w": (2 * count - delta - error) / width})
    return records


def main():
    rng = np.random.default_rng(20260907)
    records = []
    for n in range(3, 13):
        for trial in range(8):
            raw = rng.choice(np.array([-1, 1], dtype=np.int64), size=(n, n))
            matrix = np.triu(raw, 1)
            matrix += matrix.T
            records.append(check_case(matrix, rng, f"random_{n}_{trial}"))
        positive = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
        records.append(check_case(positive, rng, f"positive_clique_{n}"))
    output = {"status": "PASS", "seed": 20260907,
              "integer_paths_checked": len(records),
              "nonzero_initial_gap_paths": sum(x["P"] != x["R"] for x in records),
              "max_n": 12, "cases": records,
              "large_n_parameter_algebra": asymptotic_parameter_check()}
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
