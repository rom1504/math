"""Small exact checks of the complete-graph endpoint count and repair proof."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import json
from math import comb, factorial


def permanent(matrix):
    size = len(matrix)
    values = [0] * (1 << size)
    values[0] = 1
    for mask in range(1 << size):
        row = bin(mask).count("1")
        if row == size:
            continue
        for column in range(size):
            if not ((mask >> column) & 1):
                values[mask | (1 << column)] += values[mask] * matrix[row][column]
    return values[-1]


def exact_partition(m, ones, kernel):
    edges = list(combinations(range(m), 2))
    degree_left = [m - 1] * m
    states = {tuple(ones): 1}
    for i, j in edges:
        degree_left[i] -= 1
        degree_left[j] -= 1
        following = Counter()
        for remaining, value in states.items():
            for a, b in product((0, 1), repeat=2):
                next_i, next_j = remaining[i] - a, remaining[j] - b
                if not (0 <= next_i <= degree_left[i] and 0 <= next_j <= degree_left[j]):
                    continue
                target = list(remaining)
                target[i], target[j] = next_i, next_j
                following[tuple(target)] += value * kernel[a][b]
        states = following
    numerator = states[(0,) * m]
    denominator = 1
    for value in ones:
        denominator *= comb(m - 1, value)
    return F(numerator, denominator), numerator, denominator


def repair(m, raw, target_ones):
    # Endpoint array indexing is vertex-major; no edge-identification quotient.
    out = list(raw)
    changed = 0
    for i in range(m):
        positions = list(range(i * (m - 1), (i + 1) * (m - 1)))
        excess = sum(out[j] for j in positions) - target_ones[i]
        old, new = (1, 0) if excess > 0 else (0, 1)
        for position in positions:
            if excess == 0:
                break
            if out[position] == old:
                out[position] = new
                excess += -1 if old else 1
                changed += 1
    return tuple(out), changed


def weight(m, values, kernel):
    neighbors = [[j for j in range(m) if j != i] for i in range(m)]
    answer = 1
    for i, j in combinations(range(m), 2):
        a = values[i * (m - 1) + neighbors[i].index(j)]
        b = values[j * (m - 1) + neighbors[j].index(i)]
        answer *= kernel[a][b]
    return answer


def check_repairs():
    m = 3
    targets = [1] * m
    kernel = [[1, 2], [2, 1]]  # Positive but NOT positive semidefinite.
    preimages = Counter()
    changed_histogram = Counter()
    exact_valid_weight = 0
    for raw in product((0, 1), repeat=m * (m - 1)):
        out, changed = repair(m, raw, targets)
        assert all(sum(out[i * 2:(i + 1) * 2]) == 1 for i in range(m))
        assert changed == sum(a != b for a, b in zip(raw, out))
        assert weight(m, out, kernel) * 2**changed >= weight(m, raw, kernel)
        preimages[out] += 1
        changed_histogram[changed] += 1
        if raw == out:
            exact_valid_weight += weight(m, raw, kernel)
    assert len(preimages) == 2**m
    assert max(preimages.values()) <= sum(comb(6, r) for r in range(4))
    expectation, numerator, denominator = exact_partition(m, targets, kernel)
    assert numerator == exact_valid_weight and denominator == 2**m
    return {"raw_configurations": 64, "valid_outputs": len(preimages),
            "change_histogram": dict(sorted(changed_histogram.items())),
            "largest_preimage": max(preimages.values()), "exact_expectation": str(expectation)}


def main():
    cases = []
    for m in range(2, 8):
        ones = [(m - 1) // 2] * m
        for name, kernel, psd in (("psd", [[2, 1], [1, 2]], True),
                                  ("non_psd", [[1, 2], [2, 1]], False)):
            z, numerator, denominator = exact_partition(m, ones, kernel)
            # The constant-kernel run checks the entire denominator count.
            _, counted_denominator, _ = exact_partition(m, ones, [[1, 1], [1, 1]])
            assert counted_denominator == denominator
            row = [0] * (m - 1 - ones[0]) + [1] * ones[0]
            p_row = F(permanent([[kernel[a][b] for b in row] for a in row]), factorial(m - 1))
            if psd:
                assert z * z <= p_row**m
            cases.append({"m": m, "kernel": name, "row_ones": ones[0],
                          "partition": str(z), "finner_squared_ratio": str(z*z / p_row**m)})
    print(json.dumps({"repair": check_repairs(), "partitions": cases,
                      "limiting_pressure_both_kernels": "(1/2) log(3/2)",
                      "arithmetic": "integer and rational only",
                      "scope": "finite counting/repair regression, not an asymptotic proof"}, indent=2))


if __name__ == "__main__":
    main()
