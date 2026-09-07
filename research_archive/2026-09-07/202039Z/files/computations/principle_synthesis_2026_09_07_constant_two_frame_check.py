"""Exact bounded checks for the constant-two-frame theorem.

No asymptotic inference is drawn from this finite test.  Fractions check
the shared-edge identity, row-repair pushforward, and covariance entries.
"""
import itertools
import json
import math
from collections import defaultdict
from fractions import Fraction


def sign_mean(value):
    return (value > 0) - (value < 0)


def kappa_checks():
    answers = []
    for n in range(3, 10):
        counts = defaultdict(int)
        for bits in itertools.product((-1, 1), repeat=n - 2):
            counts[sum(bits)] += 1
        total = 2 ** (n - 2)
        lam = Fraction(math.comb(n - 2, (n - 2) // 2), total)
        for alpha, beta in itertools.product((-1, 1), repeat=2):
            expectation = Fraction(0)
            for shared in (-1, 1):
                for first, c_first in counts.items():
                    for second, c_second in counts.items():
                        expectation += Fraction(
                            c_first * c_second
                            * sign_mean(first + alpha * shared)
                            * sign_mean(second + beta * shared),
                            2 * total * total,
                        )
            assert expectation == alpha * beta * lam * lam
        answers.append({"n": n, "lambda": str(lam), "kappa": str(lam * lam)})
    return answers


def repair_checks():
    answers = []
    for n in (4, 6, 8):
        pushforward = defaultdict(Fraction)
        states_checked = 0
        for row in itertools.product((-1, 1), repeat=n):
            imbalance = sum(row)
            majority_sign = sign_mean(imbalance)
            majority = [j for j in range(n) if row[j] == majority_sign]
            r = abs(imbalance) // 2
            if r:
                choices = list(itertools.combinations(majority, r))
                p = Fraction(r, len(majority))
                a = majority_sign * p
            else:
                choices = [()]
                p = a = Fraction(0)
            means = [(1 - p) * value - a for value in row]
            observed_mean = [Fraction(0) for _ in row]
            observed_cov = [[Fraction(0) for _ in row] for _ in row]
            for chosen in choices:
                output = list(row)
                for j in chosen:
                    output[j] *= -1
                assert sum(output) == 0
                pushforward[tuple(output)] += Fraction(1, 2 ** n * len(choices))
                centered = [value - mean for value, mean in zip(output, means)]
                assert sum(value * value for value in centered) == 4 * r * (1 - p)
                for j in range(n):
                    observed_mean[j] += Fraction(output[j], len(choices))
                    for k in range(n):
                        observed_cov[j][k] += centered[j] * centered[k] / len(choices)
            assert observed_mean == means
            for j in range(n):
                for k in range(n):
                    expected = Fraction(0)
                    if r and j in majority and k in majority:
                        count = len(majority)
                        expected = 4 * p * (1 - p) * (
                            Fraction(count if j == k else 0, count - 1)
                            - Fraction(1, count - 1)
                        )
                    assert observed_cov[j][k] == expected
            states_checked += 1
        assert len(pushforward) == math.comb(n, n // 2)
        assert set(pushforward.values()) == {Fraction(1, math.comb(n, n // 2))}
        answers.append({"n": n, "input_rows": states_checked,
                        "balanced_outputs": len(pushforward), "status": "PASS"})
    return answers


if __name__ == "__main__":
    print(json.dumps({"kappa": kappa_checks(), "repair": repair_checks(),
                      "constant": (2 / math.pi + math.sqrt(2 / math.pi))
                      / (2 * math.sqrt(2)), "status": "PASS"}, indent=2))
