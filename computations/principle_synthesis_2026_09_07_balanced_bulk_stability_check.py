"""Exact finite checks for balanced-bulk stability and the Walsh boundary."""
import itertools
import json
from fractions import Fraction

import numpy as np


def walsh(order):
    return np.array([[(-1) ** bin(int(i & j)).count("1") for j in range(order)]
                     for i in range(order)], dtype=np.int64)


def build(order, codimension):
    h = walsh(order)
    transforms = np.array([h[:, np.arange(order) ^ i] for i in range(order)])
    transforms[np.arange(order), :, np.arange(order)] = 0
    seed = h.copy()
    np.fill_diagonal(seed, 0)
    support = range(1, 1 << codimension)
    for i in range(order):
        for difference in support:
            seed[i, i ^ difference] = 1
    w = np.zeros((order * order, order * order), dtype=np.int64)
    for i in range(order):
        for j in range(order):
            w[i * order:(i + 1) * order, j * order:(j + 1) * order] = (
                seed[i, j] * np.outer(transforms[i, :, j], transforms[j, :, i])
            )
    minority = np.array([(a & ((1 << codimension) - 1)) == 0
                         for a in range(order)])
    x = np.tile(1 - 2 * minority.astype(int), order)
    return w, x


def energy(w, x):
    return Fraction(int(x @ w @ x), 2)


def main():
    rng = np.random.default_rng(202609072041)
    results = []
    for order in (4, 8, 16):
        for codimension in range(1, min(3, order.bit_length() - 1) + 1):
            w, x = build(order, codimension)
            n = order * order
            assert np.all(w == w.T)
            assert np.all(np.diag(w) == 0)
            for i in range(order):
                constant = np.zeros(n, dtype=np.int64)
                constant[i * order:(i + 1) * order] = 1
                assert np.all(w @ constant == 0)
            mean = Fraction(int(x[:order].sum()), order)
            centered = [Fraction(int(value)) - mean for value in x]
            field = w @ x
            assert all(Fraction(int(field[a])) == order * centered[a] for a in range(n))
            assert np.all(x * field > 0)
            variance = sum(value * value for value in centered)
            assert energy(w, x) == order * variance / 2
            b = Fraction(7, 16)
            ratio_objective = energy(w, x) + b * sum(
                int(x[i * order:(i + 1) * order].sum()) ** 2 for i in range(order)
            )
            assert ratio_objective > b * order ** 3
            flips_checked = 0
            for unused in range(10):
                word = rng.choice((-1, 1), n)
                fields = w @ word
                for location in rng.choice(n, min(10, n), replace=False):
                    i = int(location) // order
                    old_sum = int(word[i * order:(i + 1) * order].sum())
                    changed = word.copy()
                    changed[location] *= -1
                    actual = energy(w, changed) - energy(w, word)
                    actual += b * ((old_sum - 2 * int(word[location])) ** 2 - old_sum ** 2)
                    expected = -2 * int(word[location]) * int(fields[location])
                    expected += -4 * b * int(word[location]) * old_sum + 4 * b
                    assert actual == expected
                    flips_checked += 1
            results.append({"m": order, "N": n, "minority_fraction": str(Fraction(1, 1 << codimension)),
                            "half_energy": str(energy(w, x)), "variance_ratio": "1/2",
                            "minimum_stability_margin": int(np.min(x * field)),
                            "forced_macro_edges": order * ((1 << codimension) - 1) // 2,
                            "exact_flip_checks": flips_checked, "status": "PASS"})
    print(json.dumps({"status": "PASS", "results": results}, indent=2))


if __name__ == "__main__":
    main()
