"""Exhaustive rational checks of the universal near-level entropy floor."""

from __future__ import annotations

from fractions import Fraction
import itertools
import json
import math


def main():
    reports = []
    for n in range(2, 6):
        edges = list(itertools.combinations(range(n), 2))
        words = list(itertools.product((-1, 1), repeat=n))
        count = 0
        for signing in itertools.product((-1, 1), repeat=len(edges)):
            values = [sum(a*x[i]*x[j] for a, (i, j) in zip(signing, edges))
                      for x in words]
            cap = max(abs(value) for value in values)
            index = next(i for i, value in enumerate(values) if abs(value) == cap)
            ground, sign = words[index], (1 if values[index] > 0 else -1)
            deficits = [cap-sign*value for value in values]
            assert min(deficits) == 0
            for radius in range(n+1):
                sphere = [i for i, word in enumerate(words)
                          if sum(a != b for a, b in zip(ground, word)) == radius]
                assert len(sphere) == math.comb(n, radius)
                predicted = Fraction(4*radius*(n-radius)*cap, n*(n-1))
                actual = Fraction(sum(deficits[i] for i in sphere), len(sphere))
                assert actual == predicted
                for threshold in range(1, 2*cap+1):
                    near_count = sum(cap-abs(value) <= threshold for value in values)
                    lower = max(Fraction(0), 1-predicted/threshold)*len(sphere)
                    assert near_count >= lower
            count += 1
        reports.append({"n": n, "full_signings_checked": count})
    print(json.dumps({"status": "all exact sphere and Markov checks passed", "cases": reports}, indent=2))


if __name__ == "__main__":
    main()
