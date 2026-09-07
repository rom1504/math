"""Exact Taylor pairing and ordered-four-cycle multiplicity checks."""

from itertools import combinations, product
import json
import numpy as np


def main():
    four_edges = list(combinations(range(4), 2))
    pairings = 0
    for signs in product((-1, 1), repeat=6):
        a = np.zeros((4, 4), dtype=np.int64)
        for (i, j), sign in zip(four_edges, signs):
            a[i, j] = a[j, i] = sign
        # In L_1^2/2, two orderings cancel the Taylor factor 1/2.
        coefficient_twice = sum(int(a[i, j]*a[k, l])
                                for i, j in four_edges for k, l in four_edges
                                if len({i, j, k, l}) == 4)
        matchings = int(a[0, 1]*a[2, 3]+a[0, 2]*a[1, 3]+a[0, 3]*a[1, 2])
        assert coefficient_twice == 2*matchings
        pairings += 1
    rng = np.random.default_rng(64)
    cycles = 0
    for n in range(4, 11):
        for _ in range(12):
            a = np.triu(rng.choice([-1, 1], size=(n, n)), 1)
            a = a+a.T
            edges = list(combinations(range(n), 2))
            cycle = sum(int(a[i, j]*a[k, l]*(a[i, k]*a[j, l]+a[i, l]*a[j, k]))
                        for i, j in edges for k, l in edges
                        if len({i, j, k, l}) == 4)
            trace_four = int(np.trace(a@a@a@a))
            assert 2*cycle == trace_four-n*(n-1)*(2*n-3)
            cycles += 1
    derivative_bound = (1728+2592*6+1296*54+216*648
                        +2592+3888*6+1296*54+864+1296*6)
    assert derivative_bound == 331776 == 6*55296
    print(json.dumps({"status": "exact checks passed",
                      "four_sign_pairing_checks": pairings,
                      "ordered_cycle_checks": cycles,
                      "integrated_third_derivative_bound": derivative_bound,
                      "uniform_Taylor_remainder_constant": 55296}, indent=2))


if __name__ == "__main__":
    main()
