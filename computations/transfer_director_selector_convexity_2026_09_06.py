"""Exact integer regression for the selector-convexity proof; no asymptotic claim."""

from itertools import combinations, product
import json


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def mv(a, x):
    return [dot(row, x) for row in a]


def main():
    d = 4
    ell = d-1  # exact operator bound by the maximum absolute row sum
    edges = list(combinations(range(d), 2))
    checks = 0
    for entries in product((-1, 1), repeat=len(edges)):
        a = [[0]*d for _ in range(d)]
        for (i, j), value in zip(edges, entries):
            a[i][j] = a[j][i] = value
        for n in range(1, d):
            indicators = [tuple(int(i in subset) for i in range(d))
                          for subset in combinations(range(d), n)]
            for x in product((-1, 1), repeat=d):
                for sigma in (-1, 1):
                    g = [[sigma*x[i]*a[i][j]*x[j] + int(i == j)*ell
                          for j in range(d)] for i in range(d)]
                    for u in indicators:
                        grad = mv(g, u)
                        assert dot(grad, grad) <= 4*ell*ell*n
                        twice_f_u = dot(u, grad)
                        for v in indicators:
                            diff = [ui-vi for ui, vi in zip(u, v)]
                            twice_f_v = dot(v, mv(g, v))
                            assert twice_f_u-twice_f_v <= 2*dot(grad, diff)
                            assert dot(grad, diff) <= sum(abs(gi) for gi, di in zip(grad, diff) if di)
                            assert dot(diff, mv(g, diff)) >= 0
                            checks += 1
    print(json.dumps({"status": "PASS", "evidence": "exact integer finite regression, not a replacement for the written proof",
                      "all_order4_signings": 64, "tangent_pair_checks": checks,
                      "concentration_denominator": 144*4}, indent=2))


if __name__ == "__main__":
    main()
