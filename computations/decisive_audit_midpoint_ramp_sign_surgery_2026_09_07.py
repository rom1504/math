"""Exact rational checks for the finite ramp surgery, not limit numerics."""
from fractions import Fraction as F
from itertools import combinations, product


def matrix(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def quadratic(a, x):
    return sum(a[i][j] * x[i] * x[j]
               for i in range(len(x)) for j in range(len(x)))


def gram(vectors, absolute=False):
    n = len(vectors[0])
    ans = matrix(n)
    for v in vectors:
        norm = sum(z * z for z in v)
        for i in range(n):
            for j in range(n):
                val = v[i] * v[j]
                ans[i][j] += F(abs(val) if absolute else val, norm)
    return ans


def projector_families():
    h = [(1, 1, 1, 1), (1, -1, 1, -1),
         (1, 1, -1, -1), (1, -1, -1, 1)]
    ans = []
    for mask in range(1, 16):
        ans.append([h[i] for i in range(4) if mask & (1 << i)])
    ans.extend([[(1, 0, 0, 0)], [(1, 1, 0, 0)], [(1, 1, 1, 0)],
                [(1, 0, 0, 0), (0, 1, 1, 1)],
                [(1, 1, 0, 0), (0, 0, 1, -1)]])
    return ans


def main():
    n = 4
    states = list(product((-1, 1), repeat=n))
    edges = list(combinations(range(n), 2))
    checks = 0
    means = 0
    variances = 0
    compatible = 0
    families = projector_families()
    for code in range(1 << len(edges)):
        a = matrix(n)
        for e, (i, j) in enumerate(edges):
            a[i][j] = a[j][i] = F(1 - 2 * ((code >> e) & 1))
        energies = [quadratic(a, x) / 2 for x in states]
        if max(energies) < -min(energies):
            a = [[-z for z in row] for row in a]
            energies = [-z for z in energies]
        positive, negative = max(energies), -min(energies)
        w, midpoint = (positive + negative) / 2, (positive - negative) / 2
        t = 2 * midpoint / n
        beta = max(sum(abs(sum(a[i][j] * y[j] for j in range(n)))
                       for i in range(n)) for y in states)
        assert beta <= 4 * w
        for vectors in families:
            pi = gram(vectors)
            r = len(vectors)
            for i in range(n):
                for j in range(n):
                    assert sum(pi[i][z] * pi[z][j] for z in range(n)) == pi[i][j]
            ramp = max(F(0), max(energy - w - midpoint * quadratic(pi, x) / n
                                for energy, x in zip(energies, states)))
            if all(a[i][j] * pi[i][j] >= 0 and t * abs(pi[i][j]) <= 2
                   for i, j in edges):
                probabilities = [t * a[i][j] * pi[i][j] / 2 for i, j in edges]
                total = sum(probabilities)
                assert total == t * sum(a[i][j] * pi[j][i]
                                        for i in range(n) for j in range(n)) / 4
                assert total ** 2 <= midpoint ** 2 * r / 4
                for p, (i, j) in zip(probabilities, edges):
                    assert 0 <= p <= 1
                    assert a[i][j] * (1 - 2 * p) == a[i][j] - t * pi[i][j]
                for energy, x in zip(energies, states):
                    shifted = energy - midpoint * (quadratic(pi, x) - r) / n
                    assert abs(shifted) <= w + ramp + midpoint * r / n
                compatible += 1
            for mu in (F(1, 8), F(1, 4), F(1, 2), F(1)):
                if t * mu > 1:
                    continue
                keep = [pi[i][i] <= mu for i in range(n)]
                removed = n - sum(keep)
                pi0 = [[pi[i][j] * keep[i] * keep[j] for j in range(n)] for i in range(n)]
                bigk = gram(vectors, absolute=True)
                bigk = [[bigk[i][j] * keep[i] * keep[j]
                         for j in range(n)] for i in range(n)]
                assert removed <= F(r) / mu
                trace = sum(pi0[i][i] for i in range(n))
                masked_ramp = max(F(0), max(energy - w - midpoint * quadratic(pi0, x) / n
                                           for energy, x in zip(energies, states)))
                if masked_ramp > ramp:
                    assert (masked_ramp - ramp) ** 2 <= 4 * midpoint ** 2 * F(removed, n)
                probs = [t * (bigk[i][j] + a[i][j] * pi0[i][j]) / 2
                         for i, j in edges]
                assert all(0 <= p <= t * mu <= 1 for p in probs)
                assert sum(probs) <= midpoint * r
                assert sum(4 * p * (1 - p) for p in probs) <= 4 * midpoint * r
                abar = matrix(n)
                penalty = matrix(n)
                schur = matrix(n)
                for p, (i, j) in zip(probs, edges):
                    abar[i][j] = abar[j][i] = a[i][j] * (1 - 2 * p)
                    penalty[i][j] = penalty[j][i] = a[i][j] - t * pi0[i][j]
                    schur[i][j] = schur[j][i] = a[i][j] * bigk[i][j]
                    assert abar[i][j] == penalty[i][j] - t * schur[i][j]
                schur_cap = max(abs(quadratic(schur, x) / 2) for x in states)
                # K_G<=2 is sufficient for this independent finite bound.
                assert schur_cap <= beta * mu
                for energy, x in zip(energies, states):
                    projected = energy - midpoint * quadratic(pi0, x) / n
                    assert projected >= -w
                    assert projected <= w + masked_ramp
                    assert quadratic(penalty, x) / 2 == projected + midpoint * trace / n
                    assert abs(quadratic(abar, x) / 2) <= w + masked_ramp + midpoint * r / n + t * schur_cap
                    means += 1
                checks += 1
            # Exact product-law projector variance checks. Means b*xstar.
            xstar = states[energies.index(positive)]
            for b in (F(0), F(1, 3), F(2, 3), F(1)):
                weights = [product_weight(y, xstar, b) for y in states]
                values = [quadratic(pi, y) for y in states]
                expectation = sum(q * z for q, z in zip(weights, values))
                variance = sum(q * (z - expectation) ** 2 for q, z in zip(weights, values))
                v = 1 - b * b
                assert expectation == b * b * quadratic(pi, xstar) + v * r
                assert variance <= 4 * b * b * v * n + 2 * v * v * r <= 2 * n
                variances += 1
    print("PASS:", checks, "masked flip constructions;", means,
          "exact mean-energy identities;", variances, "product-law variances;",
          compatible, "sign-compatible exact-mean constructions")


def product_weight(y, center, b):
    ans = F(1)
    for yi, ci in zip(y, center):
        ans *= (1 + b * yi * ci) / 2
    return ans


if __name__ == "__main__":
    main()
