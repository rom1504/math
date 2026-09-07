"""Exact local tree covariance checks for the sparse amplitude theorem."""
from collections import deque
from fractions import Fraction
import math
import random


def edge_tree(d, radius):
    adj = [[1], [0]]
    frontier = [0, 1]
    for _ in range(radius):
        following = []
        for v in frontier:
            for _ in range(d - 1):
                w = len(adj)
                adj.append([v])
                adj[v].append(w)
                following.append(w)
        frontier = following
    return adj


def levels_and_paths(adj, edge_signs, root):
    dist = [-1] * len(adj)
    path = [0] * len(adj)
    dist[root], path[root] = 0, 1
    queue = deque([root])
    while queue:
        v = queue.popleft()
        for w in adj[v]:
            if dist[w] < 0:
                dist[w] = dist[v] + 1
                path[w] = path[v] * edge_signs[min(v, w), max(v, w)]
                queue.append(w)
    return dist, path


def main():
    rng = random.Random(20260907)
    covariance_checks = variance_checks = wave_checks = 0
    for d in range(2, 6):
        for radius in range(2, 6):
            adj = edge_tree(d, radius)
            for repeat in range(12):
                edge_signs = {(i, j): rng.choice((-1, 1))
                              for i, row in enumerate(adj) for j in row if i < j}
                b = [Fraction(0)] + [Fraction(rng.randint(-5, 5), 7)
                                      for _ in range(radius)]
                rows = []
                level_rows = []
                for root in (0, 1):
                    dist, path = levels_and_paths(adj, edge_signs, root)
                    row = [b[k] * s if 1 <= k <= radius else Fraction(0)
                           for k, s in zip(dist, path)]
                    rows.append(row)
                    level_rows.append((dist, path))
                    expected = sum(Fraction(d * (d - 1) ** (k - 1)) * b[k] ** 2
                                   for k in range(1, radius + 1))
                    assert sum(x * x for x in row) == expected
                    variance_checks += 1
                covariance = sum(x * y for x, y in zip(*rows))
                expected = 2 * edge_signs[0, 1] * sum(
                    (d - 1) ** k * b[k] * b[k + 1] for k in range(1, radius))
                assert covariance == expected
                covariance_checks += 1
                a = [0.] + [math.sqrt(2 / (radius + 1))
                            * math.sin(k * math.pi / (radius + 1))
                            for k in range(1, radius + 1)]
                coeff = [0.] + [a[k] / math.sqrt(d * (d - 1) ** (k - 1))
                                for k in range(1, radius + 1)]
                rows = [[coeff[k] * s if 1 <= k <= radius else 0.
                         for k, s in zip(*level)] for level in level_rows]
                assert abs(sum(x * x for x in rows[0]) - 1) < 2e-13
                corr = sum(x * y for x, y in zip(*rows))
                target = (edge_signs[0, 1] * 2 * math.sqrt(d - 1) / d
                          * math.cos(math.pi / (radius + 1)))
                assert abs(corr - target) < 2e-13
                wave_checks += 1
    # Fresh coordinates: old Gaussian labels may remain in the good rows.
    rows = [[Fraction(1, 3), Fraction(2, 3), Fraction(-1, 3), Fraction(0), Fraction(0)],
            [Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(0)],
            [Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1)]]
    assert all(sum(x * y for x, y in zip(rows[i], rows[j])) == 0
               for i in range(3) for j in range(i))
    print({'status': 'PASS', 'exact_variance_checks': variance_checks,
           'exact_covariance_checks': covariance_checks,
           'normalized_wave_checks': wave_checks,
           'fresh_bad_coordinate_checks': 3})


if __name__ == '__main__':
    main()
