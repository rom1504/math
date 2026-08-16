#!/usr/bin/env python3
"""Finite checks for balanced exposure and facet-deletion responses."""

from __future__ import annotations

import itertools


def shape_distance(f, g):
    diff = [a - b for a, b in zip(f, g)]
    return (max(diff) - min(diff)) / 2


def mcshane_extension(points, values, metric):
    return tuple(min(values[j] + metric[x][p] for j, p in enumerate(points)) for x in range(len(metric)))


def balanced_exposure_checks():
    checks = 0
    for q in range(2, 11):
        metric = tuple(tuple(abs(i - j) for j in range(q)) for i in range(q))
        for gamma2 in range(1, q):
            # gamma2 represents 2*gamma and is integral.  Greedily choose a
            # non-strict gamma2-separated query set.
            code = []
            for x in range(q):
                if all(metric[x][y] >= gamma2 for y in code):
                    code.append(x)
            if len(code) % 2:
                code.pop()
            if not code:
                continue
            half = len(code) // 2
            functions = []
            for U in itertools.combinations(range(len(code)), half):
                U = set(U)
                vals = tuple(gamma2 / 2 if i in U else -gamma2 / 2 for i in range(len(code)))
                f = mcshane_extension(code, vals, metric)
                assert all(abs(f[x] - f[y]) <= metric[x][y] for x in range(q) for y in range(q))
                for i, x in enumerate(code):
                    assert f[x] == vals[i]
                functions.append(f)
                checks += 1
            for i, f in enumerate(functions):
                for g in functions[i + 1 :]:
                    assert shape_distance(f, g) >= gamma2
                    checks += 1
    return checks


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def facet_deletion_checks():
    checks = 0
    for m in range(2, 8):
        vertices = tuple(itertools.product((0, 1), repeat=m))
        # Facets of the cube: (coordinate, prescribed bit, outward normal).
        facets = []
        for i in range(m):
            plus = tuple(1 if j == i else 0 for j in range(m))
            minus = tuple(-1 if j == i else 0 for j in range(m))
            facets.append((i, 1, plus))
            facets.append((i, 0, minus))

        responses = []
        for target, (i, bit, normal) in enumerate(facets):
            beta = max(dot(v, normal) for v in vertices)
            target_vertices = {v for v in vertices if v[i] == bit}
            gap = beta - max(dot(v, normal) for v in vertices if v not in target_vertices)
            assert gap == 1
            theta = tuple(x / gap for x in normal)
            response = [max(dot(v, theta) for v in vertices)]
            for j, (coord, prescribed, _) in enumerate(facets):
                kept = [v for v in vertices if v[coord] != prescribed]
                response.append(max(dot(v, theta) for v in kept))
                expected = beta - (1 if j == target else 0)
                assert response[-1] == expected
                checks += 1
            responses.append(tuple(response))
            assert shape_distance(response, tuple(0 for _ in response)) == 0.5
        for i, f in enumerate(responses):
            for g in responses[i + 1 :]:
                assert shape_distance(f, g) == 1
                checks += 1
    return checks


def main():
    balanced = balanced_exposure_checks()
    facets = facet_deletion_checks()
    print(f"balanced McShane exposure/separation checks: {balanced}")
    print(f"cube facet-deletion response checks: {facets}")


if __name__ == "__main__":
    main()
