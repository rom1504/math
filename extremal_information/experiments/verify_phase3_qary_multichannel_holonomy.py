#!/usr/bin/env python3
"""Exhaustive F_3 audit of the scalar-closed multichannel normal form."""

from __future__ import annotations

from collections import deque
from itertools import product
import json
from pathlib import Path


Q = 3


def add(x, y):
    return tuple((a + b) % Q for a, b in zip(x, y))


def scale(a, x):
    return tuple((a * b) % Q for b in x)


def weight(x):
    return sum(a != 0 for a in x)


def rank(rows):
    rows = [list(row) for row in rows]
    if not rows:
        return 0
    ncol = len(rows[0])
    out = 0
    for col in range(ncol):
        pivot = next((i for i in range(out, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[out], rows[pivot] = rows[pivot], rows[out]
        inv = pow(rows[out][col], -1, Q)
        rows[out] = [(inv * x) % Q for x in rows[out]]
        for i in range(len(rows)):
            if i != out and rows[i][col]:
                a = rows[i][col]
                rows[i] = [(x - a * y) % Q for x, y in zip(rows[i], rows[out])]
        out += 1
    return out


def basis_vector(n, i):
    return tuple(int(j == i) for j in range(n))


def generators(d, k, vectors):
    zero_w = (0,) * d
    zero_q = (0,) * k
    kernel = []
    p_only = []
    r_only = []
    for i in range(d):
        for a in range(1, Q):
            kernel.append(scale(a, basis_vector(d, i)) + zero_q)
    for j in range(k):
        for a in range(1, Q):
            aq = scale(a, basis_vector(k, j))
            p_only.append(zero_w + aq)
            r_only.append(scale(a, vectors[j]) + aq)
    return kernel, p_only, r_only


def bfs_dist(n, gens):
    origin = (0,) * n
    dist = {origin: 0}
    queue = deque([origin])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = add(x, g)
            if y not in dist:
                dist[y] = dist[x] + 1
                queue.append(y)
    return dist


def matvec(vectors, z):
    if not vectors:
        return ()
    out = (0,) * len(vectors[0])
    for a, v in zip(z, vectors):
        out = add(out, scale(a, v))
    return out


def formula(u, vectors):
    k = len(vectors)
    return min(
        2 * weight(z) + weight(add(u, matvec(vectors, z)))
        for z in product(range(Q), repeat=k)
    )


def code(vectors):
    return {
        matvec(vectors, z)
        for z in product(range(Q), repeat=len(vectors))
    }


def dist_to_set(u, points):
    return min(weight(add(u, scale(Q - 1, c))) for c in points)


def hausdorff(c, cp):
    return max(
        max(dist_to_set(x, cp) for x in c),
        max(dist_to_set(x, c) for x in cp),
    )


def exhaustive_normal_form(d=3, k=2):
    words = list(product(range(Q), repeat=d))
    independent = [v for v in product(words, repeat=k) if rank(v) == k]
    checked_profiles = 0
    for vectors in independent:
        kernel, p_only, r_only = generators(d, k, vectors)
        dist = bfs_dist(d + k, kernel + p_only + r_only)

        # The shear fixing W sends every R letter to the corresponding P letter.
        for j, v in enumerate(vectors):
            for a in range(1, Q):
                r = scale(a, v) + scale(a, basis_vector(k, j))
                w, x = r[:d], r[d:]
                lx = matvec(vectors, x)
                sheared = add(w, scale(Q - 1, lx)) + x
                assert sheared == (0,) * d + scale(a, basis_vector(k, j))

        cv = code(vectors)
        for u in words:
            actual = dist[u + (0,) * k]
            predicted = formula(u, vectors)
            assert actual == predicted, (vectors, u, actual, predicted)
            dc = dist_to_set(u, cv)
            assert dc <= actual <= dc + 2 * k
            checked_profiles += 1
    return len(independent), checked_profiles


def grassmannian_response_audit(d=3):
    words = list(product(range(Q), repeat=d))
    nonzero = [v for v in words if any(v)]
    seen = set()
    bases = []
    for v in nonzero:
        c = frozenset(scale(a, v) for a in range(Q))
        if c not in seen:
            seen.add(c)
            bases.append((v, c))

    profiles = {
        v: {u: formula(u, (v,)) for u in words}
        for v, _ in bases
    }
    pairs = 0
    max_abs_geometry_error = 0
    for i, (v, c) in enumerate(bases):
        for vp, cp in bases[i + 1 :]:
            response_distance = max(abs(profiles[v][u] - profiles[vp][u]) for u in words)
            geometry_error = abs(response_distance - hausdorff(c, cp))
            assert geometry_error <= 2
            max_abs_geometry_error = max(max_abs_geometry_error, geometry_error)
            pairs += 1
    return len(bases), pairs, max_abs_geometry_error


def main():
    tuples, profiles = exhaustive_normal_form()
    lines, pairs, geometry_error = grassmannian_response_audit()
    result = {
        "field": 3,
        "normal_form_dimension": {"D": 3, "k": 2},
        "independent_ordered_tuples_checked": tuples,
        "kernel_endpoint_profiles_checked": profiles,
        "projective_lines_checked": lines,
        "line_pairs_checked": pairs,
        "max_abs_response_hausdorff_error": geometry_error,
        "status": "all assertions passed",
    }
    out = Path(__file__).with_name("phase3_qary_multichannel_holonomy_results.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
