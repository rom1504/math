#!/usr/bin/env python3
"""Deterministic finite checks for the 2026-09-05 campaign publication.

No files are written. Run with .venv/bin/python -B, redirect stdout only
when a writable workspace is available. Integer enumeration proves the
finite minima; floating-point evaluations only check displayed analytic
bounds numerically. The asymptotic claims rely on the accompanying proofs.
"""
import itertools
import json
import math
from collections import Counter
from fractions import Fraction


def popcount(value):
    return bin(value).count("1")


def spins(n):
    return list(itertools.product((-1, 1), repeat=n))


def mu(n):
    return Fraction(n * math.comb(n - 1, (n - 1) // 2), 2 ** (n - 1))


def rectangle_checks():
    answer = []
    for n in range(1, 5):
        for r in range(1, 5):
            ys = [(1,) + tail for tail in spins(r - 1)]
            best = n * r
            tested = 0
            for flat in spins(n * r):
                value = max(
                    sum(abs(sum(flat[i * r + j] * y[j] for j in range(r)))
                        for i in range(n))
                    for y in ys
                )
                assert value >= n * mu(r)
                assert value >= r * mu(n)
                best = min(best, value)
                tested += 1
            upper = min(
                n * float(mu(r)) + math.sqrt(2 * n * r * (r - 1) * math.log(2)),
                r * float(mu(n)) + math.sqrt(2 * n * r * (n - 1) * math.log(2)),
            )
            assert best <= upper + 1e-10
            answer.append(dict(n=n, r=r, matrices=tested, exact_minimum=best,
                               lower=str(max(n * mu(r), r * mu(n))),
                               analytic_upper_numeric=upper))
    return answer


def codewords(generators):
    code = {0}
    for g in generators:
        code |= {c ^ g for c in code}
    return sorted(code)


def code_data(generators, length):
    code = codewords(generators)
    dual = [v for v in range(1 << length)
            if all(popcount(v & c) % 2 == 0 for c in code)]
    distances = [min(popcount(v ^ c) for c in code)
                 for v in range(1 << length)]
    leaders = Counter(distances)
    assert all(count % len(code) == 0 for count in leaders.values())
    return dict(weight=dict(sorted(Counter(popcount(c) for c in code).items())),
                dual_weight=dict(sorted(Counter(popcount(c) for c in dual).items())),
                leader_weight={d: count // len(code) for d, count in sorted(leaders.items())},
                radius=max(distances))


def radial_pair_check():
    a = code_data([3, 5, 57], 6)
    b = code_data([3, 12, 48], 6)
    assert a["weight"] == b["weight"] == {0: 1, 2: 3, 4: 3, 6: 1}
    assert a["dual_weight"] == b["dual_weight"]
    assert a["leader_weight"] == {0: 1, 1: 4, 2: 3}
    assert b["leader_weight"] == {0: 1, 1: 3, 2: 3, 3: 1}
    assert a["radius"] == 2 and b["radius"] == 3
    return dict(A6=a, B6=b)


def cut_code(n):
    edges = list(itertools.combinations(range(n), 2))
    full = (1 << len(edges)) - 1
    cuts = {sum((x[i] != x[j]) << k for k, (i, j) in enumerate(edges))
            for x in spins(n)}
    return sorted(cuts | {c ^ full for c in cuts})


def code_geometry_checks():
    result = []
    for n in range(3, 6):
        edges = list(itertools.combinations(range(n), 2))
        N = len(edges)
        code = cut_code(n)
        assert len(code) == 1 << n
        distances = [min(popcount(u ^ c) for c in code)
                     for u in range(1 << N)]
        rho = max(distances)
        deepest = distances.index(rho)
        leader = min((deepest ^ c for c in code), key=popcount)
        support = [i for i in range(N) if leader >> i & 1]
        subsets = [sum(1 << support[j] for j in range(rho) if v >> j & 1)
                   for v in range(1 << rho)]
        for f in subsets:
            assert distances[leader ^ f] == rho - popcount(f)
            for g in subsets:
                assert distances[f ^ g] == popcount(f ^ g)
        for u in range(1 << N):
            q = max(abs(sum((-1 if u >> k & 1 else 1) * x[i] * x[j]
                            for k, (i, j) in enumerate(edges))) for x in spins(n))
            assert q == N - 2 * distances[u]
        result.append(dict(n=n, length=N, dimension=n, radius=rho,
                           exact_M=N - 2 * rho, leader_cube_size=len(subsets)))
    return result


def reversal_check():
    n, r = 2, 3
    edges = list(itertools.combinations(range(n + r), 2))
    tested = 0
    for signs in spins(len(edges)):
        a = dict(zip(edges, signs))
        full = max(abs(sum(a[i, j] * z[i] * z[j] for i, j in edges))
                   for z in spins(n + r))
        joint = 0
        for x in spins(n):
            for y in spins(r):
                child = sum(a[i, j] * x[i] * x[j]
                            for i, j in itertools.combinations(range(n), 2))
                child += sum(a[n + i, n + j] * y[i] * y[j]
                             for i, j in itertools.combinations(range(r), 2))
                bridge = sum(a[i, n + j] * x[i] * y[j]
                             for i in range(n) for j in range(r))
                joint = max(joint, abs(child) + abs(bridge))
        assert full == joint
        tested += 1
    return dict(split=[n, r], exact_signings_checked=tested)


def main():
    for r in range(1, 17):
        assert mu(r) == Fraction(sum(abs(sum(y)) for y in spins(r)), 1 << r)
    output = dict(
        status="All assertions passed; exact enumeration with numerical evaluation of analytic upper bounds",
        mu_formula_orders=list(range(1, 17)),
        rectangles=rectangle_checks(),
        radial_pair=radial_pair_check(),
        cut_code_geometry=code_geometry_checks(),
        bridge_reversal=reversal_check(),
    )
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
