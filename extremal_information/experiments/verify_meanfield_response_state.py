#!/usr/bin/env python3
"""Exact finite checks for the heterogeneous binary mean-field state."""

from __future__ import annotations

import itertools
import random
from fractions import Fraction


def profile(fields):
    values = sorted(fields, reverse=True)
    out = [0]
    for value in values:
        out.append(out[-1] + value)
    return tuple(out)


def response(p, lam):
    return max(value + lam * k for k, value in enumerate(p))


def maxplus_convolution(p, q):
    return tuple(
        max(p[i] + q[t - i] for i in range(max(0, t - len(q) + 1), min(len(p) - 1, t) + 1))
        for t in range(len(p) + len(q) - 1)
    )


def quadratic_profile(fields, J):
    p = profile(fields)
    return tuple(value + J * k * (k - 1) / 2 for k, value in enumerate(p))


def upper_hull(q):
    """Canonical vertices of the least concave majorant of integer samples."""
    hull = []
    for x, y in enumerate(q):
        point = (x, Fraction(y))
        hull.append(point)
        while len(hull) >= 3:
            x0, y0 = hull[-3]
            x1, y1 = hull[-2]
            x2, y2 = hull[-1]
            left_slope = (y1 - y0) / (x1 - x0)
            right_slope = (y2 - y1) / (x2 - x1)
            if left_slope > right_slope:
                break
            hull.pop(-2)
    return tuple(hull)


def nearest_grid(value, grid):
    return min(grid, key=lambda x: (abs(x - value), x))


def check_exact_state(rng):
    biconjugacy = 0
    metric = 0
    convolution = 0
    for n in range(1, 10):
        for _ in range(350):
            B = 8
            fields = tuple(rng.randrange(-B, B + 1) for _ in range(n))
            p = profile(fields)
            # Integer fields have integer supporting breakpoints in [-B,B].
            recovered = tuple(min(response(p, lam) - lam * k for lam in range(-B, B + 1)) for k in range(n + 1))
            assert recovered == p
            biconjugacy += n + 1

            other = tuple(rng.randrange(-B, B + 1) for _ in range(n))
            q = profile(other)
            d_profile = max(abs(a - b) for a, b in zip(p, q))
            d_response = max(abs(response(p, lam) - response(q, lam)) for lam in range(-B, B + 1))
            assert d_profile == d_response
            metric += 1

            m = rng.randrange(1, 8)
            tail = tuple(rng.randrange(-B, B + 1) for _ in range(m))
            assert maxplus_convolution(p, profile(tail)) == profile(fields + tail)
            convolution += 1
    return biconjugacy, metric, convolution


def check_histogram_quantization(rng):
    checks = 0
    merge_checks = 0
    B = 12
    spacing = 4
    grid = tuple(range(-B, B + 1, spacing))
    for n in range(1, 30):
        for _ in range(120):
            fields = tuple(rng.randrange(-B, B + 1) for _ in range(n))
            rounded = tuple(nearest_grid(x, grid) for x in fields)
            p, q = profile(fields), profile(rounded)
            assert max(abs(a - b) for a, b in zip(p, q)) <= spacing * n / 2
            for lam in range(-B, B + 1):
                assert abs(response(p, lam) - response(q, lam)) <= spacing * n / 2
                checks += 1

            cut = rng.randrange(n + 1)
            left, right = fields[:cut], fields[cut:]
            hist = lambda xs: tuple(sum(nearest_grid(x, grid) == g for x in xs) for g in grid)
            assert tuple(a + b for a, b in zip(hist(left), hist(right))) == hist(fields)
            merge_checks += 1
    return checks, merge_checks


def check_quadratic_roofs(rng):
    strict = 0
    chord = 0
    sharpness = 0
    composition = 0
    roof_congruence = 0

    # Strict roof collapse: (0,0,J) and (0,a,J) share the endpoint chord.
    for J in range(2, 20):
        for a in range(1, (J - 1) // 2 + 1):
            q0 = quadratic_profile((0, 0), J)
            qa = quadratic_profile((a, -a), J)
            assert q0 != qa
            for lam2 in range(-4 * J, 4 * J + 1):
                lam = Fraction(lam2, 2)
                assert response(q0, lam) == response(qa, lam)
            strict += 1

    # The sharp size-n uniform threshold J >= 4B/n puts every interior
    # point below the endpoint chord.  Just below it, the extremal two-level
    # field list violates the chord at every nontrivial split.
    B = 6
    for n in range(2, 10):
        J = Fraction(4 * B, n)
        for _ in range(250):
            fields = tuple(rng.randrange(-B, B + 1) for _ in range(n))
            q = quadratic_profile(fields, J)
            endpoint = q[-1]
            for k in range(n + 1):
                assert n * q[k] <= k * endpoint
                chord += 1
            for lam in range(-8 * B, 8 * B + 1):
                assert response(q, lam) == max(0, endpoint + lam * n)

        J_below = J - Fraction(1, 100 * n)
        k = n // 2
        fields = (B,) * k + (-B,) * (n - k)
        q = quadratic_profile(fields, J_below)
        assert n * q[k] > k * q[-1]
        sharpness += 1

    # Raw bilinear merge is exact and associative.
    for _ in range(1200):
        J = rng.randrange(-5, 8)
        blocks = [tuple(rng.randrange(-6, 7) for _ in range(rng.randrange(1, 6))) for _ in range(3)]
        qa, qb, qc = (quadratic_profile(block, J) for block in blocks)

        def merge(q, r):
            return tuple(
                max(q[i] + r[t - i] + J * i * (t - i)
                    for i in range(max(0, t - len(r) + 1), min(len(q) - 1, t) + 1))
                for t in range(len(q) + len(r) - 1)
            )

        direct = quadratic_profile(tuple(itertools.chain.from_iterable(blocks)), J)
        assert merge(merge(qa, qb), qc) == direct
        assert merge(qa, merge(qb, qc)) == direct
        composition += 2

    # Equality of child concave roofs is a congruence for every same-J append.
    alphabet = (-2, -1, 0, 1, 2)
    for J in range(1, 6):
        for n in range(1, 4):
            groups = {}
            for fields in itertools.product(alphabet, repeat=n):
                q = quadratic_profile(fields, J)
                groups.setdefault(upper_hull(q), []).append(q)
            for group in groups.values():
                if len(group) < 2:
                    continue
                representatives = group[: min(4, len(group))]
                for qa, qb in itertools.combinations(representatives, 2):
                    for m in range(1, 3):
                        for tail in itertools.product(alphabet, repeat=m):
                            qc = quadratic_profile(tail, J)

                            def merge(q, r):
                                return tuple(
                                    max(q[i] + r[t - i] + J * i * (t - i)
                                        for i in range(max(0, t - len(r) + 1), min(len(q) - 1, t) + 1))
                                    for t in range(len(q) + len(r) - 1)
                                )

                            assert upper_hull(merge(qa, qc)) == upper_hull(merge(qb, qc))
                            roof_congruence += 1
    return strict, chord, sharpness, composition, roof_congruence


def main():
    rng = random.Random(20260816)
    biconjugacy, metric, convolution = check_exact_state(rng)
    quantized, hist_merge = check_histogram_quantization(rng)
    strict, chord, sharpness, quadratic, roof_congruence = check_quadratic_roofs(rng)
    print(f"linear-field biconjugacy coordinates: {biconjugacy}")
    print(f"contextual metric identities: {metric}")
    print(f"max-plus/sorted-union identities: {convolution}")
    print(f"quantized response bounds: {quantized}")
    print(f"histogram merge identities: {hist_merge}")
    print(f"strict quadratic roof collapses: {strict}")
    print(f"sharp-threshold chord inequalities: {chord}")
    print(f"below-threshold chord obstructions: {sharpness}")
    print(f"quadratic associative merges: {quadratic}")
    print(f"same-roof continuation congruences: {roof_congruence}")


if __name__ == "__main__":
    main()
