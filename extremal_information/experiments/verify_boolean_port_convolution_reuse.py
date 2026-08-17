#!/usr/bin/env python3
"""Exact finite checks for Boolean-port convolution and reuse laws."""

from fractions import Fraction
from itertools import product


def group(p):
    return [(1,) + tail for tail in product((-1, 1), repeat=p - 1)]


def mul(x, y):
    return tuple(a * b for a, b in zip(x, y))


def kernel(x):
    return Fraction(abs(sum(x)), len(x))


def convolution(mu, nu, points):
    out = {z: Fraction(0) for z in points}
    for x in points:
        for y in points:
            out[mul(x, y)] += mu[x] * nu[y]
    return out


def response(mu, points):
    return {
        eps: sum(mu[s] * kernel(mul(s, eps)) for s in points)
        for eps in points
    }


def action(nu, f, points):
    return {
        eps: sum(nu[t] * f[mul(t, eps)] for t in points)
        for eps in points
    }


def metric(mu, nu, points):
    rmu = response(mu, points)
    rnu = response(nu, points)
    return max(abs(rmu[x] - rnu[x]) for x in points)


def delta(points, at):
    return {x: Fraction(int(x == at)) for x in points}


def main():
    checks = 0
    for p in range(2, 6):
        points = group(p)
        size = len(points)
        e = points[-1]  # product ordering puts the all-plus word last
        uniform = {x: Fraction(1, size) for x in points}

        # Use several exact probability laws, including all point masses.
        measures = [delta(points, x) for x in points]
        measures += [uniform]
        if size >= 2:
            measures += [
                {
                    x: Fraction(int(x == points[0]) + int(x == points[1]), 2)
                    for x in points
                }
            ]

        for mu in measures:
            for nu in measures:
                lhs = response(convolution(mu, nu, points), points)
                rhs = action(nu, response(mu, points), points)
                assert lhs == rhs
                checks += 1

        assert convolution(uniform, uniform, points) == uniform

        # CR.14 is an exact supremum identity, not just a bound at the
        # all-plus query.  Uniform response is constant, and the point-mass
        # response stays in [0,1] and attains 1.
        uniform_response = response(uniform, points)
        constant = next(iter(uniform_response.values()))
        assert set(uniform_response.values()) == {constant}
        assert constant <= Fraction(1, 2)
        assert metric(uniform, delta(points, e), points) == 1 - constant

        # The order-two subgroup is idempotent, whereas diagonal sample reuse
        # always squares to the identity.
        a = (1, -1) + (1,) * (p - 2)
        subgroup = {
            x: Fraction(int(x == e) + int(x == a), 2) for x in points
        }
        assert convolution(subgroup, subgroup, points) == subgroup
        assert metric(subgroup, delta(points, e), points) == Fraction(
            sum(x != y for x, y in zip(a, e)), p
        )

        # Verify the Bernoulli parity formula for several t and powers.
        for t in (Fraction(1, 10), Fraction(1, 4), Fraction(2, 5)):
            base = {
                x: (1 - t if x == e else t if x == a else Fraction(0))
                for x in points
            }
            power = delta(points, e)
            for length in range(1, 7):
                power = convolution(power, base, points)
                q = (1 - (1 - 2 * t) ** length) / 2
                expected = {
                    x: (1 - q if x == e else q if x == a else Fraction(0))
                    for x in points
                }
                assert power == expected
                # This checks the sharp L-factor response calculation, not
                # only the underlying Bernoulli parity distribution.
                assert metric(power, delta(points, e), points) == q * metric(
                    delta(points, a), delta(points, e), points
                )
                assert metric(base, delta(points, e), points) == t * metric(
                    delta(points, a), delta(points, e), points
                )

        # Exact Doeblin test: lambda=alpha*u+(1-alpha)*point mass.
        alpha = Fraction(1, 3)
        lam = {
            x: alpha * uniform[x] + (1 - alpha) * int(x == a)
            for x in points
        }
        for mu in measures:
            for nu in measures:
                contracted = metric(
                    convolution(mu, lam, points),
                    convolution(nu, lam, points),
                    points,
                )
                assert contracted <= (1 - alpha) * metric(mu, nu, points)

        # A declared three-leaf occurrence tree has exactly the iterated
        # convolution marginal.  This explicitly checks the sampling
        # interpretation, not just the response identity.
        if p == 3:
            leaf_laws = (measures[-1], measures[1], uniform)
            enumerated = {x: Fraction(0) for x in points}
            for x, y, z in product(points, repeat=3):
                enumerated[mul(mul(x, y), z)] += (
                    leaf_laws[0][x] * leaf_laws[1][y] * leaf_laws[2][z]
                )
            iterated = convolution(
                convolution(leaf_laws[0], leaf_laws[1], points),
                leaf_laws[2],
                points,
            )
            assert enumerated == iterated

            # Integer row counts multiply under tensor product.  Normalizing
            # the exact pair-count histogram gives precisely convolution.
            left_counts = {x: i + 1 for i, x in enumerate(points)}
            right_counts = {x: 2 * len(points) - i for i, x in enumerate(points)}
            left_total = sum(left_counts.values())
            right_total = sum(right_counts.values())
            pair_counts = {x: 0 for x in points}
            for x, y in product(points, repeat=2):
                pair_counts[mul(x, y)] += left_counts[x] * right_counts[y]
            assert sum(pair_counts.values()) == left_total * right_total
            left_law = {x: Fraction(left_counts[x], left_total) for x in points}
            right_law = {x: Fraction(right_counts[x], right_total) for x in points}
            expected_law = convolution(left_law, right_law, points)
            assert all(
                Fraction(pair_counts[x], left_total * right_total) == expected_law[x]
                for x in points
            )

    print(f"PASS: {checks} exact convolution response identities plus reuse tests")


if __name__ == "__main__":
    main()
