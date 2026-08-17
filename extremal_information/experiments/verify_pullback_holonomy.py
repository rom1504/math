#!/usr/bin/env python3
"""Finite checks for pullback refinement, contraction, and rotation holonomy."""

from decimal import Decimal, ROUND_FLOOR, getcontext
from fractions import Fraction
from itertools import product


getcontext().prec = 90
ONE = Decimal(1)
ALPHA = (Decimal(5).sqrt() - ONE) / 2


def cantor_decoder(x):
    if x <= Fraction(1, 3):
        return 3 * x
    if x <= Fraction(2, 3):
        return 2 - 3 * x
    return 3 * x - 2


def cantor_encode(x, bit):
    return (x + 2 * bit) / 3


def cantor_memory_checks():
    checks = 0
    # Depth eight already checks 32,640 distinct pairs.  The proof itself is
    # symbolic; keeping this finite regression modest makes it cheap to rerun.
    for depth in range(1, 9):
        states = {}
        for word in product((0, 1), repeat=depth):
            x = Fraction(0)
            for bit in word:
                x = cantor_encode(x, bit)
            states[word] = x
            decoded = x
            for bit in reversed(word):
                # The current leading ternary digit is the newest encoder.
                assert decoded <= Fraction(1, 3) if bit == 0 \
                    else decoded >= Fraction(2, 3)
                decoded = cantor_decoder(decoded)
            assert decoded == 0
            checks += 1

        words = tuple(states)
        for i, word in enumerate(words):
            for other in words[i + 1:]:
                x, y = states[word], states[other]
                separation = Fraction(0)
                for _ in range(depth):
                    separation = max(separation, abs(x - y))
                    x, y = cantor_decoder(x), cantor_decoder(y)
                assert separation >= Fraction(1, 3)
                checks += 1
    return checks


def helly_tie_checks():
    checks = 0
    for size in range(3, 12):
        # Use c_i=1/size.  Every proper subsystem is a forest and can be
        # integrated; the full residual vector always sums to -1.
        c = [Fraction(1, size)] * size
        r = [Fraction(-1, size)] * size
        increments = [c_i + r_i for c_i, r_i in zip(c, r)]
        assert sum(increments) == 0
        x = [Fraction(0)]
        for increment in increments[:-1]:
            x.append(x[-1] + increment)
        residuals = [x[(i + 1) % size] - x[i] - c[i]
                     for i in range(size)]
        assert residuals == r
        assert sum(residuals) == -1
        assert max(abs(value) for value in residuals) == Fraction(1, size)
        checks += 1
    return checks


def frac_decimal(x):
    return x - x.to_integral_value(rounding=ROUND_FLOOR)


def circle_distance(x, y=Decimal(0)):
    z = abs(frac_decimal(x) - frac_decimal(y))
    return min(z, ONE - z)


def rotation_word(x, length):
    cut = ONE - ALPHA
    return tuple(int(frac_decimal(x + t * ALPHA) >= cut)
                 for t in range(length))


def rotation_refinement_checks():
    checks = 0
    for horizon in range(0, 51):
        # Observations at times 0,...,horizon have the Sturmian boundaries
        # -j alpha, j=0,...,horizon+1.
        boundaries = sorted(frac_decimal(-j * ALPHA)
                            for j in range(horizon + 2))
        assert all(boundaries[i + 1] - boundaries[i] > Decimal("1e-70")
                   for i in range(len(boundaries) - 1))
        samples = []
        for i, left in enumerate(boundaries):
            right = boundaries[(i + 1) % len(boundaries)]
            if i + 1 == len(boundaries):
                right += ONE
            samples.append(frac_decimal((left + right) / 2))
        words = {rotation_word(x, horizon + 1) for x in samples}
        assert len(words) == horizon + 2
        checks += 1

    # Algebraic norm gives ||m alpha|| > 1/(3m).  It forces a macroscopic
    # phase displacement within at most 2m repetitions.
    for m in range(1, 1001):
        z = m * ALPHA
        nearest = Decimal(int(z + Decimal("0.5")))
        assert abs(z - nearest) > ONE / (3 * m)
        assert any(circle_distance(k * z) >= Decimal("0.25")
                   for k in range(1, 2 * m + 1))
        checks += 1
    return checks


def nearest_grid(x, denominator):
    scaled = x * denominator
    q = scaled.numerator // scaled.denominator
    if 2 * (scaled - q) >= 1:
        q += 1
    return Fraction(q, denominator)


def contraction_checks():
    # Two switched affine contractions x -> x/2 and x -> x/2+1/2.
    denominator = 64
    eta = Fraction(1, 2 * denominator)
    rho = Fraction(1, 2)
    bound = eta / (1 - rho)
    checks = 0
    for depth in range(1, 13):
        for word in product((0, 1), repeat=depth):
            x = Fraction(1, 3)
            q = nearest_grid(x, denominator)
            for letter in word:
                x = x / 2 + Fraction(letter, 2)
                q = nearest_grid(q / 2 + Fraction(letter, 2), denominator)
                assert abs(x - q) <= bound
            checks += 1
    return checks


def main():
    cantor = cantor_memory_checks()
    helly = helly_tie_checks()
    rotation = rotation_refinement_checks()
    contraction = contraction_checks()
    print(f"Cantor contextual-packing checks: {cantor}")
    print(f"Helly tie-cycle checks: {helly}")
    print(f"rotation pullback/holonomy checks: {rotation}")
    print(f"contractive net checks: {contraction}")


if __name__ == "__main__":
    main()
