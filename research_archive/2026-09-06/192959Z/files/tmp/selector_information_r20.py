#!/usr/bin/env python3
"""Exact Q(sqrt(17)) audit for the Wave 20 A9 ground-channel certificate."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np

from check_response_dual_r16 import A9, spins


@dataclass(frozen=True)
class Quad:
    """The exact number a + b*sqrt(17)."""

    a: Fraction
    b: Fraction = Fraction(0)

    def __add__(self, other: object) -> "Quad":
        z = as_quad(other)
        return Quad(self.a + z.a, self.b + z.b)

    __radd__ = __add__

    def __neg__(self) -> "Quad":
        return Quad(-self.a, -self.b)

    def __sub__(self, other: object) -> "Quad":
        return self + (-as_quad(other))

    def __rsub__(self, other: object) -> "Quad":
        return as_quad(other) - self

    def __mul__(self, other: object) -> "Quad":
        z = as_quad(other)
        return Quad(self.a * z.a + 17 * self.b * z.b,
                    self.a * z.b + self.b * z.a)

    __rmul__ = __mul__

    def inverse(self) -> "Quad":
        den = self.a * self.a - 17 * self.b * self.b
        assert den != 0
        return Quad(self.a / den, -self.b / den)

    def __truediv__(self, other: object) -> "Quad":
        return self * as_quad(other).inverse()

    def __rtruediv__(self, other: object) -> "Quad":
        return as_quad(other) / self

    def sign(self) -> int:
        """Return the exact sign, comparing rational squares when needed."""
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if (self.a > 0) == (self.b > 0):
            return 1 if self.a > 0 else -1
        comparison = self.a * self.a - 17 * self.b * self.b
        if comparison == 0:
            return 0
        if self.a > 0:
            return 1 if comparison > 0 else -1
        return -1 if comparison > 0 else 1

    def value(self) -> float:
        return float(self.a) + float(self.b) * math.sqrt(17)


def as_quad(value: object) -> Quad:
    if isinstance(value, Quad):
        return value
    return Quad(Fraction(value))


def ground_patterns() -> tuple[dict[tuple[int, ...], list[tuple[int, np.ndarray]]], int]:
    n = len(A9)
    x_states = spins(n)
    child_norms = []
    for omitted in range(n):
        keep = [i for i in range(n) if i != omitted]
        energies = np.einsum("bi,ij,bj->b", x_states[:, keep],
                             A9[np.ix_(keep, keep)], x_states[:, keep])
        child_norms.append(int(np.max(np.abs(energies))))
    assert child_norms == [24] * n

    patterns: dict[tuple[int, ...], list[tuple[int, np.ndarray]]] = {}
    full_norm = 24
    for sigma in (-1, 1):
        for x in x_states:
            pattern = []
            for omitted in range(n):
                keep = [i for i in range(n) if i != omitted]
                energy = sigma * int(x[keep] @ A9[np.ix_(keep, keep)] @ x[keep])
                if energy == child_norms[omitted]:
                    pattern.append(omitted)
            if pattern:
                patterns.setdefault(tuple(pattern), []).append((sigma, x))
    assert len(patterns) == 22
    assert sum(map(len, patterns.values())) == 25
    return patterns, full_norm


def main() -> None:
    patterns, full_norm = ground_patterns()
    a = Quad(Fraction(21, 130), Fraction(-4, 130))
    b = Quad(Fraction(1, 20), Fraction(1, 20))
    c = Quad(Fraction(15, 52), Fraction(-1, 52))
    assert a.sign() > 0 and b.sign() > 0 and c.sign() > 0

    masses = {
        (3, 5, 7): a,
        (4, 7, 8): a,
        (0, 6, 7): b,
        (1, 2, 7): b,
        (0, 4, 8): c,
        (2, 3, 5): c,
    }
    assert set(masses) <= set(patterns)
    assert sum(masses.values(), Quad(Fraction(0))) == Quad(Fraction(1))

    # Every supported coverage pattern has a representative that is also a
    # full A9 ground, as claimed in the ledger.
    for pattern in masses:
        assert any(sigma * int(x @ A9 @ x) == full_norm
                   for sigma, x in patterns[pattern])

    z = [sum((mass for pattern, mass in masses.items() if i in pattern),
             Quad(Fraction(0))) for i in range(9)]
    u = Quad(Fraction(22, 65), Fraction(2, 65))
    v = Quad(Fraction(9, 20), Fraction(-1, 20))
    w = Quad(Fraction(11, 26), Fraction(1, 26))
    assert z == [u, b, u, v, v, v, b, w, v]

    pi = [Fraction(t, 25) for t in (4, 2, 4, 0, 4, 4, 2, 5, 0)]
    active = set(masses)
    for pattern in patterns:
        h = sum((Quad(pi[i]) / z[i] for i in pattern if pi[i]),
                Quad(Fraction(0)))
        gap = Quad(Fraction(1)) - h
        assert gap.sign() >= 0
        assert (gap.sign() == 0) == (pattern in active)

    information = -(8 * math.log(u.value()) + 4 * math.log(b.value())
                    + 8 * math.log(v.value()) + 5 * math.log(w.value())) / 25
    assert abs(information - 1.022686788870189) < 1e-12
    print({"patterns": len(patterns), "supported": len(masses),
           "information_nats": information})
    print("PASS: exact A9 selector-information KKT certificate")


if __name__ == "__main__":
    main()
