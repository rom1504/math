#!/usr/bin/env python3
"""Finite checks for the orbit-query large-deviation theorem."""

from __future__ import annotations

import math


def binary_tail(n: int, a: float) -> int:
    # Sum of n Rademachers is at least a*n.
    lo = math.ceil((1 + a) * n / 2)
    return sum(math.comb(n, j) for j in range(lo, n + 1))


def kl(p: float, q: float) -> float:
    return p * math.log(p / q) + (1 - p) * math.log((1 - p) / (1 - q))


def main() -> None:
    for a in (0.2, 0.4, 0.6):
        rate = kl((1 + a) / 2, 0.5)
        previous = None
        for n in (20, 40, 80, 160):
            count = binary_tail(n, a)
            empirical = -(math.log(count) - n * math.log(2)) / n
            assert empirical <= rate + 0.2
            previous = empirical
        assert previous is not None and abs(previous - rate) < 0.05
        print(f"a={a:.1f}: empirical rate {previous:.6f}, Cramer rate {rate:.6f}")
    print("orbit-query large-deviation checks: PASS")


if __name__ == "__main__":
    main()
