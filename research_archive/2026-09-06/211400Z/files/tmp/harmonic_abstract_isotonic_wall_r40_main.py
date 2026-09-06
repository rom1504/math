#!/usr/bin/env python3
"""Numerical audit of the two-bit abstract isotonic wall."""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.legendre import leggauss


def binary_kl(p: float) -> float:
    q = 1.0 - p
    ans = 0.0
    if p > 0.0:
        ans += p * math.log(2.0 * p)
    if q > 0.0:
        ans += q * math.log(2.0 * q)
    return ans


def metrics(L: float, order: int = 400) -> tuple[float, float, float]:
    nodes, weights = leggauss(order)
    K = 0.0
    # Integrate in u=sL, where the mass lives on a fixed scale.
    for node, weight in zip(nodes, weights):
        u = 0.5 * L * (float(node) + 1.0)
        du = 0.5 * L * float(weight)
        s = u / L
        eu = math.exp(u) if u < 700.0 else math.inf
        if math.isinf(eu):
            h = 0.0
            a = 0.0
        else:
            h = eu / (1.0 + eu) ** 2
            a = 1.0 / (1.0 + 3.0 * eu)
        a1 = 1.0 / (1.0 + 3.0 * math.exp(L)) if L < 700.0 else 0.0
        w = 0.5 * (a + a1)
        v = L * L * h
        K += du / L * (4.0 / 3.0) * v * math.sqrt(
            max(0.0, s * (1.0 - s) * w * (1.0 - w))
        )
    p = 1.0 / (1.0 + math.exp(-L))
    context_mass = (1.0 + math.exp(L)) / (1.0 + 3.0 * math.exp(L))
    endpoint = 2.0 * context_mass * binary_kl(p)
    z = 0.25 * (1.0 + 3.0 * math.exp(L))
    mu_high = 3.0 * math.exp(L) / (1.0 + 3.0 * math.exp(L))
    entropy = mu_high * L - math.log(z)
    return K * K, endpoint, entropy


def main() -> None:
    ratios = []
    final = None
    for L in (4.0, 8.0, 16.0, 32.0, 64.0):
        K2, endpoint, entropy = metrics(L)
        final = (K2, endpoint, entropy)
        ratios.append(K2 / endpoint)
        print({"L": L, "K2": K2, "K2/L": K2 / L,
               "endpoint": endpoint, "entropy": entropy,
               "ratio": K2 / endpoint})
    assert ratios[-1] > 5.0 * ratios[1]
    assert final is not None
    assert abs(final[1] - 2.0 * math.log(2.0) / 3.0) < 1e-12
    assert abs(final[2] - math.log(4.0 / 3.0)) < 1e-12
    print("PASS abstract harmonic isotonic wall")


if __name__ == "__main__":
    main()
