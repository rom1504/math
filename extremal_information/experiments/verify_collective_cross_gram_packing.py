#!/usr/bin/env python3
"""Finite checks for collective_cross_gram_packing_and_response_modulus.md."""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np

from verify_bcx_two_port_holonomy import regular_hadamard


def signs(p: int) -> list[np.ndarray]:
    return [np.asarray(x, dtype=np.int64) for x in product((-1, 1), repeat=p)]


def projective_hamming(s: np.ndarray, t: np.ndarray) -> int:
    h = int(np.count_nonzero(s != t))
    return min(h, len(s) - h)


def collective_distance(
    g: np.ndarray, r: np.ndarray, gp: np.ndarray, rp: np.ndarray
) -> float:
    p = g.shape[0]
    best = 0.0
    for epsilon in signs(p):
        for sigma in (-1, 1):
            delta = (g - gp) + sigma * (r - rp)
            best = max(best, abs(float(epsilon @ delta @ epsilon)) / (p * p))
    return best


def greedy_projective_code(p: int, radius: int) -> list[np.ndarray]:
    # First coordinate +1 chooses one representative from each projective class.
    candidates = [s for s in signs(p) if s[0] == 1]
    code: list[np.ndarray] = []
    for s in candidates:
        if all(projective_hamming(s, t) >= radius for t in code):
            code.append(s)
    return code


def verify_rank_one_packing() -> None:
    p, eta = 8, 0.25
    radius = int(eta * p)
    code = greedy_projective_code(p, radius)
    assert len(code) > 1

    q, h, _ = regular_hadamard(2)
    n = h.shape[0]
    w = np.ones(n, dtype=np.int64)
    assert np.array_equal(h @ w, q * w)

    pairs = []
    for s in code:
        ports = np.stack([int(si) * w for si in s])
        g = ports @ ports.T / n
        r = ports @ h @ ports.T / (q * n)
        assert np.array_equal(g, np.outer(s, s))
        assert np.array_equal(r, g)
        assert np.array_equal(np.diag(g), np.ones(p))
        pairs.append((s, g, r))

    minimum = float("inf")
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            s, g, r = pairs[i]
            t, gp, rp = pairs[j]
            hamming = int(np.count_nonzero(s != t))
            expected = 8 * hamming * (p - hamming) / (p * p)
            observed = collective_distance(g, r, gp, rp)
            assert abs(observed - expected) < 1e-12
            assert observed >= 8 * eta * (1 - eta) - 1e-12
            minimum = min(minimum, observed)

    print(
        "rank-one projective packing:",
        f"p={p}, radius={radius}, size={len(code)}, min_dq={minimum:.6f}",
    )


def psi(kappa: float, a: float, b: float) -> float:
    """Convex one-variable trust response, including its t=0 hard limit."""

    if a == 0.0 and kappa * kappa * b <= 4.0:
        return 0.5 + kappa * kappa * b / 4.0

    def derivative(t: float) -> float:
        dangerous = 0.0 if a == 0.0 else kappa * kappa * a / (2 * t * t)
        safe = kappa * kappa * b / (2 * (t + 2) ** 2)
        return 0.5 - dangerous - safe

    lo = 0.0
    hi = max(1.0, kappa * sqrt(max(a, 1e-300)))
    while derivative(hi) < 0:
        hi *= 2
    for _ in range(160):
        mid = (lo + hi) / 2
        if derivative(mid) < 0:
            lo = mid
        else:
            hi = mid
    t = hi
    return 0.5 + t / 2 + kappa * kappa * (a / t + b / (t + 2)) / 2


def verify_response_modulus() -> None:
    rng = np.random.default_rng(20260817)
    for kappa in (0.4, 1.0, 1.7):
        for delta in (1e-4, 0.01, 0.2):
            e = delta / 2
            bound = kappa * sqrt(e) + kappa * kappa * e / 4
            for _ in range(500):
                a, b = rng.random(2)
                ap = float(np.clip(a + rng.uniform(-e, e), 0, 1))
                bp = float(np.clip(b + rng.uniform(-e, e), 0, 1))
                gap = abs(psi(kappa, a, b) - psi(kappa, ap, bp))
                assert gap <= bound + 2e-12
    print("uniform trust-response modulus: PASS")


def verify_hard_case_and_contextual_gap() -> None:
    kappa = 1.0
    for c in (1e-6, 1e-4, 0.01, 0.1):
        b = (1 - c) / 2
        gap = psi(kappa, c, b) - psi(kappa, 0.0, b)
        lower = kappa * sqrt(c * (1 - kappa * kappa * b / 4))
        upper = kappa * sqrt(c)
        assert lower <= gap + 2e-12
        assert gap <= upper + 2e-12

    # A rank-one Hamming pair has an exact directed response-table gap.
    p = 8
    s = np.ones(p, dtype=np.int64)
    t = s.copy()
    t[:2] = -1
    a_s = float((s @ s) ** 2 / p**2)
    a_t = float((s @ t) ** 2 / p**2)
    observed = abs(psi(kappa, a_s, 0) - psi(kappa, a_t, 0))
    assert abs(observed - 2 * kappa * 2 / p) < 1e-12
    print(
        "hard-case and contextual response:",
        f"sqrt_ratio={gap/sqrt(c):.6f}, rank_one_gap={observed:.6f}",
    )


def main() -> None:
    verify_rank_one_packing()
    verify_response_modulus()
    verify_hard_case_and_contextual_gap()
    print("collective cross-Gram packing checks: PASS")


if __name__ == "__main__":
    main()
