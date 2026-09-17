#!/usr/bin/env python3
"""Numerical audit of the Wave-56 mesoscopic common-core asymptotics.

All calculations are in logs.  For each (n,p,ell), the script computes the
exact Johnson shell law, the exact common-core transition law, the smallest
number of top-intersection partners carrying lambda_2-h mass, and the best
one-threshold rectangle B(s) K_ell(s).
"""

from __future__ import annotations

import math


def lbinom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def logadd(x: float, y: float) -> float:
    if x == -math.inf:
        return y
    if y == -math.inf:
        return x
    if x < y:
        x, y = y, x
    return x + math.log1p(math.exp(y - x))


def entropy(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log(x) - (1 - x) * math.log(1 - x)


def audit(n: int, p0: float, ell: int) -> None:
    m = round(p0 * n)
    p = m / n
    lo = max(0, 2 * m - n)
    logN = lbinom(n, m)
    logd = lbinom(m, ell) + lbinom(n - ell, m - ell)

    lam1 = ell * (n - m) / (m * (n - ell))
    lam2 = (
        ell * (ell - 1) * (n - m) * (n - m - 1)
        / (m * (m - 1) * (n - ell) * (n - ell - 1))
    )
    h = math.exp(-lbinom(n - ell, m - ell))
    target = lam2 - h

    # log A_j and log K(j), with j=m (the self-loop) omitted below.
    rows: list[tuple[int, float, float, float]] = []
    for j in range(lo, m + 1):
        logA = lbinom(m, j) + lbinom(n - m, m - j)
        logK = lbinom(j, ell) - logd if j >= ell else -math.inf
        rows.append((j, logA, logK, logA + logK))

    # Descending shell tail.  weighted_tail is exact K mass, count_tail is
    # exact number of partners.  Locate first crossing of lambda_2-h.
    log_count_tail = -math.inf
    log_weight_tail = -math.inf
    crossing = None
    best_rect = (-math.inf, None)
    for j, logA, logK, logW in reversed(rows[:-1]):
        log_count_tail = logadd(log_count_tail, logA)
        log_weight_tail = logadd(log_weight_tail, logW)
        log_rect = log_count_tail + logK
        if log_rect > best_rect[0]:
            best_rect = (log_rect, j)
        if crossing is None and log_weight_tail >= math.log(target):
            crossing = (j, log_count_tail, log_weight_tail)

    assert crossing is not None
    jcap, logr, logw = crossing
    a = (1 - p) / p
    predicted_cost = 0.5 * a * a * ell * ell / n
    alpha = ell / n
    beta_star = alpha + (p - alpha) ** 2 / (1 - alpha)
    vstar = p * entropy((p - beta_star) / p) + (1 - p) * entropy(
        (p - beta_star) / (1 - p)
    )
    entropy_cost = n * (entropy(p) - vstar)
    exact_cost = logN - logr
    predicted_rect = p / (math.sqrt(2 * math.pi) * (1 - p)) * math.sqrt(n) / ell
    ratio_constant = math.exp(best_rect[0]) / predicted_rect
    mu = ell + (m - ell) ** 2 / (n - ell)
    var = (m - ell) ** 2 * (n - m) ** 2 / ((n - ell) ** 2 * (n - ell - 1))
    D = (1 - lam1) + n * (lam1 - lam2)

    print(
        f"n={n:7d} ell={ell:7d} gamma={math.log(ell)/math.log(n):.4f} "
        f"lambda2={lam2:.3e} D/(a ell)={D/(a*ell):.5f}"
    )
    print(
        f"  mu/n={mu/n:.8f} beta*={beta_star:.8f} var/n={var/n:.8f} "
        f"cap-shell={jcap}"
    )
    print(
        f"  log(N/r_cap)={exact_cost:.5f} quad={predicted_cost:.5f} "
        f"I(beta*)={entropy_cost:.5f} exact/I={exact_cost/entropy_cost:.5f} "
        f"captured={math.exp(logw):.3e}"
    )
    print(
        f"  max rectangle={math.exp(best_rect[0]):.3e} at j={best_rect[1]}, "
        f"pred={predicted_rect:.3e}, ratio={ratio_constant:.5f}, "
        f"rectangle/lambda2={math.exp(best_rect[0])/lam2:.3e}"
    )


if __name__ == "__main__":
    # Large enough to see the mesoscopic constants while keeping O(n) shell
    # enumeration cheap.  The exponents straddle 5/6.
    for exponent in (0.76, 0.80, 5 / 6, 0.87):
        n0 = 200_000
        audit(n0, 0.6, max(2, round(n0**exponent)))
