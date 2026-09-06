#!/usr/bin/env python3
"""Exact Wave 51 audit of a cubic completion-tail lower bound.

All selectors, local oriented states, and outside completions are exhausted
for the stored small exact minimizers.  The order-ten matrix is one
deterministic exact-minimizer sample, not an order-ten classification.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def cubic_max(aa: float, bb: float, r: float) -> float:
    """max_(r<=u<=aa) (aa-u)(u-r)(u+bb), by its unique interior root."""
    disc = (aa - bb + r) ** 2 + 3 * (bb * (aa + r) - aa * r)
    # Roundoff can only make a theoretically nonnegative discriminant tiny negative.
    disc = max(0.0, disc)
    u = (aa - bb + r + math.sqrt(disc)) / 3
    u = min(aa, max(r, u))
    return (aa - u) * (u - r) * (u + bb)


def audit(a: np.ndarray, name: str, m: int, t: float = 0.0) -> dict:
    n = len(a)
    q = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    p32 = p**1.5

    count = 0
    actual_sum = 0.0
    second_sum = 0.0
    cubic_sum = 0.0
    bonami_sum = 0.0
    linear_dominance_sum = 0.0
    combined_sum = 0.0
    negative_actual_sum = 0.0
    negative_count = 0
    cubic_positive_count = 0
    cubic_only_count = 0
    linear_only_count = 0
    linear_far_count = 0
    max_mu3_error = 0.0
    max_numerator_error = 0.0
    max_fourth_ratio = 0.0
    min_actual_minus_bound = math.inf

    for s in itertools.combinations(range(n), m):
        s = tuple(s)
        outside = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        outside_a = a[np.ix_(outside, outside)]
        qs = qnorm(child)
        completion_states = list(itertools.product((-1, 1), repeat=n - m))
        for sigma in (-1, 1):
            for y in projective_spins(m):
                e = sigma * int(y @ child @ y)
                h = qs - p32 * q - t
                g = (1 - p2) * e - h
                if g > 1e-12:
                    continue
                r = -g / p2
                aa = q + e
                bb = q - e
                if not (0 <= r < aa and aa > 0):
                    continue

                cross = a[np.ix_(outside, s)] @ y
                bvec = 2 * sigma * cross
                # The oriented quadratic coefficient matrix is sigma*A[T].
                bmat = sigma * outside_a
                v_formula = 4 * int(cross @ cross) + 2 * (n - m) * (n - m - 1)
                mu3_formula = 6 * int(bvec @ bmat @ bvec) + 8 * int(
                    np.trace(bmat @ bmat @ bmat)
                )

                zvals = []
                for wt in completion_states:
                    w = np.asarray(wt, dtype=np.int64)
                    zvals.append(int(bvec @ w + w @ bmat @ w))
                zvals = np.asarray(zvals, dtype=np.float64)
                actual = float(np.mean(zvals <= -r + 1e-12))
                v = float(np.mean(zvals**2))
                mu3 = float(np.mean(zvals**3))
                mu4 = float(np.mean(zvals**4))
                if v > 0:
                    max_fourth_ratio = max(max_fourth_ratio, mu4 / (v * v))
                max_mu3_error = max(max_mu3_error, abs(mu3 - mu3_formula))
                assert abs(v - v_formula) < 1e-12

                d = q * q - e * e - v
                n3 = mu3 + 2 * e * v - r * d
                poly = (zvals + aa) * (zvals + r) * (zvals - bb)
                max_numerator_error = max(max_numerator_error, abs(float(np.mean(poly)) - n3))
                m3 = cubic_max(aa, bb, r)
                cubic = max(0.0, n3) / m3 if m3 > 0 else 0.0
                second = max(0.0, v - r * bb) / ((aa - r) * (aa + bb))
                bonami = max(0.0, 1 / 18 - r / math.sqrt(v)) ** 2 if v > 0 else 0.0
                linvar = float(bvec @ bvec)
                quadvar = float(2 * (n - m) * (n - m - 1))
                linear_dominance = 0.0
                if linvar > 0:
                    linsd = math.sqrt(linvar)
                    # Each theta in this finite rational grid gives a rigorous
                    # Paley-Zygmund/Cantelli bound; maximization preserves it.
                    for j in range(1, 1000):
                        theta = j / 1000
                        vcut = theta * linsd - r
                        if vcut <= 0:
                            continue
                        candidate = 0.5 * max(
                            0.0,
                            (1 - theta * theta) ** 2 / 3
                            - quadvar / (quadvar + vcut * vcut),
                        )
                        linear_dominance = max(linear_dominance, candidate)
                combined = max(second, cubic, bonami, linear_dominance)

                min_actual_minus_bound = min(min_actual_minus_bound, actual - combined)
                count += 1
                negative_count += 1
                actual_sum += actual
                negative_actual_sum += actual
                second_sum += second
                cubic_sum += cubic
                bonami_sum += bonami
                linear_dominance_sum += linear_dominance
                combined_sum += combined
                cubic_positive_count += int(cubic > 1e-15)
                cubic_only_count += int(cubic > 1e-15 and second <= 1e-15)
                linear_only_count += int(linear_dominance > 1e-15 and second <= 1e-15)
                linear_far_count += int(linear_dominance > 1e-15 and r * bb >= v - 1e-12)

    out = {
        "name": name,
        "n": n,
        "m": m,
        "q": q,
        "states": count,
        "actual": actual_sum / count,
        "second": second_sum / count,
        "cubic": cubic_sum / count,
        "bonami": bonami_sum / count,
        "linear_dominance": linear_dominance_sum / count,
        "combined": combined_sum / count,
        "cubic_positive_fraction": cubic_positive_count / count,
        "cubic_only_fraction": cubic_only_count / count,
        "linear_only_fraction": linear_only_count / count,
        "linear_far_fraction": linear_far_count / count,
        "min_actual_minus_bound": min_actual_minus_bound,
        "max_mu3_error": max_mu3_error,
        "max_numerator_error": max_numerator_error,
        "max_fourth_ratio": max_fourth_ratio,
    }
    print(
        f"{name:10s} n={n:2d} m={m:2d} q={q:2d} states={count:6d} "
        f"actual={out['actual']:.8f} H2={out['second']:.8f} "
        f"H3={out['cubic']:.8f} HB={out['bonami']:.8f} "
        f"HLQ={out['linear_dominance']:.8f} "
        f"Hcomb={out['combined']:.8f} mu4/V2<={out['max_fourth_ratio']:.4f} "
        f"H3>0={out['cubic_positive_fraction']:.4f} "
        f"H3-only={out['cubic_only_fraction']:.4f} "
        f"HLQ-only={out['linear_only_fraction']:.4f} "
        f"HLQ-far={out['linear_far_fraction']:.4f}"
    )
    assert count > 0
    assert min_actual_minus_bound >= -1e-12
    assert max_mu3_error < 1e-12
    assert max_numerator_error < 1e-9
    assert max_fourth_ratio <= 81 + 1e-12
    return out


def pairing_audit() -> None:
    """Verify the exact w <-> -w projective-pair CDF formula."""
    a = sample_order_ten(0)
    n, m = len(a), 6
    q = qnorm(a)
    checks = 0
    for s in itertools.combinations(range(n), m):
        s = tuple(s)
        outside = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        bmat = a[np.ix_(outside, outside)]
        for y in projective_spins(m):
            bvec = 2 * (a[np.ix_(outside, s)] @ y)
            # Several integer thresholds, including negative ones.
            for h in (-q, -q // 2, 0, q // 2, q):
                vals = []
                formula = []
                for wt in projective_spins(n - m):
                    w = np.asarray(wt, dtype=np.int64)
                    qv = int(w @ bmat @ w)
                    lv = int(bvec @ w)
                    vals.extend((qv + lv, qv - lv))
                    formula.append((int(qv - abs(lv) <= h) + int(qv + abs(lv) <= h)) / 2)
                assert abs(np.mean(np.asarray(vals) <= h) - np.mean(formula)) < 1e-12
                checks += 1
    print(f"pairing identities checked: {checks}")


def main() -> None:
    results = []
    for a, name, ms in (
        (A6, "A6", (3, 4, 5)),
        (A8, "A8", (4, 5, 6, 7)),
        (A9, "A9", (5, 6, 7, 8)),
        (sample_order_ten(0), "A10_seed0", (5, 6, 7, 8, 9)),
    ):
        for m in ms:
            results.append(audit(a, name, m))
    pairing_audit()
    assert any(x["combined"] > x["second"] + 1e-12 for x in results)
    print("PASS cubic completion-tail audit")


if __name__ == "__main__":
    main()
