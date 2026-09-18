#!/usr/bin/env python3
"""Wave 50 audit of the completion-CDF/Cantelli reduction.

All states, completions, and selectors are exhausted inside each supplied
matrix.  The order-ten matrices are deterministic exact-minimizer samples,
not an exhaustive classification of order ten.
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


def audit(a: np.ndarray, name: str, m: int, t: float = 0.0) -> dict:
    n = len(a)
    q = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    p32 = p ** 1.5

    actual_sum = 0.0
    cantelli_sum = 0.0
    two_branch_sum = 0.0
    support_nonpositive_sum = 0.0
    positive_g_mass = 0.0
    strong_local_mass = 0.0
    favorable_from_positive_g = 0.0
    favorable_from_nonpositive_g = 0.0
    conditional_checks = 0
    min_actual_minus_bound = math.inf
    max_variance_identity_error = 0.0
    max_parseval_excess = -math.inf
    min_abs_g = math.inf
    selectors = list(itertools.combinations(range(n), m))

    for s in selectors:
        s = tuple(s)
        outside = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        qs = qnorm(child)
        y_states = list(projective_spins(m))
        completion_states = list(itertools.product((-1, 1), repeat=n - m))
        for sigma in (-1, 1):
            for y in y_states:
                raw_local = int(y @ child @ y)
                e = sigma * raw_local
                y_gap = qs - p32 * q
                g = (1.0 - p2) * e - y_gap + t
                min_abs_g = min(min_abs_g, abs(g))

                cross = a[np.ix_(outside, s)] @ y
                variance_formula = 4 * int(cross @ cross) + 2 * (n - m) * (n - m - 1)
                full_energies = []
                for w_tuple in completion_states:
                    x = np.empty(n, dtype=np.int64)
                    x[list(s)] = y
                    x[list(outside)] = np.asarray(w_tuple, dtype=np.int64)
                    full_energies.append(sigma * int(x @ a @ x))
                full_energies = np.asarray(full_energies, dtype=np.int64)
                mean_energy = float(np.mean(full_energies))
                variance_enum = float(np.mean((full_energies - mean_energy) ** 2))
                max_variance_identity_error = max(
                    max_variance_identity_error,
                    abs(variance_enum - variance_formula),
                )
                max_parseval_excess = max(
                    max_parseval_excess,
                    e * e + variance_formula - q * q,
                )
                assert abs(mean_energy - e) < 1e-12

                delta = qs - e
                b_allow = (p32 - p2) * q
                if delta <= b_allow + t + 1e-12:
                    strong_local_mass += 1.0
                threshold = q + (b_allow + t - delta) / p2
                actual = float(np.mean(full_energies <= threshold + 1e-12))
                if g > 0:
                    bound = g * g / (g * g + p2 * p2 * variance_formula)
                    support_bound = 0.0
                    positive_g_mass += 1.0
                    favorable_from_positive_g += actual
                else:
                    bound = 0.0
                    # Here favorable means Z <= -r for the centered completion
                    # increment Z=E_full-e.  Its support is [-aa,bb].
                    r = -g / p2
                    aa = q + e
                    bb = q - e
                    if 0 <= r < aa and aa > 0:
                        support_bound = max(0.0, variance_formula - r * bb) / (
                            (aa - r) * (aa + bb)
                        )
                        # Equivalent near-cap form.  The Parseval slack is
                        # nonnegative because every full energy is in [-q,q].
                        d_slack = q * q - e * e - variance_formula
                        h = y_gap - t
                        denominator = 2 * q * (p2 * q + e - h)
                        simplified = max(
                            0.0,
                            (q - e) * (p2 * q + e - h) - p2 * d_slack,
                        ) / denominator
                        assert abs(support_bound - simplified) < 1e-12
                    else:
                        support_bound = 0.0
                    support_nonpositive_sum += support_bound
                    favorable_from_nonpositive_g += actual
                two_branch = max(bound, support_bound)
                min_actual_minus_bound = min(min_actual_minus_bound, actual - bound)
                min_actual_minus_bound = min(min_actual_minus_bound, actual - two_branch)
                actual_sum += actual
                cantelli_sum += bound
                two_branch_sum += two_branch
                conditional_checks += 1

    scale = float(conditional_checks)
    result = {
        "name": name,
        "n": n,
        "m": m,
        "q": q,
        "Z0": actual_sum / scale,
        "cantelli": cantelli_sum / scale,
        "cantelli_over_Z0": cantelli_sum / actual_sum,
        "two_branch": two_branch_sum / scale,
        "two_branch_over_Z0": two_branch_sum / actual_sum,
        "support_nonpositive": support_nonpositive_sum / scale,
        "positive_g_mass": positive_g_mass / scale,
        "strong_local_mass": strong_local_mass / scale,
        "Z0_share_positive_g": favorable_from_positive_g / actual_sum,
        "Z0_share_nonpositive_g": favorable_from_nonpositive_g / actual_sum,
        "min_actual_minus_bound": min_actual_minus_bound,
        "max_variance_identity_error": max_variance_identity_error,
        "max_parseval_excess": max_parseval_excess,
        "min_abs_g": min_abs_g,
    }
    print(
        f"{name:12s} n={n:2d} m={m:2d} q={q:2d} "
        f"Z0={result['Z0']:.9f} H+={result['cantelli']:.9f} "
        f"H2={result['two_branch']:.9f} H2/Z={result['two_branch_over_Z0']:.4f} "
        f"mass(g>0)={result['positive_g_mass']:.6f} "
        f"old={result['strong_local_mass']:.6f} "
        f"Zshare(g<=0)={result['Z0_share_nonpositive_g']:.4f}"
    )
    assert min_actual_minus_bound >= -1e-12
    assert max_variance_identity_error < 1e-12
    assert max_parseval_excess <= 0
    assert abs(result["Z0_share_positive_g"] + result["Z0_share_nonpositive_g"] - 1) < 1e-12
    return result


def main() -> None:
    results = []
    for a, name, ms in (
        (A6, "A6", (3, 4, 5)),
        (A8, "A8", (4, 5, 6, 7)),
        (A9, "A9", (5, 6, 7, 8)),
        (sample_order_ten(0), "A10_seed_0", (5, 6, 7, 8, 9)),
    ):
        for m in ms:
            results.append(audit(a, name, m))
    assert all(r["cantelli"] > 0 for r in results)
    print("PASS completion-CDF Cantelli audit")


if __name__ == "__main__":
    main()
