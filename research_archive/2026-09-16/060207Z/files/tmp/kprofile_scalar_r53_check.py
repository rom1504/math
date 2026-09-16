#!/usr/bin/env python3
"""Deterministic checks for the Wave 53 scalar K-profile memo.

This verifies the exact finite-dimensional K-functional formula and the
energy/max-coordinate lower bound, then searches stored exact minimizers for
pointwise profile failures.  The finite examples are algebraic obstructions,
not asymptotic evidence about saved populations.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins
from far_tail_kfunctional_r52_check import exact_kfunctional, qnorm


def clipped_formula(beta: np.ndarray, h: int) -> float:
    """Supremum of the explicit clipped dual family."""

    b = np.abs(np.asarray(beta, dtype=float))
    if not np.any(b):
        return 0.0
    # The supremum is attained at N(a)=h unless the l2 constraint is inactive.
    support = int(np.count_nonzero(b))
    if h >= support:
        return float(np.sum(b))
    lo = 0.0
    hi = float(np.max(b))
    while float(np.sum(np.minimum(1.0, b / hi) ** 2)) > h:
        hi *= 2.0
    for _ in range(120):
        mid = (lo + hi) / 2.0
        n_a = float(np.sum(np.minimum(1.0, b / mid) ** 2))
        if n_a > h:
            lo = mid
        else:
            hi = mid
    a = hi
    z = np.minimum(1.0, b / a)
    return float(b @ z)


def audit_random_vectors() -> None:
    rng = np.random.default_rng(5301)
    checked = 0
    sharp_cases = 0
    for k in range(1, 18):
        for h in range(1, k + 2):
            for _ in range(30):
                beta = rng.integers(-9, 10, size=k).astype(float)
                kval = exact_kfunctional(beta, h)
                kclip = clipped_formula(beta, h)
                assert abs(kval - kclip) < 1e-9
                energy = float(beta @ beta)
                maximum = float(np.max(np.abs(beta)))
                if maximum == 0:
                    lower = 0.0
                else:
                    lower = min(energy / maximum, math.sqrt(h * energy))
                assert kval + 1e-9 >= lower

                # Prefix/head plus diffuse-tail dual certificates, including
                # all cases in which the proportional tail is box-feasible.
                b = np.sort(np.abs(beta))[::-1]
                for s in range(min(h, k) + 1):
                    tail = b[s:]
                    tail_energy = float(tail @ tail)
                    budget = h - s
                    if budget == 0:
                        certificate = float(np.sum(b[:s]))
                    elif tail_energy == 0:
                        certificate = float(np.sum(b[:s]))
                    elif float(np.max(tail)) ** 2 * budget <= tail_energy + 1e-12:
                        certificate = float(np.sum(b[:s])) + math.sqrt(
                            budget * tail_energy
                        )
                    else:
                        continue
                    assert kval + 1e-9 >= certificate
                checked += 1

            # Flat vectors show sharpness in both regimes whenever possible.
            for support in range(1, k + 1):
                beta = np.zeros(k)
                beta[:support] = 7.0
                kval = exact_kfunctional(beta, h)
                energy = float(beta @ beta)
                maximum = 7.0
                lower = min(energy / maximum, math.sqrt(h * energy))
                assert abs(kval - lower) < 1e-9
                sharp_cases += 1
    print(f"random/exact formula checks={checked}; flat sharpness checks={sharp_cases}")


def best_failure(a: np.ndarray, name: str, m: int) -> dict[str, object] | None:
    """Find a g<0 state with high-energy necessary condition but K<2r."""

    n = len(a)
    q = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    hfun = max(1, min(3, n - m - 1))
    best: dict[str, object] | None = None
    negative_states = 0
    energy_possible_states = 0
    profile_failures = 0

    for s0 in itertools.combinations(range(n), m):
        s = tuple(s0)
        tset = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        cross = a[np.ix_(tset, s)]
        qchild = qnorm(child)
        hmargin = qchild - p ** 1.5 * q  # t=0 in ledger notation
        for sigma in (-1, 1):
            for y0 in projective_spins(m):
                y = np.asarray(y0, dtype=np.int64)
                e = sigma * int(y @ child @ y)
                g = (1 - p2) * e - hmargin
                if g >= 0:
                    continue
                negative_states += 1
                r = -g / p2
                u = cross @ y
                beta = 2 * sigma * u
                energy = float(beta @ beta)
                maximum = float(np.max(np.abs(beta)))
                kval = exact_kfunctional(beta, hfun)
                # Necessary energy condition for K>=2r, but actual K fails.
                if energy + 1e-9 < 4 * r * r / hfun:
                    continue
                energy_possible_states += 1
                if kval + 1e-9 >= 2 * r:
                    continue
                profile_failures += 1
                ratio = kval / (2 * r)
                candidate = {
                    "name": name,
                    "n": n,
                    "m": m,
                    "q": q,
                    "H": hfun,
                    "S": s,
                    "sigma": sigma,
                    "y": tuple(int(x) for x in y),
                    "Qchild": qchild,
                    "e": e,
                    "hmargin": hmargin,
                    "r": r,
                    "u": tuple(int(x) for x in u),
                    "Ebeta": energy,
                    "Mbeta": maximum,
                    "K": kval,
                    "ratio": ratio,
                    "energy_ratio": energy * hfun / (4 * r * r),
                }
                if best is None or ratio < float(best["ratio"]):
                    best = candidate
    if best is not None:
        best["negative_states"] = negative_states
        best["energy_possible_states"] = energy_possible_states
        best["profile_failures"] = profile_failures
    return best


def audit_stored_minimizers() -> None:
    found = 0
    for a, name, m in (
        (A6, "A6", 3),
        (A8, "A8", 4),
        (A9, "A9", 5),
        (sample_order_ten(0), "A10_seed0", 5),
    ):
        witness = best_failure(a, name, m)
        if witness is None:
            print(f"{name}: no requested witness")
            continue
        found += 1
        print(
            f"{name}: n={witness['n']} m={witness['m']} q={witness['q']} "
            f"H={witness['H']} S={witness['S']} sigma={witness['sigma']} "
            f"y={witness['y']} Q_S={witness['Qchild']} e={witness['e']} "
            f"r={witness['r']:.12f} u={witness['u']} "
            f"Ebeta={witness['Ebeta']:.0f} Mbeta={witness['Mbeta']:.0f} "
            f"K={witness['K']:.12f} K/(2r)={witness['ratio']:.9f} "
            f"HE/(4r^2)={witness['energy_ratio']:.9f} "
            f"counts(negative,energy-possible,fail)="
            f"({witness['negative_states']},{witness['energy_possible_states']},"
            f"{witness['profile_failures']})"
        )
        if name == "A8":
            exact_r = (124 - 70 * math.sqrt(2)) / 3
            assert witness["S"] == (0, 1, 2, 5)
            assert witness["u"] == (0, 0, 4, 4)
            assert witness["Qchild"] == 12 and witness["e"] == 4
            assert abs(float(witness["r"]) - exact_r) < 1e-12
            # For u=(0,0,4,4), support 2 <= H=3, hence K_H(u)=8.
            # The energy-only requirement passes, but the max-coordinate
            # requirement fails: sqrt(HG)>r>G/||u||_infty=8.
            assert math.sqrt(3 * 32) > exact_r > 8
            assert abs(float(witness["K"]) - 16) < 1e-12
    assert found >= 1


def exponent_model_audit() -> None:
    """Check exponent compatibility of the abstract concentrated model."""

    # d=n^(1/2-c+eta), H=n^(3/4-c), r~n*d.  The choices below make
    # r/T -> infinity, d/H -> 0, and r/q -> 0 while cross energy is above
    # the energy-only threshold.  This is only an exponent/model audit.
    for c in (0.02, 0.10, 0.20, 0.249):
        eta = min(c, 0.25) / 2
        d_exp = 0.5 - c + eta
        h_exp = 0.75 - c
        r_exp = 1 + d_exp
        target_exp = 1.5 - c
        q_exp = 1.5
        gext_exp = 2 + d_exp
        energy_need_exp = 2 * r_exp - h_exp
        assert r_exp > target_exp
        assert d_exp < h_exp
        assert r_exp < q_exp
        assert gext_exp > energy_need_exp
    print("abstract concentrated-model exponent checks passed")


def main() -> None:
    audit_random_vectors()
    audit_stored_minimizers()
    exponent_model_audit()
    print("PASS Wave 53 scalar K-profile audit")


if __name__ == "__main__":
    main()
