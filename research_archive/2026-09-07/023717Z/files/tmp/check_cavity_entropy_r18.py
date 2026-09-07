#!/usr/bin/env python3
"""Independent audit of the Wave-18 cavity/Hellinger identities.

The computation uses the full oriented cube (sigma,x).  This is a twofold
lift of the oriented-projective convention in the ledger, so all partition
sums are doubled and every ratio is unchanged.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, product
import math

import numpy as np
from scipy.special import logsumexp

from check_selected_child_soft_r17 import A8, A9


def states(n: int):
    for sigma in (-1, 1):
        for x in product((-1, 1), repeat=n):
            yield sigma, x


def energy(a: np.ndarray, sigma: int, x: tuple[int, ...]) -> int:
    return sigma * sum(
        2 * int(a[i, j]) * x[i] * x[j]
        for i in range(len(a)) for j in range(i + 1, len(a))
    )


def data(a: np.ndarray):
    records = [(sigma, x, energy(a, sigma, x)) for sigma, x in states(len(a))]
    q = max(e for _, _, e in records)
    deficits = Counter(q - e for _, _, e in records)
    return q, records, deficits


def child(a: np.ndarray, i: int) -> np.ndarray:
    ids = [j for j in range(len(a)) if j != i]
    return a[np.ix_(ids, ids)]


def log_deficit_sum(deficits: Counter[int], beta: float) -> float:
    terms = [math.log(count) - beta * deficit for deficit, count in deficits.items()]
    return float(logsumexp(terms))


def audit_one(a: np.ndarray, beta: float) -> tuple[int, list[tuple[int, float]]]:
    r = len(a)
    q, records, deficits = data(a)
    logd = log_deficit_sum(deficits, beta)
    by_state = {(sigma, x): e for sigma, x, e in records}
    out = []
    for i in range(r):
        c = child(a, i)
        qc, _, deficits_c = data(c)
        d = q - qc

        # Every unordered i-edge is represented by its endpoint x_i=+1.
        edge_average_deficits = Counter()
        edge_bc_terms = []
        for sigma, x, e in records:
            if x[i] != 1:
                continue
            xf = list(x)
            xf[i] = -1
            xf = tuple(xf)
            ef = by_state[(sigma, xf)]
            du, dv = q - e, q - ef
            assert (du + dv) % 2 == 0
            edge_average_deficits[(du + dv) // 2] += 1
            edge_bc_terms.append(-beta * (du + dv) / 2)

        shifted_child = Counter({d + delta: count for delta, count in deficits_c.items()})
        assert edge_average_deficits == shifted_child

        logdc = log_deficit_sum(deficits_c, beta)
        kappa_deficit = d + (logd - math.log(2) - logdc) / beta
        log_bc_edges = math.log(2) + float(logsumexp(edge_bc_terms)) - logd
        assert abs(log_bc_edges + beta * kappa_deficit) < 2e-12

        # Direct Gibbs Bhattacharyya coefficient over directed states.
        logz = logsumexp([beta * e for _, _, e in records])
        bc = 0.0
        probs = {
            (sigma, x): math.exp(beta * e - logz)
            for sigma, x, e in records
        }
        for sigma, x, _ in records:
            xf = list(x)
            xf[i] *= -1
            bc += math.sqrt(probs[(sigma, x)] * probs[(sigma, tuple(xf))])
        assert abs(math.log(bc) + beta * kappa_deficit) < 2e-12
        out.append((d, kappa_deficit))
    return q, out


def all_restrictions(name: str, a: np.ndarray) -> None:
    total = 0
    for r in range(2, len(a) + 1):
        for ids in combinations(range(len(a)), r):
            b = a[np.ix_(ids, ids)]
            for beta in (0.1, 0.5, 1.0):
                _, values = audit_one(b, beta)
                total += len(values)
    print(f"{name}: PASS {total} restriction-coordinate-temperature audits")


def hamming_ball_example(r: int, k: int, beta: float) -> tuple[float, float, float]:
    # Uniform on two disjoint antipodal Hamming balls of radius k.
    ball = sum(math.comb(r, j) for j in range(k + 1))
    internal_one_direction_one_ball = sum(math.comb(r - 1, j) for j in range(k))
    bc = 2 * internal_one_direction_one_ball / ball
    kappa = -math.log(bc) / beta
    entropy_per_coordinate = math.log(2 * ball) / r
    return bc, kappa, entropy_per_coordinate


def main() -> None:
    for name, a in (("A8", A8), ("A9", A9)):
        for beta in (0.1, 0.5, 1.0):
            q, values = audit_one(a, beta)
            print(name, "beta", beta, "Q", q,
                  "d", [d for d, _ in values],
                  "kappa", [round(k, 9) for _, k in values])
    all_restrictions("A8", A8)
    all_restrictions("A9", A9)

    print("antipodal Hamming-ball entropy-only obstruction, beta=1")
    for r in (16, 64, 256, 1024, 4096):
        k = math.isqrt(r)
        bc, kappa, entropy_rate = hamming_ball_example(r, k, 1.0)
        print(r, k, f"BC={bc:.12g}", f"kappa={kappa:.9f}",
              f"kappa/sqrt(r)={kappa/math.sqrt(r):.9f}",
              f"log|S|/r={entropy_rate:.9f}")


if __name__ == "__main__":
    main()
