#!/usr/bin/env python3
"""Small-order checks for the three pre-ledger candidates.

This deliberately gauge-fixes the star at vertex 0 to +1, which is valid
because switching does not change either the hard norm or the spin pressure.
It exhausts every remaining signing through order 7.
"""

from itertools import combinations
import math
import numpy as np


def spin_table(n):
    # Quotient the irrelevant global spin flip by fixing x_0=+1.
    tails = ((np.arange(1 << (n - 1))[:, None]
              >> np.arange(n - 1)) & 1)
    x = np.ones((1 << (n - 1), n), dtype=np.int8)
    x[:, 1:] = 1 - 2 * tails
    return x


def edges(n):
    return list(combinations(range(n), 2))


def energy_table(n):
    x = spin_table(n)
    ee = edges(n)
    return np.stack([x[:, i] * x[:, j] for i, j in ee], axis=1).astype(np.int16)


def gauge_fixed_signings(n, batch=4096):
    ee = edges(n)
    free = [k for k, (i, j) in enumerate(ee) if i != 0]
    total = 1 << len(free)
    for lo in range(0, total, batch):
        ids = np.arange(lo, min(lo + batch, total), dtype=np.uint64)
        aa = np.ones((len(ids), len(ee)), dtype=np.int8)
        bits = ((ids[:, None] >> np.arange(len(free), dtype=np.uint64)) & 1)
        aa[:, free] = 1 - 2 * bits.astype(np.int8)
        yield aa


def logmeanexp(rows):
    maxima = rows.max(axis=1)
    return maxima + np.log(np.exp(rows - maxima[:, None]).mean(axis=1))


def exact_search(n, betas):
    prod = energy_table(n)
    best_m = math.inf
    best_a = []
    best_p = {beta: math.inf for beta in betas}
    best_p_a = {beta: None for beta in betas}
    for aa in gauge_fixed_signings(n):
        vals = aa.astype(np.int16) @ prod.T
        hard = np.abs(vals).max(axis=1)
        local = hard.min()
        if local < best_m:
            best_m = int(local)
            best_a = [a.copy() for a in aa[hard == local]]
        elif local == best_m:
            best_a.extend(a.copy() for a in aa[hard == local])
        for beta in betas:
            pressure = logmeanexp(beta * np.abs(vals) / math.sqrt(n)) / n
            k = int(np.argmin(pressure))
            if pressure[k] < best_p[beta]:
                best_p[beta] = float(pressure[k])
                best_p_a[beta] = aa[k].copy()
    return best_m, best_a, best_p, best_p_a


def induced_norm(a, n, subset):
    ee = edges(n)
    amap = {edge: int(a[k]) for k, edge in enumerate(ee)}
    subset = tuple(subset)
    b = np.array([amap[tuple(sorted((subset[i], subset[j])))]
                  for i, j in edges(len(subset))], dtype=np.int16)
    return int(np.max(np.abs(energy_table(len(subset)) @ b)))


def main():
    betas = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    results = {}
    for n in range(2, 8):
        results[n] = exact_search(n, betas)
        hard, optimizers, pressure, _ = results[n]
        print(f"n={n} M={hard} hard_optimizers_gauge={len(optimizers)}")
        print("  P:", " ".join(f"b={b:g}:{pressure[b]:.9f}" for b in betas))

    print("\nC1 exact 2/3-power defects (positive means exact subadditivity fails):")
    for m in range(2, 6):
        for r in range(2, 6):
            if m + r in results:
                defect = (results[m + r][0] ** (2 / 3)
                          - results[m][0] ** (2 / 3)
                          - results[r][0] ** (2 / 3))
                print(f"  {m}+{r}: {defect:+.9f}")

    print("\nC2 best principal normalized gap among hard-optimal parents:")
    for n in range(3, 8):
        parent = results[n][0] / n ** 1.5
        opts = results[n][1]
        for m in range(2, n):
            best = math.inf
            worst_over_parents = -math.inf
            for a in opts:
                child = min(induced_norm(a, n, S)
                            for S in combinations(range(n), m)) / m ** 1.5
                best = min(best, child)
                worst_over_parents = max(worst_over_parents, child)
            print(f"  n={n} m={m}: parent={parent:.9f} "
                  f"best_parent_child={best:.9f} worst_parent_child={worst_over_parents:.9f}")

    print("\nC3 same-beta additive pressure defects N*P_N-m*P_m-r*P_r:")
    for m in range(2, 6):
        for r in range(m, 6):
            if m + r in results:
                fields = []
                for beta in betas:
                    defect = ((m + r) * results[m + r][2][beta]
                              - m * results[m][2][beta]
                              - r * results[r][2][beta])
                    fields.append(f"b={beta:g}:{defect:+.6f}")
                print(f"  {m}+{r}: " + " ".join(fields))


if __name__ == "__main__":
    main()
