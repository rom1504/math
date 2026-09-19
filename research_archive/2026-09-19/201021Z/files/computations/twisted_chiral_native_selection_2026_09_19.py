#!/usr/bin/env python3
"""Exact small-seed audit of a bilinear-ground-pair switching rule."""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def inspect(matrix):
    matrix = np.asarray(matrix, dtype=np.int64)
    n = len(matrix)
    spins = np.array(list(itertools.product([-1, 1], repeat=n)), dtype=np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
    bilinear = spins @ matrix @ spins.T
    beta = int(np.max(np.abs(bilinear)))
    i, j = np.nonzero(np.abs(bilinear) == beta)
    candidates = {tuple((spins[x] * spins[y] * spins[x, 0] * spins[y, 0]).tolist())
                  for x, y in zip(i, j)}
    caps = {}
    for switch in spins[spins[:, 0] == 1]:
        b = matrix * switch[:, None] * switch[None, :]
        cap = int(np.max(np.abs(energies[:, None] - energies[None, :] + spins @ b @ spins.T)))
        caps[tuple(switch.tolist())] = cap
    selected = min(caps[s] for s in candidates)
    unrestricted = min(caps.values())
    return {
        "order": n, "matrix": matrix.tolist(), "beta": beta,
        "seed_cap": int(np.max(np.abs(energies))),
        "bilinear_ground_product_count": len(candidates),
        "ground_product_core_minimum": selected,
        "unrestricted_switch_core_minimum": unrestricted,
        "ground_product_best_switch": list(next(s for s in candidates if caps[s] == selected)),
        "unrestricted_best_switch": list(next(s for s in caps if caps[s] == unrestricted)),
    }


def main():
    records = []
    counts = {}
    for n in range(3, 7):
        edges = list(itertools.combinations(range(1, n), 2))
        worst = None
        failures = []
        for word in range(1 << len(edges)):
            a = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
            for bit, (i, j) in enumerate(edges):
                a[i, j] = a[j, i] = 1 if word & (1 << bit) else -1
            result = inspect(a)
            if worst is None or result["ground_product_core_minimum"] * worst["beta"] > worst["ground_product_core_minimum"] * result["beta"]:
                worst = result
            if 2 * result["ground_product_core_minimum"] > 3 * result["beta"]:
                failures.append(result)
        counts[n] = 1 << len(edges)
        records.append({"order": n, "worst": worst, "three_halves_failures": failures[:1],
                        "three_halves_failure_count": len(failures)})
        print(n, counts[n], worst["ground_product_core_minimum"], worst["beta"], len(failures), flush=True)
    bipolar = []
    for m in range(2, 9):
        j = np.ones((m, m), dtype=np.int64)
        a = np.block([[j - np.eye(m, dtype=np.int64), j],
                      [j, -j + np.eye(m, dtype=np.int64)]])
        spins = np.array(list(itertools.product([-1, 1], repeat=2*m)), dtype=np.int64)
        fields = spins @ a
        values = np.abs(fields).sum(axis=1)
        assert int(values.max()) == 2*m*m
        maximizing = np.flatnonzero(values == values.max())
        assert len(maximizing) == 4
        assert np.all(fields[maximizing] != 0)
        products = spins[maximizing] * np.sign(fields[maximizing])
        assert np.all(products == products[:, :1])
        x = np.ones(2*m, dtype=np.int64)
        y = np.concatenate([np.ones(m, dtype=np.int64), -np.ones(m, dtype=np.int64)])
        core_witness = int((x @ a @ x - y @ a @ y)//2 + x @ a @ y)
        assert core_witness == 4*m*m - 2*m
        # The exact spectral upper is strictly below the next integer.
        assert 8*m*m*(2*m*m - 2*m + 1) < (core_witness + 1)**2
        bipolar.append({"half_order": m, "beta": 2*m*m,
                        "bilinear_maximizing_second_spin_count": len(maximizing),
                        "all_ground_products_constant": True,
                        "untwisted_core_cap": core_witness,
                        "spectral_integer_upper_pass": True})
    path = ROOT / "computations/results/twisted_chiral_native_selection_2026_09_19.json"
    path.write_text(json.dumps({"normalization": "H=x^T A x/2; core means no bridge matching", "root_gauged_seed_counts": counts, "records": records,
                               "bipolar_ground_rule_asymptotic_counterfamily_checks": bipolar}, indent=2) + "\n")


if __name__ == "__main__":
    main()
