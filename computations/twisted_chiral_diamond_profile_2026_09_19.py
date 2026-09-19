#!/usr/bin/env python3
"""Exact diamond norm from 3^n partial spins and 2^n subset profiles."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def analyze(matrix, label):
    matrix = np.array(matrix, dtype=np.int16)
    n = len(matrix)
    states = (((np.arange(3**n)[:, None] // (3**np.arange(n))) % 3) - 1).astype(np.int16)
    supports = ((states != 0) * (1 << np.arange(n))).sum(axis=1)
    energy = np.einsum("bi,ij,bj->b", states, matrix, states) // 2
    upper = np.zeros(1 << n, dtype=np.int16)
    lower = upper.copy()
    np.maximum.at(upper, supports, energy)
    np.maximum.at(lower, supports, -energy)
    masks = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(np.int16)
    best, witness = -1, None
    full = (1 << n) - 1
    for start in range(0, len(states), 256):
        fields = np.abs(states[start:start + 256] @ matrix)
        score = fields @ (1 - masks).T + lower[None, :] + upper[full ^ supports[start:start + 256], None]
        value = int(score.max())
        if value > best:
            i, j = np.unravel_index(score.argmax(), score.shape)
            best, witness = value, (start + int(i), int(j))
    b_index, c_support = witness
    b = states[b_index]
    a_support = full ^ supports[b_index]
    a = states[np.flatnonzero((supports == a_support) & (energy == upper[a_support]))[0]]
    c = states[np.flatnonzero((supports == c_support) & (energy == -lower[c_support]))[0]]
    fields = matrix @ b
    d = np.where(masks[c_support] == 0, np.where(fields >= 0, 1, -1), 0)
    assert np.all(np.abs(a) + np.abs(b) == 1)
    assert np.all(np.abs(c) + np.abs(d) == 1)
    value = int((a @ matrix @ a - c @ matrix @ c) // 2 + b @ matrix @ d)
    assert value == best
    beta = int(np.abs(states[supports == full] @ matrix).sum(axis=1).max())
    width_twice = int(upper[full] + lower[full])
    assert 3 * best <= 4 * beta
    assert 2 * best <= 2 * width_twice + beta
    result = {
        "label": label, "matrix": matrix.tolist(), "order": n,
        "seed_cap": int(max(upper[full], lower[full])),
        "positive_maximum": int(upper[full]), "negative_maximum": int(lower[full]),
        "beta": beta, "diamond_cap": best,
        "diamond_witness": [z.tolist() for z in [a, b, c, d]],
        "partial_spin_count": 3**n, "subset_count": 1 << n,
    }
    print(label, result["seed_cap"], beta, best, flush=True)
    return result


def main():
    records = []
    for n in range(3, 11):
        orbit_file = ROOT / f"computations/results/m{n}_minimizer_orbits.json"
        if orbit_file.exists():
            data = json.loads(orbit_file.read_text())
            seeds = [(f"m{n}_class{item['class']}", item["representative_matrix"]) for item in data["classes"]]
        else:
            data = json.loads((ROOT / f"computations/results/exact_m{n}.json").read_text())
            seeds = [(f"m{n}_stored", data["matrix"])]
        records.extend(analyze(a, label) for label, a in seeds)
    conference = json.loads((ROOT / "computations/results/conference_order10_gf9.json").read_text())
    records.append(analyze(conference["conference_matrix"], "m10_conference"))
    output = {
        "normalization": "H_A=x^T A x/2; Delta=max_diamonds |H(a)-H(c)+b^T A d|",
        "exact_formula": "max_(b in {0,+-1}^n,J subset[n]) U(supp(b)^c)+V(J)+sum_(i notinJ)|(Ab)_i|",
        "all_integer_checks_pass": True,
        "records": records,
    }
    path = ROOT / "computations/results/twisted_chiral_diamond_profile_2026_09_19.json"
    path.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
