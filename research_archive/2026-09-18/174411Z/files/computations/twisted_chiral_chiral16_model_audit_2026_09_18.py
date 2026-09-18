#!/usr/bin/env python3
"""Independent construction/constraint audits for the chiral16 CP-SAT model."""

import json
from pathlib import Path

import numpy as np

from twisted_chiral_chiral16_relaxation_2026_09_18 import build_constraints, normalize


ROOT = Path(__file__).resolve().parents[1]


def all_spins(n):
    ids = np.arange(1 << n, dtype=np.uint64)
    return 2*((ids[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int64)-1


def direct_histogram(a, b, d):
    c = b + np.diag(d)
    parent = np.block([[a, c], [c, -a]])
    x = all_spins(2*len(a))
    energy = np.einsum("bi,ij,bj->b", x, parent, x)//2
    hist = np.bincount(energy+100, minlength=201)
    return hist, int(abs(energy).max())


def main():
    rng = np.random.default_rng(202609181735)
    records = []
    for case in range(8):
        a, b = (rng.choice([-1, 1], size=(8, 8)) for _ in range(2))
        a = np.triu(a, 1); a += a.T
        b = np.triu(b, 1); b += b.T
        d = rng.choice([-1, 1], size=8)
        aa, bb = normalize(a, b, d)
        old, old_cap = direct_histogram(a, b, d)
        new, new_cap = direct_histogram(aa, bb, np.ones(8, dtype=np.int64))
        assert np.array_equal(old, new)
        records.append({"case": case, "cap": old_cap, "histogram_invariant": True,
                        "A": a.tolist(), "B": b.tolist(), "d": d.tolist()})
    hint = json.loads((ROOT / "computations/results/twisted_chiral_2026_09_18_profile_width_m8_class0_checkpoint.json").read_text())
    a, b = normalize(np.asarray(hint["child_matrix"]), np.asarray(hint["bridge_hollow_matrix"]), np.asarray(hint["d"]))
    constraints = []
    for target in [28, 30, 48]:
        variables, rows, bounds = build_constraints(8, target)
        for case in range(8):
            if case == 0:
                aa, bb = a, b
            else:
                aa, bb = (rng.choice([-1, 1], size=(8, 8)) for _ in range(2))
                aa = np.triu(aa, 1); aa += aa.T; aa[0, 1:] = aa[1:, 0] = 1
                bb = np.triu(bb, 1); bb += bb.T
            bits = np.array([((aa if kind == "A" else bb)[i, j]+1)//2 for kind, i, j in variables])
            lhs = rows.astype(np.int64) @ bits
            feasible = bool(np.all(lhs >= bounds[:, 0]) and np.all(lhs <= bounds[:, 1]))
            _, cap = direct_histogram(aa, bb, np.ones(8, dtype=np.int64))
            assert feasible == (cap <= target)
            constraints.append({"target": target, "case": case, "cap": cap, "constraints_feasible": feasible})
    result = {"seed": 202609181735, "gauge_tests": records, "constraint_tests": constraints,
              "method": "all 65536 parent spins, independent direct quadratic form; histogram equality under gauges"}
    destination = ROOT / "computations/results/twisted_chiral_chiral16_model_audit_2026_09_18.json"
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"gauge_tests_passed": len(records), "constraint_tests_passed": len(constraints)}))


if __name__ == "__main__":
    main()
