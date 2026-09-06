"""Finite integer replay of the root's same-order switched-clique repair."""
import importlib.util
import json
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPINS = {}


def spins(n):
    if n not in SPINS:
        x = np.ones((1 << (n-1), n), dtype=np.int16)
        masks = np.arange(len(x), dtype=np.uint32)
        for k in range(1, n):
            x[:, k] = 1 - 2*((masks >> (k-1)) & 1).astype(np.int16)
        SPINS[n] = x
    return SPINS[n]


def extrema(a):
    if not len(a):
        return 0, 0, np.array([], dtype=np.int16), np.array([], dtype=np.int16)
    x = spins(len(a))
    twice = (x @ a * x).sum(axis=1, dtype=np.int64)
    assert np.all(twice % 2 == 0)
    energy = twice // 2
    return int(energy.max()), -int(energy.min()), x[energy.argmax()], x[energy.argmin()]


def verify(name, matrix):
    a = np.asarray(matrix, dtype=np.int16)
    n = len(a)
    p, r, x, y = extrema(a)
    if p < r:
        a = -a
        p, r, x, y = extrema(a)
    delta = p-r
    k = 0
    while comb(k, 2) < delta:
        k += 1
    local = x * (a @ x)
    assert np.all(local >= 0) and int(local.sum()) == 2*p
    low = np.flatnonzero(n*local <= 8*p)
    # KG>3/2, so the following rational Frobenius test is a sufficient
    # certificate that U=all vertices meets the root's spectral-core bound.
    full_core_valid = n**3*(n-1) <= 144*p*p
    size_condition = k <= n//4
    theorem_admissible = full_core_valid and size_condition
    if theorem_admissible:
        assert len(low) >= k
        selected = low[:k]
    else:
        selected = np.argsort(local, kind="stable")[:k]
    sub = a[np.ix_(selected, selected)]
    sp, sr, _, _ = extrema(sub)
    u = max(sp, sr)
    rowmass = int(local[selected].sum())
    changed = a.copy()
    if k:
        switched = y[selected]
        changed[np.ix_(selected, selected)] = -switched[:, None]*switched[None, :] + np.eye(k, dtype=np.int16)
    pp, rr, _, _ = extrema(changed)
    c, b = comb(k, 2), k//2
    assert p-rowmass-u+b <= pp <= p+u+b
    assert r+c-u <= rr <= r+c+u
    assert max(pp, rr) >= p-u
    assert max(pp, rr) <= p+u+k
    assert abs(pp-rr) <= rowmass+2*u+k
    if theorem_admissible:
        assert n*rowmass <= 8*k*p
        assert n*u <= 6*k*p  # stronger than 4 KG k Q/n here
    if delta:
        assert c-delta < k-1 and k*k <= 16*delta
    return {
        "name": name, "order": n, "matrix_after_global_orientation": a.tolist(),
        "P": p, "R": r, "delta": delta, "k": k,
        "selected": selected.tolist(), "old_block_cap": u,
        "positive_ground_rowmass_on_selected": rowmass,
        "new_P": pp, "new_R": rr,
        "full_core_frobenius_certificate": full_core_valid,
        "k_le_floor_n_over_4": size_condition,
        "root_theorem_admissible": theorem_admissible,
        "status": "All exact two-sided finite bounds pass",
    }


def main():
    records = []
    data = json.loads((ROOT / "computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json").read_text())
    for case in data["cases"]:
        records.append(verify("external_upper_"+str(case["order"]), case["matrix"]))
    spec = importlib.util.spec_from_file_location("stored_minimizers", ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    for name, matrix, _ in old.cases():
        records.append(verify(name, matrix))
    rng = np.random.default_rng(260906)
    for n in range(12, 17):
        for trial in range(8):
            a = np.zeros((n, n), dtype=np.int16)
            for i in range(n):
                for j in range(i+1, n):
                    a[i, j] = a[j, i] = 2*int(rng.integers(2))-1
            records.append(verify(f"random_n{n}_trial{trial}", a))
    result = {
        "scope": "All finite block inequalities tested; spectral-core theorem asserted only on marked admissible records",
        "seed": 260906, "cases": records,
        "finite_bound_tests": len(records),
        "root_theorem_admissible_tests": sum(case["root_theorem_admissible"] for case in records),
        "nontrivial_admissible_repairs": sum(case["root_theorem_admissible"] and case["delta"]>0 for case in records),
    }
    path = ROOT / "computations/results/transfer_adversary_same_order_orientation_repair_2026_09_06.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key:value for key,value in result.items() if key != "cases"}))
    print("PASS: exact same-order repair extrema and selected-core bounds")


if __name__ == "__main__":
    main()
