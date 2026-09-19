#!/usr/bin/env python3
"""Integer fourth-Gram diagonal optimization and explicit tensor witnesses."""
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def optimum_diagonal(D):
    q = len(D)
    square = D @ D
    c = np.diag(square @ D)
    off = square - np.diag(np.diag(square))
    omega0 = int(np.sum(off * off))
    order = np.argsort(c, kind="stable")
    prefix = 0
    best = None
    for k in range(q + 1):
        if k:
            prefix += int(c[order[k - 1]])
        s = 2 * k - q
        omega = omega0 + 2 * q * (q - 2) + 2 * s * s + 4 * (2 * prefix - int(c.sum()))
        h = -np.ones(q, dtype=np.int64)
        h[order[:k]] = 1
        H = D + np.diag(h)
        gram_error = H @ H - q * np.eye(q, dtype=np.int64)
        assert omega == int(np.sum(gram_error * gram_error))
        if best is None or omega < best[0]:
            best = (omega, h.copy())
    return best, c, omega0


def main():
    rng = np.random.default_rng(2026091918)
    tested = 0
    tensor_cases = []
    for q in range(2, 9):
        for trial in range(12):
            X = np.triu(rng.choice([-1, 1], (q, q)), 1)
            D = X + X.T
            (best, hbest), c, omega0 = optimum_diagonal(D)
            values = []
            for htuple in itertools.product([-1, 1], repeat=q):
                h = np.array(htuple, dtype=np.int64)
                H = D + np.diag(h)
                G = H @ H
                error = G - q * np.eye(q, dtype=np.int64)
                direct = int(np.sum(error * error))
                formula = omega0 + 2*q*(q-2) + 2*int(h.sum())**2 + 4*int(c @ h)
                assert formula == direct
                values.append(direct)
                tested += 1
            assert best == min(values)
            if trial == 0:
                H = D + np.diag(hbest)
                K = np.kron(H, H)
                z = H.reshape(-1)
                trace4 = int(np.trace(H @ H @ H @ H))
                assert int(z @ K @ z) == trace4 == q**3 + best
                hollow = K - np.diag(np.diag(K))
                assert int(z @ hollow @ z) == trace4 - int(np.trace(H))**2
                tensor_cases.append({"q": q, "H": H.tolist(), "trace4": trace4,
                                     "trace": int(np.trace(H)), "integer_vector_identity_pass": True})

    parents = []
    seen = set()
    source = ROOT / "computations/results/twisted_chiral_2026_09_19_joint_target38.json"
    for row in json.loads(source.read_text())["records"]:
        if row["matrix_sha256"] in seen:
            continue
        seen.add(row["matrix_sha256"])
        D = np.array(row["parent_matrix"], dtype=np.int64)
        q = len(D)
        (omega, h), c, omega0 = optimum_diagonal(D)
        n = q // 2
        assert np.array_equal(c[:n], -c[n:])
        chiral_formula = omega0 + 2*q*(q-2) - 4*int(np.abs(c).sum())
        assert omega == chiral_formula
        H = D + np.diag(h)
        gram = H @ H
        K = np.kron(H, H)
        z = H.reshape(-1)
        assert int(z @ K @ z) == q**3 + omega
        parents.append({"source": str(source.relative_to(ROOT)), "matrix_sha256": row["matrix_sha256"],
                        "parent_cap": row["cap"], "order": q,
                        "triangle_row_sums_c": c.tolist(), "omega0": omega0,
                        "minimum_diagonal_gram_defect": omega, "minimizing_diagonal": h.tolist(),
                        "rho_lower_numerator": q**3 + omega, "rho_lower_denominator": q**3,
                        "tensor_square_hollow_witness_energy_numerator": q**3 + omega - int(h.sum())**2,
                        "tensor_square_hollow_witness_energy_denominator": 2,
                        "all_diagonals_covered_by_sorted_prefix_proof": True,
                        "chiral_linear_formula_pass": True,
                        "scope": "ordinary full-sign diagonal completions and literal tensor powers only"})
    out = {"status": "integer certificates; not a convergence theorem",
           "seed": 2026091918, "diagonal_formula_checks": tested,
           "random_hollow_inputs": 84, "tensor_identity_cases": tensor_cases,
           "order20_parents": parents}
    path = ROOT / "computations/results/twisted_chiral_tensor_gram_2026_09_19.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"diagonal_formula_checks": tested,
                      "parents": [{k: r[k] for k in ["matrix_sha256", "minimum_diagonal_gram_defect", "rho_lower_numerator", "rho_lower_denominator"]} for r in parents]}))


if __name__ == "__main__":
    main()
