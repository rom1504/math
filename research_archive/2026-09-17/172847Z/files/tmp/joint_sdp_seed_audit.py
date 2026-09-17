#!/usr/bin/env python3
"""Numerically audit degree-two SOS widths of saved seed signings."""

from __future__ import annotations


import json
from pathlib import Path

import cvxpy as cp
import numpy as np


ROOT = Path("/home/math/quadra")


def load(path: str, key: str) -> np.ndarray:
    return np.asarray(json.loads((ROOT / path).read_text())[key], dtype=float)


def endpoint(matrix: np.ndarray, maximize: bool) -> tuple[float, str]:
    n = len(matrix)
    gram = cp.Variable((n, n), symmetric=True)
    constraints = [gram >> 0, cp.diag(gram) == 1]
    expression = cp.sum(cp.multiply(matrix, gram)) / 2
    problem = cp.Problem(cp.Maximize(expression) if maximize else cp.Minimize(expression), constraints)
    value = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-9,
        tol_feas=1e-9,
        tol_gap_rel=1e-9,
        max_iter=1000,
    )
    return float(value), str(problem.status)


def exact_cap(matrix: np.ndarray) -> int:
    n = len(matrix)
    best = 0
    integer = matrix.astype(np.int64)
    for code in range(1 << (n - 1)):
        spin = np.ones(n, dtype=np.int64)
        for i in range(1, n):
            if code & (1 << (i - 1)):
                spin[i] = -1
        best = max(best, abs(int(spin @ integer @ spin) // 2))
    return best


def main() -> None:
    m8 = json.loads((ROOT / "computations/results/m8_minimizer_orbits.json").read_text())
    cases = [
        ("M6_conference", load("computations/results/exact_m6.json", "matrix")),
        ("M8_class0", np.asarray(m8["classes"][0]["representative_matrix"], dtype=float)),
        ("M8_class1", np.asarray(m8["classes"][1]["representative_matrix"], dtype=float)),
        ("M10_nonconference", load("computations/results/exact_m10.json", "matrix")),
        ("GF9_conference_n10", load("computations/results/conference_order10_gf9.json", "conference_matrix")),
        ("M12_witness", load("computations/results/heuristic_m12_seed20260731.json", "matrix")),
        ("M14_conference", load("computations/results/conference_completion_m13.json", "conference_matrix")),
    ]
    rows=[]
    for name,matrix in cases:
        upper,su=endpoint(matrix,True)
        lower,sl=endpoint(matrix,False)
        n=len(matrix)
        cap=exact_cap(matrix)
        rows.append({
            "name":name,"n":n,"cap":cap,"U1":upper,"L1":lower,
            "W1":(upper-lower)/2,"mu1":(upper+lower)/2,
            "R1":max(upper,-lower),"width_gap":(upper-lower)/2-cap,
            "universal_floor":n*np.sqrt(n-1)/2,"statuses":[su,sl],
            "conference_error":float(np.max(np.abs(matrix@matrix-(n-1)*np.eye(n)))),
        })
        print(json.dumps(rows[-1],sort_keys=True),flush=True)
    print(json.dumps({"schema":"joint-sdp-seed-audit-v1","classification":"numerical CLARABEL SDP; caps exact by exhaustive enumeration","rows":rows},indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
