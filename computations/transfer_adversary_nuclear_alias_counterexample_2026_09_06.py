#!/usr/bin/env python3
"""Finite cube checks of the scoped Hadamard alias counterexample.

Enumerating eta suffices: a product of any fixed number of independent
Boolean colors is again an independent Boolean eta vector. This does not
test or claim near-optimality of a hollow signing.
"""

import json
import math
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.linalg import hadamard


def main():
    records = []
    for n in (4, 8, 16):
        ids = np.arange(1 << n, dtype=np.uint64)
        eta = 1.0 - 2.0 * ((ids[:, None] >> np.arange(n, dtype=np.uint64)) & 1)
        h = hadamard(n).astype(float) / math.sqrt(n)
        x = eta @ h.T
        u = (x**3 - (3 - 2/n)*x) / math.sqrt(6)
        y = u @ h.T
        v = (n - 1)*(n - 2)/n**2
        cov_u = u.T @ u / len(eta)
        alias = y.T @ u / len(eta)
        assert np.max(np.abs(cov_u-v*np.eye(n))) < 2e-12
        assert np.max(np.abs(alias-v*h)) < 2e-12
        assert np.max(np.abs(x.T @ u / len(eta))) < 2e-12
        response = np.tanh(x+y)
        dual_cross = float(np.mean(response*y))
        # Any proposed diagonal coefficient from tanh'''/sqrt(6) is
        # bounded by 2/sqrt(6), whether or not it is evaluated here.
        diagonal_bound = 2 / math.sqrt(6*n)
        records.append({
            "n": n,
            "eta_patterns": len(eta),
            "exact_variance": v,
            "alias_row_l2": v,
            "actual_H_dual_cross": dual_cross,
            "any_bounded_diagonal_dual_cost": diagonal_bound,
            "floating_lower_bound_expression": dual_cross-diagonal_bound,
        })
    print(json.dumps({"status": "PASS", "cases": records}, indent=2))


if __name__ == "__main__":
    main()
