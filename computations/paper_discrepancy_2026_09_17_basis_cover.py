"""Exact finite checks for entangled tensor-Hadamard near-ground words."""

from __future__ import annotations

import json

import numpy as np


def main():
    h4 = np.ones((4, 4), dtype=np.int64) - 2*np.eye(4, dtype=np.int64)
    basis4 = np.array([[1, 1, 1, 1], [1, 1, -1, -1],
                      [1, -1, 1, -1], [1, -1, -1, 1]], dtype=np.int64)
    hp = np.ones((1, 1), dtype=np.int64)
    basis = np.ones((1, 1), dtype=np.int64)
    report = []
    for a in range(1, 5):
        hp = np.kron(h4, hp)
        basis = np.kron(basis4, basis)
        p = len(hp)
        eigenvalues = (basis @ hp * basis).sum(axis=1) // p
        eigenvalue = max(set(eigenvalues.tolist()), key=lambda value: np.sum(eigenvalues == value))
        selected = basis[eigenvalues == eigenvalue]
        assert len(selected) >= p//2
        n = 4*p
        u, v = selected[:2]
        x = np.concatenate((u, -u, v, -v))
        # Apply H4 tensor Hp without constructing its n-by-n matrix.
        action = (h4 @ x.reshape(4, p) @ hp.T).reshape(n)
        assert np.array_equal(action, -2*eigenvalue*x)
        diagonal = (-1)**(a+1)
        energy_twice = int(x @ action)-diagonal*n
        cap_twice = n*(2**(a+1)+1)
        deficit_twice = cap_twice-abs(energy_twice)
        assert deficit_twice in (0, 2*n)
        full_basis = np.kron(basis4, basis)
        max_absolute_correlation = int(np.abs(full_basis @ x).max())
        nearest_distance = (n-max_absolute_correlation)//2
        assert nearest_distance == n//4
        report.append({"n": n, "common_eigenspace_basis_size": len(selected),
                       "selected_pair_count": len(selected)**2,
                       "absolute_energy_deficit": deficit_twice//2,
                       "canonical_basis_nearest_distance": nearest_distance})
    print(json.dumps({"status": "all exact integer checks passed", "cases": report}, indent=2))


if __name__ == "__main__":
    main()
