#!/usr/bin/env python3
"""Verify the rejected candidate and packaged lazy equality instance."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from audit_dependent_4lift import B_MINUS
from bounded_chiral_equality_search import (
    HERE,
    inherited_chiral_rows,
)
from search_chiral_equality_constraints import necessary_system


def main() -> None:
    rejected = np.load(HERE / "chiral_equality_rejected_candidate.npz")
    b = rejected["matrix"].astype(np.int64)
    spin = rejected["violating_spin"].astype(np.int64)
    energy = int(rejected["violating_energy"])
    assert b.shape == (48, 48)
    assert np.array_equal(b, b.T)
    assert np.all(np.diag(b) == 0)
    assert set(b[~np.eye(48, dtype=bool)]) == {-1, 1}
    assert int(spin @ b @ spin) == energy == 440

    (
        _,
        _,
        _,
        edge_index,
        base_constraints,
        _,
    ) = necessary_system(B_MINUS, add_one_fibre_cuts=False)
    edge_values = np.asarray(
        [(b[u, v] + 1) // 2 for u, v in edge_index.edges],
        dtype=np.float64,
    )
    lhs = base_constraints.A @ edge_values
    assert np.all(lhs >= base_constraints.lb - 1e-9)
    assert np.all(lhs <= base_constraints.ub + 1e-9)

    sym_a, sym_l, sym_u, inherited = inherited_chiral_rows(
        edge_index, np.eye(4, dtype=np.int8)
    )
    sym_lhs = sym_a @ edge_values
    assert np.all(sym_lhs >= sym_l - 1e-9)
    assert np.all(sym_lhs <= sym_u + 1e-9)
    assert np.array_equal(inherited.T @ b @ inherited, -b)

    instance = np.load(HERE / "chiral_equality_lazy_instance.npz")
    shape = tuple(map(int, instance["constraint_matrix_shape"]))
    matrix = csr_matrix(
        (
            instance["constraint_matrix_data"],
            instance["constraint_matrix_indices"],
            instance["constraint_matrix_indptr"],
        ),
        shape=shape,
    )
    assert shape == (1508, 1128)
    assert len(instance["lower"]) == shape[0]
    assert len(instance["upper"]) == shape[0]
    assert instance["lazy_spins"].shape == (1, 48)
    assert np.array_equal(instance["lazy_spins"][0], spin)
    assert matrix.shape[1] == len(edge_index.edges)

    print("rejected inherited-chiral candidate: verified")
    print("candidate satisfies 1507 master constraints")
    print("exact separating energy:", energy)
    print("packaged lazy instance:", shape, "with one violated-spin cut")
    print("status:", instance["reason"].item())


if __name__ == "__main__":
    main()
