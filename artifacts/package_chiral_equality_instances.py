#!/usr/bin/env python3
"""Package the unrestricted and inherited-chiral equality masters."""

from __future__ import annotations

import numpy as np

from audit_dependent_4lift import B_MINUS
from bounded_chiral_equality_search import HERE
from search_chiral_equality_constraints import necessary_system


def main() -> None:
    (
        q,
        positives,
        negatives,
        edge_index,
        constraints,
        labels,
    ) = necessary_system(B_MINUS, add_one_fibre_cuts=True)
    matrix = constraints.A.tocsr()
    np.savez_compressed(
        HERE / "chiral_equality_unrestricted_master.npz",
        constraint_matrix_data=matrix.data,
        constraint_matrix_indices=matrix.indices,
        constraint_matrix_indptr=matrix.indptr,
        constraint_matrix_shape=np.asarray(matrix.shape),
        lower=constraints.lb,
        upper=constraints.ub,
        edges=np.asarray(edge_index.edges, dtype=np.int16),
        labels=np.asarray(labels),
        seed=B_MINUS,
        seed_q=q,
        positive_extremizers=np.asarray(positives, dtype=np.int8),
        negative_extremizers=np.asarray(negatives, dtype=np.int8),
        target=320,
        variable_convention=np.asarray(
            "z_uv in {0,1}; lifted edge sign b_uv=2*z_uv-1"
        ),
    )
    print(
        "saved unrestricted master:",
        matrix.shape,
        "nonzeros",
        matrix.nnz,
    )


if __name__ == "__main__":
    main()
