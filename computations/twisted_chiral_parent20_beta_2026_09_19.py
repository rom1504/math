#!/usr/bin/env python3
"""Exact direct response norms for the retained inequivalent order20 parents."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def inspect(record):
    matrix = np.asarray(record["parent_matrix"], dtype=np.int16)
    n = len(matrix)
    assert np.array_equal(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0)
    assert np.all(np.abs(matrix + np.eye(n, dtype=np.int16)) == 1)
    beta, positive, negative = -1, -1, -1
    beta_count = positive_count = negative_count = 0
    product_hist = np.zeros(n + 1, dtype=np.int64)
    field_zero_hist = np.zeros(n + 1, dtype=np.int64)
    witness = None
    for first in range(0, 1 << (n - 1), 4096):
        words = np.arange(first, min(first + 4096, 1 << (n - 1)), dtype=np.int64)
        spins = np.ones((len(words), n), dtype=np.int16)
        spins[:, 1:] = 1 - 2*((words[:, None] >> np.arange(n - 1)) & 1)
        fields = spins @ matrix
        values = np.abs(fields).sum(axis=1)
        energy = (spins * fields).sum(axis=1) // 2
        top, bottom = int(energy.max()), int(-energy.min())
        if top > positive:
            positive, positive_count = top, 0
        if top == positive:
            positive_count += int(np.count_nonzero(energy == positive))
        if bottom > negative:
            negative, negative_count = bottom, 0
        if bottom == negative:
            negative_count += int(np.count_nonzero(energy == -negative))
        local = int(values.max())
        if local > beta:
            beta, beta_count = local, 0
            product_hist[:] = 0
            field_zero_hist[:] = 0
            idx = int(values.argmax())
            x = spins[idx].copy()
            y = np.where(fields[idx] >= 0, 1, -1).astype(np.int16)
            assert int(x @ matrix @ y) == beta
            witness = [x.tolist(), y.tolist()]
        if local == beta:
            selected = values == beta
            beta_count += int(np.count_nonzero(selected))
            # Order20 row sums are odd, so a maximizing response has no zero fields.
            assert np.all(fields[selected] != 0)
            products = spins[selected] * np.sign(fields[selected])
            sizes = (products == 1).sum(axis=1)
            product_hist += np.bincount(sizes, minlength=n + 1)
            field_zero_hist += np.bincount((fields[selected] == 0).sum(axis=1), minlength=n + 1)
    assert positive == negative == int(record["cap"])
    x, y = witness
    scalar = sum(int(x[i]) * int(matrix[i, j]) * int(y[j])
                 for i in range(n) for j in range(n))
    assert scalar == beta
    return {
        "matrix_sha256_from_source": record["matrix_sha256"],
        "matrix": matrix.tolist(), "order": n,
        "seed_cap": max(positive, negative), "positive_maximum": positive,
        "negative_maximum": negative,
        "positive_extremizer_projective_count": positive_count,
        "negative_extremizer_projective_count": negative_count,
        "beta": beta, "beta_to_quadratic_cap_ratio": beta/max(positive, negative),
        "beta_maximizing_response_projective_count": beta_count,
        "beta_maximizer_product_positive_coordinate_histogram": product_hist.tolist(),
        "beta_maximizer_zero_field_histogram": field_zero_hist.tolist(),
        "bilinear_witness": witness, "scalar_witness_value": scalar,
        "projective_spins_enumerated": 1 << (n - 1),
    }


def main():
    data = json.loads((ROOT / "computations/results/twisted_chiral_2026_09_19_joint_target38.json").read_text())
    unique = {}
    for item in data["records"]:
        unique.setdefault(item["matrix_sha256"], item)
    results = [inspect(item) for item in unique.values()]
    for item in results:
        print(item["matrix_sha256_from_source"], item["seed_cap"], item["beta"],
              item["beta_maximizing_response_projective_count"], flush=True)
    output = {
        "normalization": "H=x^T A x/2; beta=max_x ||Ax||_1",
        "method": "all projective spins, direct integer matrix multiplication; independent scalar witness check",
        "records": results,
    }
    (ROOT / "computations/results/twisted_chiral_parent20_beta_2026_09_19.json").write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
