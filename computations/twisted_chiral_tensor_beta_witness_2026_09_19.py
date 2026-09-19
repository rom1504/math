#!/usr/bin/env python3
"""A lower witness, not an optimum claim, for beta(H12 tensor-cubed).

The finite conference-six beta deficit does not persist in this natural
Hadamard tensor family.  The retained witness is checked by dense integer
matrix multiplication independently of the tensor response used to find it.
"""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "computations/results/twisted_chiral_tensor_beta_witness_2026_09_19.json"
C = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, 1, 1, -1, -1],
    [1, 1, 0, -1, -1, 1],
    [1, 1, -1, 0, 1, -1],
    [1, -1, -1, 1, 0, 1],
    [1, -1, 1, -1, 1, 0],
], dtype=np.int16)
I = np.eye(6, dtype=np.int16)
H = np.block([[C - I, C + I], [C + I, -C + I]])


def tensor_response(words):
    words = np.einsum("ai,bijk->bajk", H, words, optimize=True)
    words = np.einsum("aj,bijk->biak", H, words, optimize=True)
    return np.einsum("ak,bijk->bija", H, words, optimize=True)


def small_product_witness(factor, seed, steps, target):
    matrix = np.kron(H, factor).astype(np.int64)
    n = len(matrix)
    rng = np.random.default_rng(seed)
    words = rng.choice(np.array([-1, 1], dtype=np.int16), size=(1024, n))
    best, best_word = -1, None
    for _ in range(steps):
        fields = words @ matrix
        values = np.abs(fields).sum(axis=1)
        row = int(np.argmax(values))
        if int(values[row]) > best:
            best, best_word = int(values[row]), words[row].copy()
        responses = np.where(fields >= 0, 1, -1).astype(np.int16)
        updated = np.where(responses @ matrix >= 0, 1, -1).astype(np.int16)
        for row in np.flatnonzero(np.all(updated == words, axis=1)):
            updated[row, rng.choice(n, size=2, replace=False)] *= -1
        words = updated
    x = np.where(matrix @ best_word >= 0, 1, -1)
    assert int(x @ matrix @ best_word) == best == target
    # Deliberately separate plain-integer reconstruction for these small seeds.
    scalar_score = sum(int(x[i]) * int(matrix[i, j]) * int(best_word[j])
                       for i in range(n) for j in range(n))
    assert scalar_score == target
    return {"order": n, "sylvester_factor_order": len(factor),
            "bilinear_lower_bound": best, "normalized_lower_bound": best / n**1.5,
            "x": x.tolist(), "y": best_word.tolist(),
            "scalar_integer_reconstruction_pass": True,
            "optimum_claimed": False}


def main():
    assert np.array_equal(C @ C, 5 * I)
    assert np.array_equal(H, H.T)
    assert np.array_equal(H @ H, 12 * np.eye(12, dtype=np.int16))
    assert np.trace(H) == 0
    signs12 = np.array(list(itertools.product([-1, 1], repeat=12)), dtype=np.int16)
    beta12 = int(np.abs(signs12 @ H).sum(axis=1).max())
    assert beta12 == 36

    rng = np.random.default_rng(2026091902)
    words = rng.choice(np.array([-1, 1], dtype=np.int16), size=(48, 12, 12, 12))
    best, best_word = -1, None
    for _ in range(100):
        fields = tensor_response(words)
        values = np.abs(fields).sum(axis=(1, 2, 3))
        row = int(np.argmax(values))
        if int(values[row]) > best:
            best, best_word = int(values[row]), words[row].copy()
        responses = np.where(fields >= 0, 1, -1).astype(np.int16)
        updated = np.where(tensor_response(responses) >= 0, 1, -1).astype(np.int16)
        stuck = np.all(updated == words, axis=(1, 2, 3))
        for row in np.flatnonzero(stuck):
            updated[row].reshape(-1)[rng.choice(1728, size=8, replace=False)] *= -1
        words = updated

    # Independent dense scalar-integer reconstruction of the reported witness.
    matrix = np.kron(np.kron(H, H), H).astype(np.int64)
    y = best_word.reshape(-1).astype(np.int64)
    dense_fields = matrix @ y
    assert np.array_equal(dense_fields, tensor_response(best_word[None])[0].reshape(-1))
    x = np.where(dense_fields >= 0, 1, -1)
    score = int(x @ matrix @ y)
    assert score == best == 66592
    assert score / 1728**1.5 > 0.927

    # H tensor H has a flat Boolean eigenvector vec(H), so its beta is 12^3.
    square = np.kron(H, H).astype(np.int64)
    flat = H.reshape(-1).astype(np.int64)
    assert np.array_equal(square @ flat, 12 * flat)
    h2 = np.array([[1, 1], [1, -1]], dtype=np.int16)
    h4 = np.kron(h2, h2)
    bent4 = np.array([1, 1, 1, -1])
    assert np.array_equal(h4 @ bent4, 2 * bent4)
    small_products = [
        small_product_witness(h2, 2026091904, 15, 112),
        small_product_witness(h4, 2026091903, 30, 324),
    ]

    data = {
        "scope": "Explicit bilinear lower witness only; no tensor-cube optimum claimed.",
        "conference6": C.tolist(),
        "hadamard12": H.tolist(),
        "beta_hadamard12_exact": beta12,
        "order": 1728,
        "tensor_power": 3,
        "bilinear_lower_bound": score,
        "normalized_lower_bound": score / 1728**1.5,
        "x": x.tolist(),
        "y": y.tolist(),
        "dense_integer_reconstruction_pass": True,
        "tensor_square_flat_eigenvector_pass": True,
        "extension": "Product with vec(H12) extends this normalized lower bound to every odd tensor power >=3. Even powers have normalized beta exactly1. Hollowing changes beta by at most the order.",
        "small_sylvester_products": small_products,
        "sylvester_extension": "Product with the flat H4 sign eigenvector extends the order24 bound to all odd Sylvester exponents, and the order48 bound to all positive even exponents.",
    }
    OUTPUT.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({k: v for k, v in data.items() if k not in {"conference6", "hadamard12", "x", "y", "small_sylvester_products"}}, indent=2))
    for item in small_products:
        print(json.dumps({k: v for k, v in item.items() if k not in {"x", "y"}}))


if __name__ == "__main__":
    main()
