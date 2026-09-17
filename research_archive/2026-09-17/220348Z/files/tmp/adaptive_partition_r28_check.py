#!/usr/bin/env python3
"""Exact checks for the Wave 28 adaptive-partition collision memo."""

from __future__ import annotations

import itertools

import numpy as np


A8 = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, -1, 1, 1, -1, -1],
        [1, 1, 0, 1, -1, 1, -1, -1],
        [1, -1, 1, 0, -1, -1, -1, 1],
        [1, 1, -1, -1, 0, -1, 1, -1],
        [1, 1, 1, -1, -1, 0, 1, 1],
        [1, -1, -1, -1, 1, 1, 0, 1],
        [1, -1, -1, 1, -1, 1, 1, 0],
    ],
    dtype=np.int64,
)


def spin(word: str) -> np.ndarray:
    return np.array([1 if char == "+" else -1 for char in word], dtype=np.int64)


def qnorm(matrix: np.ndarray) -> int:
    n = len(matrix)
    return max(
        abs(int(x @ matrix @ x))
        for tail in itertools.product((-1, 1), repeat=n - 1)
        for x in (np.array((1,) + tail, dtype=np.int64),)
    )


def row_square(matrix: np.ndarray, x: np.ndarray) -> int:
    return int(x @ matrix @ matrix @ x)


def signature_blocks(family: list[np.ndarray]) -> list[tuple[int, ...]]:
    classes: dict[tuple[int, ...], list[int]] = {}
    base = family[0]
    for i in range(len(base)):
        signature = tuple(int(x[i] * base[i]) for x in family[1:])
        classes.setdefault(signature, []).append(i)
    return [tuple(block) for block in classes.values()]


def minimal_coset_rows(
    matrix: np.ndarray, family: list[np.ndarray]
) -> list[int]:
    base = family[0]
    blocks = signature_blocks(family)
    rows: list[int] = []
    for block_signs in itertools.product((-1, 1), repeat=len(blocks)):
        x = base.copy()
        for block_sign, block in zip(block_signs, blocks):
            x[list(block)] *= block_sign
        rows.append(row_square(matrix, x))
    return rows


def main() -> None:
    # Three individually cap-64 cuts.  Each is an exact completion of the
    # displayed order-four child ground, but their ternary product has row 72.
    family = [spin("+------+"), spin("+-----+-"), spin("+-+---+-")]
    selectors = [(0, 1, 2, 6), (0, 1, 2, 7), (0, 1, 2, 3)]

    assert [row_square(A8, x) for x in family] == [64, 64, 40]
    for x, selector in zip(family, selectors):
        idx = np.array(selector, dtype=np.int64)
        child = A8[np.ix_(idx, idx)]
        assert abs(int(x[idx] @ child @ x[idx])) == qnorm(child) == 8

    ternary_product = family[0] * family[1] * family[2]
    assert "".join("+" if value == 1 else "-" for value in ternary_product) == "+-+----+"
    assert row_square(A8, ternary_product) == 72

    blocks = signature_blocks(family)
    assert blocks == [(0, 1, 3, 4, 5), (2,), (6, 7)]
    rows = minimal_coset_rows(A8, family)
    assert sorted(set(rows)) == [40, 64, 72]
    assert max(rows) == 72

    print("PASS adaptive_partition_r28_check")
    print("child rows:", [row_square(A8, x) for x in family])
    print("signature blocks:", blocks)
    print("minimal containing coset row values:", sorted(set(rows)))
    print("ternary-product row:", row_square(A8, ternary_product))


if __name__ == "__main__":
    main()
