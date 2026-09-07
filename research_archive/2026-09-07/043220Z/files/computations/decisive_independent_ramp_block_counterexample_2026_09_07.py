"""Exact finite checks for the deterministic ramp rank/coherence obstruction."""

from fractions import Fraction as Q
from itertools import product
import json
import numpy as np


def main():
    f = np.kron(np.array([[1, 1], [1, -1]], dtype=np.int64),
                np.array([[1, 1], [1, -1]], dtype=np.int64))
    signs = np.array([-1, -1, 1, -1], dtype=np.int64)
    balanced = signs[:, None]*f*signs[None, :]
    assert np.array_equal(balanced, balanced.T)
    assert np.array_equal(balanced@balanced.T, 4*np.eye(4, dtype=np.int64))
    assert balanced.sum() == 0
    hb = f
    hl = np.kron(balanced, np.kron(f, f))
    b, ell = len(hb), len(hl)
    n = b*ell
    root_n, root_ell = 16, 8
    assert root_n**2 == n and root_ell**2 == ell and ell == 4*root_n
    assert np.array_equal(hl@hl.T, ell*np.eye(ell, dtype=np.int64))
    assert hl.sum() == 0
    core = np.kron(hb-np.diag(np.diag(hb)), hl)
    a = core+np.kron(np.eye(b, dtype=np.int64),
                    np.ones((ell, ell), dtype=np.int64)-np.eye(ell, dtype=np.int64))
    assert np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0)
    assert np.all(np.abs(a+np.eye(n, dtype=np.int64)) == 1)
    e0 = n*(ell-1)//2
    block_words = 0
    for z in product((-1, 1), repeat=b):
        x = np.repeat(z, ell)
        assert x@core@x == 0
        assert x@a@x == 2*e0
        block_words += 1

    # Exact sample comparison with a non-block coordinate projector.
    r = 3
    coordinate_diagonal = np.array([1]*r+[0]*(n-r), dtype=np.int64)
    average_projection = sum(
        np.repeat(z, ell)@(coordinate_diagonal*np.repeat(z, ell))
        for z in product((-1, 1), repeat=b)
    )//block_words
    assert average_projection == r <= ell*r

    # The block projector, represented by ell*Pi, is sign-compatible.
    projector_numerator = np.kron(np.eye(b, dtype=np.int64),
                                 np.ones((ell, ell), dtype=np.int64))
    assert np.all(a*projector_numerator >= 0)
    assert np.array_equal(projector_numerator@projector_numerator,
                          ell*projector_numerator)

    full_repair = np.kron(hb, hl)
    assert np.array_equal(full_repair@full_repair.T, n*np.eye(n, dtype=np.int64))
    repair = full_repair-np.diag(np.diag(full_repair))
    edits = int(np.count_nonzero(np.triu(a != repair, 1)))
    assert edits <= e0
    assert np.all((a != repair)[core != 0] == 0)
    c = Q(ell-2-2*root_n-2*root_ell, ell-1+root_n+root_ell)
    assert c == Q(14, 87)
    print(json.dumps({"status": "exact checks passed", "n": n,
                      "blocks": b, "block_size": ell,
                      "all_block_constant_words": block_words,
                      "block_constant_energy": e0,
                      "finite_zero_error_capture_fraction": str(c),
                      "explicit_sparse_repair_edits": edits,
                      "block_projector_sign_compatible": True}, indent=2))


if __name__ == "__main__":
    main()
