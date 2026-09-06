"""Exact counterexample to a proposed universal-quotient rigidity implication.

A Boolean spectral vector of the all-one seed specialization need not have
a Boolean orbit under all seed sign vectors, even when the conjugation
preserves full signs for every constant-diagonal seed.
"""

import argparse
import itertools
import json

import numpy as np

from continued_convergence_symmetric_switching_2026_09_06 import walsh
from continued_convergence_universal_overlap_2026_09_06 import all_seed_parents


def canonical_walsh_audit():
    records = []
    parity = lambda value: bin(value).count("1") % 2
    for d in (2, 4, 6, 8):
        n = 1 << d
        h = walsh(n)
        weights = np.array([bin(x).count("1") for x in range(n)], dtype=np.int64)
        q = (weights * (weights - 1) // 2) % 2
        linear_parity = (d // 2 + 1) % 2
        linears = [ell for ell in range(n) if parity(ell) == linear_parity]
        functions = np.array([
            1 - 2 * ((q + np.array([parity(ell & x) for x in range(n)])) % 2)
            for ell in linears
        ], dtype=np.int64)
        transformed = functions @ h
        root = 1 << (d // 2)
        orientations = transformed[:, 0] // root
        assert np.all(np.abs(orientations) == 1)
        assert np.array_equal(transformed, root * orientations[:, None] * functions)
        # Checking direction pairs suffices for all affine planes because
        # a quadratic second derivative is independent of the affine offset.
        admissible_pairs = 0
        for a, b in itertools.combinations(range(1, n), 2):
            beta = parity(a & b) ^ (parity(a) & parity(b))
            assert beta == int(q[a] ^ q[b] ^ q[a ^ b])
            if beta == 0:
                admissible_pairs += 1
                assert np.all(functions[:, 0] * functions[:, a] * functions[:, b] * functions[:, a ^ b] == 1)
        records.append({"binary_dimension": d, "order": n,
                        "common_projective_spectral_vectors": len(linears),
                        "admissible_direction_spaces": admissible_pairs // 3})
    return records


def verify(seed_order):
    s = 16
    h = walsh(s)
    quad = [0, 1, 6, 7]
    local_vector = np.array([1, -1, -1, 1], dtype=np.int64)
    outer_vector = np.array([1, -1, -1, -1, -1, 1, 1, 1,
                             -1, 1, 1, 1, -1, 1, 1, 1], dtype=np.int64)
    assert np.array_equal(h @ outer_vector, 4 * outer_vector)
    two_o = 2 * np.eye(s, dtype=np.int64)
    two_o[np.ix_(quad, quad)] -= np.outer(local_vector, local_vector)
    assert np.array_equal(two_o @ two_o.T, 4 * np.eye(s, dtype=np.int64))
    assert np.all(np.isin(two_o @ h, (-2, 2)))
    assert np.all(np.isin(two_o @ h @ two_o.T, (-4, 4)))

    # R = P_plus ⊗ I + P_minus ⊗ O on the first two seed fibres.
    # It is also the rank-one reflector I - ww^T/4 on eight coordinates.
    w = np.zeros(seed_order * s, dtype=np.int64)
    w[quad] = local_vector
    w[np.array(quad) + s] = -local_vector
    four_r = 4 * np.eye(seed_order * s, dtype=np.int64) - np.outer(w, w)
    assert np.array_equal(four_r @ four_r.T, 16 * np.eye(seed_order * s, dtype=np.int64))
    parents = all_seed_parents(seed_order, s)
    sixteen_new = four_r @ parents @ four_r.T
    assert np.all(np.isin(sixteen_new, (-16, 16)))
    assert np.array_equal(sixteen_new[0], 16 * parents[0])

    y = np.kron(np.ones(seed_order, dtype=np.int64), outer_vector)
    assert np.array_equal(four_r @ y, 4 * y)
    assert np.array_equal(parents[0] @ y, 4 * seed_order * y)
    seed_spin = np.ones(seed_order, dtype=np.int64)
    seed_spin[1] = -1
    bad_image_num = four_r @ np.kron(seed_spin, outer_vector)
    assert not np.all(np.abs(bad_image_num) == 4)

    return {
        "seed_order": seed_order,
        "outer_order": s,
        "all_seed_count": len(parents),
        "outer_quad": quad,
        "outer_reflector_vector": local_vector.tolist(),
        "boolean_outer_plus_eigenvector": outer_vector.tolist(),
        "twice_outer_image": (two_o @ outer_vector).tolist(),
        "failing_seed_spin": seed_spin.tolist(),
        "four_times_failing_image": bad_image_num.tolist(),
        "eight_coordinate_reflector_vector": w.tolist(),
        "all_one_specialization_unchanged": True,
        "canonical_walsh_audit": canonical_walsh_audit(),
        "scope": "Refutes J-Boolean-vector implies common-Boolean-orbit. Other common quotients may remain; no cap or convergence claim.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-order", type=int, default=4, choices=[2, 3, 4])
    parser.add_argument("--output")
    args = parser.parse_args()
    result = verify(args.seed_order)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
