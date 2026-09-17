"""Exact small block-dictionary response enumeration and spectral replay."""

from fractions import Fraction
from itertools import product
import json
import math
import numpy as np


def words(n):
    return list(product((-1, 1), repeat=n))


def block_code(k, p):
    # Fix c_0=1 to remove the sole duplicate representation.
    result = []
    for tail in words(k - 1):
        c = (1,) + tail
        for v in words(p):
            result.append(tuple(a * b for a in c for b in v))
    return result


def dephased_matrices(k, p):
    # Row/column sign switching preserves the uniform rank-one query law.
    for signs in words((k - 1) * (p - 1)):
        out = np.ones((k, p), dtype=np.int16)
        if signs:
            out[1:, 1:] = np.asarray(signs, dtype=np.int16).reshape(k - 1, p - 1)
        yield out.reshape(-1)


def exact_dictionary_case(shapes):
    block_codes = [block_code(k, p) for k, p in shapes]
    code = np.asarray([sum(parts, ()) for parts in product(*block_codes)], dtype=np.int16)
    matrix_blocks = [list(dephased_matrices(k, p)) for k, p in shapes]
    minimum = None
    minimum_count = 0
    total = 0
    rank_one_value = Fraction(int(np.abs(code.sum(axis=1)).sum()), len(code))
    for pieces in product(*matrix_blocks):
        physical = np.concatenate(pieces)
        numerator = int(np.abs(code @ physical).sum())
        value = Fraction(numerator, len(code))
        if minimum is None or value < minimum:
            minimum = value
            minimum_count = 1
        elif value == minimum:
            minimum_count += 1
        total += 1
    assert minimum <= rank_one_value
    size = sum(k * p for k, p in shapes)
    eta = Fraction(sum(k * p * (k + p) for k, p in shapes), size * size)
    return {"shapes": shapes, "physical_switching_classes": total,
            "code_size": len(code), "minimum_exact_mean": str(minimum),
            "independent_rank_one_mean": str(rank_one_value),
            "minimum_classes": minimum_count,
            "eta": str(eta), "rank_one_is_finite_minimizer": minimum == rank_one_value}


def spectral_audit():
    rng = np.random.default_rng(20260917)
    pointwise = 0
    moments = 0
    for q in range(2, 7):
        spins = np.asarray(words(q), dtype=float)
        for repetition in range(12):
            raw = rng.choice((-1, 1), size=(q, q))
            matrix = np.triu(raw, 1)
            matrix = matrix + matrix.T
            for polarity in (-1, 1):
                signed = polarity * matrix
                curvature = float(np.linalg.eigvalsh(signed)[-1])
                row = spins * (spins @ signed)
                expectation = np.maximum(curvature - row, 0) ** 2
                assert np.all(expectation.mean(axis=0) <= curvature ** 2 + q - 1 + 1e-10)
                moments += q
                for _ in range(8):
                    gaussian = rng.normal(size=q)
                    sigma = rng.uniform(0.05, 5, size=q)
                    z = np.where(gaussian >= 0, 1.0, -1.0)
                    baseline = 0.5 * z @ signed @ z + sigma @ np.abs(gaussian)
                    row_at_z = z * (signed @ z)
                    certificate = baseline + 2 * np.maximum(curvature - row_at_z - sigma * np.abs(gaussian), 0).sum()
                    true_value = np.max(0.5 * np.einsum('bi,ij,bj->b', spins, signed, spins) + spins @ (sigma * gaussian))
                    assert true_value <= certificate + 1e-10
                    pointwise += 1
    return {"pointwise_spectral_max_certificates": pointwise,
            "exact_sign_law_row_moment_checks": moments}


def main():
    shapes = [[(2, 2)], [(2, 3)], [(3, 3)], [(3, 4)], [(4, 4)],
              [(4, 5)], [(2, 2), (2, 2)], [(2, 3), (3, 2)],
              [(1, 2), (2, 2), (2, 1)]]
    cases = [exact_dictionary_case(item) for item in shapes]
    constants = {}
    for r in (1, 2, 3, 4, 8, 16, 64):
        beta = math.exp(math.log(2) + math.lgamma((r + 1) / 2)
                        - 0.5 * math.log(math.pi * r) - math.lgamma(r / 2))
        constants[str(r)] = {"beta": beta, "cap_threshold_two_beta_over_three": 2 * beta / 3}
    print(json.dumps({"status": "PASS", "exact_dictionary_cases": cases,
                      "spectral_audit": spectral_audit(), "equal_block_constants": constants},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
