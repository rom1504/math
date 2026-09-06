"""Exact signed-coefficient tensor-rank identity tests; no Boolean closure claim."""
import json
from pathlib import Path
import numpy as np


def tensor_power(v, power):
    answer = np.array([1], dtype=np.int64)
    for _ in range(power):
        answer = np.kron(answer, v)
    return answer


def main():
    rng = np.random.default_rng(26090651)
    checks = 0
    mixed_sign_cases = 0
    zero_diagonal_cases = 0
    minimum_rank_slack = None
    for n in range(1, 9):
        for repetition in range(20):
            vectors = rng.integers(-2, 3, size=(n, 3), dtype=np.int64)
            if repetition % 4 == 0:
                vectors[0] = 0
                zero_diagonal_cases += 1
            weights = rng.integers(-3, 4, size=n, dtype=np.int64)
            mixed_sign_cases += bool(np.any(weights < 0) and np.any(weights > 0))
            gram = vectors @ vectors.T
            extra = rng.integers(-2, 3, size=(n, 2), dtype=np.int64)
            positive_increment = extra @ extra.T
            for power in range(1, 5):
                odd_degree = 2*power-1
                columns = np.asarray([tensor_power(v, power) for v in vectors])
                operator = columns.T @ (weights[:, None]*columns)
                frobenius_square = int((operator*operator).sum())
                trace = int(np.trace(operator))
                predicted_trace = int((weights * np.diag(gram)**power).sum())
                schur_congruence = weights[:, None] * gram**odd_degree * weights[None, :]
                kernel_trace = int((gram * schur_congruence.T).sum())
                assert trace == predicted_trace
                assert kernel_trace == frobenius_square
                slack = n*frobenius_square-trace*trace
                assert slack >= 0
                # Signed diagonal congruence of the odd Schur power is PSD;
                # test the exact trace comparison against a PSD increment.
                assert int((positive_increment*schur_congruence.T).sum()) >= 0
                minimum_rank_slack = slack if minimum_rank_slack is None else min(minimum_rank_slack, slack)
                checks += 1
    result = {
        "scope": "Exact finite algebra only; no old-field Boolean mixed-covariance comparison",
        "seed": 26090651, "tensor_rank_checks": checks,
        "vector_dimensions": 3, "numbers_of_rows": [1, 8],
        "odd_degrees": [1, 3, 5, 7],
        "mixed_signed_diagonal_inputs": mixed_sign_cases,
        "inputs_with_zero_Gram_diagonal": zero_diagonal_cases,
        "minimum_nonnegative_rank_slack": minimum_rank_slack,
    }
    path = Path(__file__).resolve().parent / "results/transfer_adversary_gaussian_return_rank_2026_09_06.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))
    print("PASS: signed tensor-rank identity, rank inequality, and PSD trace comparison")


if __name__ == "__main__":
    main()
