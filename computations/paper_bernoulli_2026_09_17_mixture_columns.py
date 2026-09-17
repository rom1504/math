"""Checks for the explicit finite-mixture tensor-column extension lemma.

Exact finite-law identities and limiting conditional Gaussian means are
tested separately. No near-extreme cover realization is claimed.
"""

import itertools
import json
import math
from fractions import Fraction

import numpy as np
from scipy.special import logsumexp

from paper_bernoulli_2026_09_17_basis_columns import hadamard


def main():
    rng = np.random.default_rng(2026091704)
    k, p = 4, 4
    n = k * p
    hk, hp = hadamard(k), hadamard(p)
    inner_signs = np.array(list(itertools.product((-1, 1), repeat=p)))
    base_columns = np.array([np.kron(row, g) for row in hk for g in inner_signs])
    base_centers = np.kron(hk, hp)
    permutation = rng.permutation(n)
    switches = rng.choice((-1, 1), size=n)
    inverse_permutation = np.argsort(permutation)

    columns = np.concatenate((base_columns, base_columns[:, permutation] * switches))
    centers = np.concatenate((base_centers, base_centers[:, permutation] * switches))
    assert np.all(columns.sum(axis=0) == 0)
    assert np.array_equal(columns.T @ columns, len(columns) * np.eye(n, dtype=int))
    exact_means = [Fraction(int(np.abs(columns @ f).sum()), len(columns)) for f in centers]

    # Repeating every physical position t times gives an exact sequence
    # whose conditional Gaussian limiting means can be read from the
    # finitely many mode coefficient vectors. Both groupings are repeated.
    gaussian_means = []
    kappa = math.sqrt(2 / math.pi)
    for f in centers:
        original_coordinates = [f, (f * switches)[inverse_permutation]]
        means_by_frame = []
        for old_f in original_coordinates:
            coefficients = hk @ old_f.reshape(k, p)
            assert np.max(np.abs(coefficients)) <= k
            assert np.sum(coefficients ** 2) == k * n
            conditional_mean = kappa * np.mean(np.linalg.norm(coefficients, axis=1))
            assert conditional_mean <= kappa * math.sqrt(n) + 1e-12
            means_by_frame.append(conditional_mean)
        gaussian_means.append(sum(means_by_frame) / 2 / math.sqrt(n))
    theorem_coefficient = 3 * kappa / 4
    assert max(gaussian_means) <= theorem_coefficient + 1e-12

    for _ in range(2000):
        theta = rng.normal(size=n) * 10 ** rng.uniform(-2, 0.5)
        log_mgf = logsumexp(columns @ theta) - math.log(len(columns))
        assert log_mgf <= k * np.dot(theta, theta) / 2 + 1e-11

    # Deterministic balancing acts on the SUM of the conditional column
    # Grams. The two physical coordinate groupings need not agree.
    balanced_checks = []
    for q in range(1, 65):
        frame_counts = (q // 2, q - q // 2)
        aggregate_gram = np.zeros((n, n), dtype=int)
        for frame, frame_count in enumerate(frame_counts):
            mode_counts = np.full(k, frame_count // k, dtype=int)
            mode_counts[:frame_count % k] += 1
            assert mode_counts.sum() == frame_count
            for mode, multiplicity in enumerate(mode_counts):
                linear_map = np.kron(hk[mode, :, None], np.eye(p, dtype=int))
                if frame:
                    linear_map = linear_map[permutation] * switches[:, None]
                aggregate_gram += int(multiplicity) * (linear_map @ linear_map.T)
        assert np.array_equal(np.diag(aggregate_gram), np.full(n, q))
        largest = float(np.linalg.eigvalsh(aggregate_gram)[-1])
        assert largest <= q + 2 * k + 1e-10
        if q % (2 * k) == 0:
            assert np.array_equal(aggregate_gram, q * np.eye(n, dtype=int))
        balanced_checks.append(largest - q)

    # Independent rational verification of the fourth-moment exception
    # counting algebra used in the two-frame applicability warning.
    c4 = Fraction(11, 27)
    for edges in range(1, 101):
        exceptions = edges // 2
        direct = edges + c4 * edges * (edges - 1) + 2 * (1 - c4) * exceptions
        assert direct <= c4 * edges ** 2 + 2 * (1 - c4) * edges

    slopes = [{"frames": r, "mode_order": kk,
               "zero_cover_slope": kappa * (1 - (1 - 1 / math.sqrt(kk)) / r)}
              for r, kk in ((2, 4), (3, 8), (4, 16), (5, 256))]
    assert all(item["zero_cover_slope"] < 1.5 * 0.4333221116640807 for item in slopes)
    print(json.dumps({
        "status": "PASS",
        "n": n,
        "mixture_law_atoms_with_multiplicity": len(columns),
        "center_count_with_multiplicity": len(centers),
        "covariance": "I exactly",
        "largest_exact_finite_center_mean": str(max(exact_means)),
        "largest_replicated_limit_normalized_mean": max(gaussian_means),
        "theorem_coefficient": theorem_coefficient,
        "mgf_checks": 2000,
        "balanced_aggregate_gram_checks": len(balanced_checks),
        "largest_balanced_gram_excess_over_q": max(balanced_checks),
        "mean_response_slopes": slopes,
        "scope": "Explicit column laws and conditional cover theorem only."
    }, indent=2))


if __name__ == "__main__":
    main()
