#!/usr/bin/env python3
"""Exact/numerical checks for hypercube_flatness_literature_scout.md.

The deterministic part verifies the two-block hollow signing, its complete
spectrum, its exposed eigenvector, and the limiting flatness formula.  The
random part is a seeded sanity check of the elementary Haar/Gaussian baseline;
it is not evidence for the imported Bourgade--Yau theorem.
"""

from __future__ import annotations

from math import cos, hypot, pi, sqrt

import numpy as np


def two_block_signing(m: int) -> np.ndarray:
    """Return the 2m by 2m signing from (HF.13)."""

    ones = np.ones((m, m), dtype=np.int64)
    eye = np.eye(m, dtype=np.int64)
    return np.block([[ones - eye, ones], [ones, -(ones - eye)]])


def verify_exact_block_algebra() -> int:
    checks = 0
    for m in range(2, 13):
        signing = two_block_signing(m)
        n = 2 * m

        assert np.array_equal(signing, signing.T)
        assert np.array_equal(np.diag(signing), np.zeros(n, dtype=np.int64))
        off_diagonal = signing[~np.eye(n, dtype=bool)]
        assert set(off_diagonal.tolist()) == {-1, 1}

        # The exact block quotient K satisfies K^2=r^2 I over the integers.
        quotient = np.asarray(((m - 1, m), (m, -(m - 1))), dtype=np.int64)
        radius_squared = (m - 1) ** 2 + m**2
        assert np.array_equal(
            quotient @ quotient,
            radius_squared * np.eye(2, dtype=np.int64),
        )

        radius = sqrt(radius_squared)
        expected = np.asarray(
            [-radius] + [-1.0] * (m - 1) + [1.0] * (m - 1) + [radius]
        )
        observed = np.linalg.eigvalsh(signing.astype(float))
        assert np.allclose(observed, expected, atol=2e-12, rtol=2e-12)

        ratio = (radius - (m - 1)) / m
        alpha = sqrt(2 / (1 + ratio**2))
        beta = ratio * alpha
        exposed = np.concatenate(
            (np.full(m, alpha), np.full(m, beta))
        )
        assert abs(float(exposed @ exposed) - n) < 2e-12
        residual = signing @ exposed - radius * exposed
        assert np.linalg.norm(residual) < 2e-11

        l1_ratio = float(np.abs(exposed).sum()) / n
        formula = (1 + ratio) / sqrt(2 * (1 + ratio**2))
        assert abs(l1_ratio - formula) < 2e-14
        checks += 10

    return checks


def verify_flatness_limit() -> tuple[int, float]:
    target = cos(pi / 8)
    previous_error = float("inf")
    checks = 0
    final_ratio = 0.0

    for exponent in range(1, 21):
        m = 2**exponent
        radius = hypot(m - 1, m)
        ratio = (radius - (m - 1)) / m
        l1_ratio = (1 + ratio) / sqrt(2 * (1 + ratio**2))
        error = abs(l1_ratio - target)
        assert error < previous_error
        previous_error = error
        final_ratio = l1_ratio
        checks += 1

    assert abs(final_ratio - target) < 2e-7
    assert abs((1 - final_ratio) - (1 - target)) < 2e-7
    return checks + 2, final_ratio


def verify_gaussian_baseline() -> tuple[int, float, float]:
    """Seeded numerical check of (HF.12), in bounded-memory batches."""

    rng = np.random.default_rng(20260817)
    dimension = 4096
    trials = 512
    batch = 32
    samples: list[np.ndarray] = []

    for start in range(0, trials, batch):
        count = min(batch, trials - start)
        gaussian = rng.normal(size=(count, dimension))
        # For u=sqrt(N)g/||g||, ||u||_1/N is the expression below.
        ratios = np.abs(gaussian).sum(axis=1) / (
            sqrt(dimension) * np.linalg.norm(gaussian, axis=1)
        )
        samples.append(ratios)

    values = np.concatenate(samples)
    mean = float(values.mean())
    standard_deviation = float(values.std(ddof=1))
    target = sqrt(2 / pi)

    # These loose deterministic thresholds only catch normalization/regression
    # errors.  They are not probabilistic claims in the accompanying note.
    assert abs(mean - target) < 1.5e-3
    assert standard_deviation < 1.2e-2
    assert abs((1 - mean) - (1 - target)) < 1.5e-3
    return 3, mean, standard_deviation


def main() -> None:
    checks = verify_exact_block_algebra()
    limit_checks, block_ratio = verify_flatness_limit()
    checks += limit_checks
    gaussian_checks, gaussian_mean, gaussian_sd = verify_gaussian_baseline()
    checks += gaussian_checks

    print(
        "two-block signing:",
        f"l1/n={block_ratio:.10f}",
        f"target={cos(pi / 8):.10f}",
    )
    print(
        "Gaussian/Haar baseline:",
        f"mean={gaussian_mean:.10f}",
        f"target={sqrt(2 / pi):.10f}",
        f"sample_sd={gaussian_sd:.10f}",
    )
    print(f"hypercube-flatness scout checks passed: {checks}")


if __name__ == "__main__":
    main()
