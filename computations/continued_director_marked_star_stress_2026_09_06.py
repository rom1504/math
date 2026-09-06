#!/usr/bin/env python3
"""Test one precise missing centered-coefficient feedback energy term.

No convergence or asymptotic inference is made from these data. Integer
full-cube cap replay labels a stored order-14 signing; larger cases are
actual sign constructions, not asserted minimizers. The two independent
Monte Carlo samples avoid reusing the coefficient-fit sample as its test.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing
from continued_feedback_hidden_third_return_2026_09_06 import construct


def run(a: np.ndarray, name: str, samples: int, seed: int, exact: bool) -> dict:
    n = len(a)
    assert np.array_equal(a, a.T) and not np.any(a.diagonal())
    assert np.all(np.abs(a[np.triu_indices(n, 1)]) == 1)
    b = a / np.sqrt(n - 1)
    q = b @ b
    cov_d = (1 - 3 / (n - 1)) * np.eye(n) + 2 * q / (n - 1)
    gamma = np.diag(q @ b)
    rng = np.random.default_rng(seed)
    batch = 512

    def seed_batches():
        count = (1 << n) if exact else samples
        for start in range(0, count, batch):
            size = min(batch, count - start)
            if exact:
                labels = np.arange(start, start + size, dtype=np.uint64)
                yield 1.0 - 2.0 * ((labels[:, None] >> np.arange(n, dtype=np.uint64)) & 1)
            else:
                yield rng.choice((-1.0, 1.0), size=(size, n))

    def fields(s):
        g = s @ b
        d = s * (g * g - 1) / np.sqrt(2)
        y = d @ b
        v = y @ b
        z = (y**3 - 3 * y) @ b / np.sqrt(6)
        return g, d, y, v, z

    mean_cos = np.zeros(n)
    total = 0
    for s in seed_batches():
        v = fields(s)[3]
        mean_cos += np.cos(v).sum(axis=0)
        total += len(s)
    mean_cos /= total

    acc = np.zeros((4, 2))
    actual_cov_d = np.zeros((n, n)) if exact else None
    tested_cap = 0
    for s in seed_batches():
        g, d, y, v, z = fields(s)
        residual = v - d - gamma * y
        centered_cross = np.mean((np.sin(v) @ b) * (np.cos(v) - mean_cos) * z, axis=1)
        phi = np.sin(g) * np.cos(y)
        regression = np.mean(residual * phi, axis=1)
        data = np.stack((centered_cross, regression,
                         np.mean(residual * y, axis=1),
                         np.mean(residual * d, axis=1)))
        acc[:, 0] += data.sum(axis=1)
        acc[:, 1] += (data * data).sum(axis=1)
        if exact:
            actual_cov_d += d.T @ d
            sint = s.astype(np.int64)
            twice_energy = np.sum((sint @ a) * sint, axis=1)
            assert np.all(twice_energy % 2 == 0)
            tested_cap = max(tested_cap, int(np.max(np.abs(twice_energy)) // 2))
    means = acc[:, 0] / total
    se = np.sqrt(np.maximum(acc[:, 1] / total - means * means, 0) / total)
    names = ("centered_coefficient_cross", "bounded_old_test_regression",
             "Y_test_regression", "D_test_regression")
    result = {
        "status": "finite exact enumeration" if exact else "Monte Carlo falsification only",
        "family": name, "n": n, "samples_per_pass": total, "seed": seed,
        "operator_norm": float(np.linalg.norm(b, 2)),
        "gamma_mean": float(gamma.mean()), "gamma_mean_square": float(np.mean(gamma**2)),
        "statistics": {key: {"mean": float(value),
                             "conditional_test_se": 0.0 if exact else float(error)}
                       for key, value, error in zip(names, means, se)},
        "warning": "Monte Carlo SE excludes coefficient-pilot uncertainty; no asymptotic claim",
    }
    if exact:
        result["integer_replayed_cap"] = tested_cap
        result["exact_covariance_identity_max_float_error"] = float(
            np.max(np.abs(actual_cov_d / total - cov_d)))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=20000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    stored = json.loads(Path("computations/results/heuristic_m14_from_conference.json").read_text())
    cases = [(np.asarray(stored["matrix"], dtype=np.int64), "stored_order14_cap21", True)]
    cases += [(make_signing(k)[0], f"Steiner_parameter{k}", False) for k in (3, 4, 5)]
    cases += [(construct(m, 9062026 + m), f"hidden_third_return_{m}", False)
              for m in (4, 16, 64)]
    output = []
    for index, (a, name, exact) in enumerate(cases):
        record = run(a, name, args.samples, 901372 + index, exact)
        output.append(record)
        print(json.dumps(record), flush=True)
    if args.output:
        args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
