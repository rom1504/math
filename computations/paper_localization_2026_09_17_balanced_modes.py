#!/usr/bin/env python3
"""Finite exact-law diagnostics for deterministic balanced mode bridges.

Integer covariance identities and finite-law expectations are exact up to
the final floating-point norm evaluation. MGF and asymptotic entropy values
are numerical diagnostics only; the canonical artifact contains the proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp


def hadamard(k: int) -> np.ndarray:
    h = np.ones((1, 1), dtype=np.int64)
    while h.shape[0] < k:
        h = np.block([[h, h], [h, -h]])
    assert h.shape == (k, k)
    assert np.array_equal(h @ h.T, k * np.eye(k, dtype=np.int64))
    return h


def words(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def entropy(t: float) -> float:
    if t == 0 or t == 1:
        return 0.0
    return -t * math.log(t) - (1 - t) * math.log1p(-t)


def finite_case(k: int, p: int, ell: int, q: int,
                rng: np.random.Generator) -> dict:
    h = hadamard(k)
    n0, n = k * p, k * p + ell
    counts = np.array([q // k + int(a < q % k) for a in range(k)], dtype=np.int64)
    mode_gram = sum((int(counts[a]) * np.outer(h[a], h[a]) for a in range(k)),
                    np.zeros((k, k), dtype=np.int64))
    assert np.array_equal(h @ mode_gram @ h.T, k * k * np.diag(counts))
    assert np.linalg.eigvalsh(mode_gram).max() <= q + k + 1e-10
    if q % k == 0:
        assert np.array_equal(mode_gram, q * np.eye(k, dtype=np.int64))
    atoms = words(p + ell)
    largest_mgf_excess = -math.inf
    largest_mean_excess = -math.inf
    maximum_effective_modes = 0.0
    for trial in range(18):
        if trial < k:
            x0 = np.outer(h[trial], rng.choice((-1, 1), size=p))
        else:
            x0 = rng.choice((-1, 1), size=(k, p))
        leftovers = rng.choice((-1, 1), size=ell)
        coefficients = h @ x0
        energies = np.sum(coefficients ** 2, axis=1)
        assert int(energies.sum()) == k * n0
        effective = float(np.sqrt(energies).sum() ** 2 / (k * n0))
        assert 1 - 1e-12 <= effective <= k + 1e-12
        maximum_effective_modes = max(maximum_effective_modes, effective)
        means, log_mgfs = [], []
        variance_proxy = float(counts @ energies + q * ell)
        assert variance_proxy <= (q + k) * n + 1e-9
        t = float(rng.uniform(-4, 4) / math.sqrt(max(variance_proxy, 1)))
        for a in range(k):
            coeff = np.concatenate((coefficients[a], leftovers))
            responses = atoms @ coeff
            absolute = np.abs(responses)
            mean = float(absolute.mean())
            means.append(mean)
            log_mgfs.append(float(logsumexp(t * (absolute - mean)) - math.log(len(atoms))))
        bridge_mean = float(counts @ means)
        upper_mean = (q + k) * math.sqrt(n0 * effective / k) + q * math.sqrt(ell)
        mean_excess = bridge_mean - upper_mean
        assert mean_excess <= 1e-9
        largest_mean_excess = max(largest_mean_excess, mean_excess)
        mgf_excess = float(counts @ log_mgfs - t * t * variance_proxy / 2)
        assert mgf_excess <= 1e-10, (k, p, ell, q, trial, mgf_excess)
        largest_mgf_excess = max(largest_mgf_excess, mgf_excess)
    return {"k": k, "p": p, "ell": ell, "n": n, "q": q,
            "mode_counts": counts.tolist(), "finite_atoms_per_column": len(atoms),
            "is_exactly_isotropic_in_aggregate": q % k == 0,
            "largest_centered_absolute_mgf_excess": largest_mgf_excess,
            "largest_mean_bound_excess": largest_mean_excess,
            "largest_effective_mode_count_tested": maximum_effective_modes}


def entropy_bound_examples() -> list[dict]:
    result = []
    for k in (64, 256, 1024, 4096, 16384, 65536):
        p = k
        u = max(1, int(k ** 0.25))
        s = math.ceil(math.sqrt(u * k))
        error = u / s + s / k
        if error >= 0.5:
            continue
        log_choose = float(gammaln(k + 1) - gammaln(s + 1) - gammaln(k - s + 1))
        normalized_bound = (log_choose / (k * p)
                            + (s / (2 * k)) * math.log1p(k / s)
                            + entropy(error))
        result.append({"k": k, "p": p, "cutoff_u": u, "top_modes_s": s,
                       "prediction_error_bound": error,
                       "normalized_log_cardinality_bound": normalized_bound,
                       "normalized_center_mean_coefficient_bound": math.sqrt(u / k)})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/balanced_modes_audit.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(2026091704)
    cases = [finite_case(k, p, ell, q, rng)
             for k, p, ell in [(2, 2, 0), (2, 3, 1), (4, 2, 0),
                               (4, 3, 3), (8, 2, 0), (8, 2, 3)]
             for q in [1, k - 1, k, k + 1, 2 * k + 3]]
    result = {"status": "PASS", "seed": 2026091704,
              "finite_cases": cases, "entropy_bound_examples": entropy_bound_examples(),
              "scope": "Finite-law diagnostics; analytic theorem and all-order limits are proved separately."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "finite_cases": len(cases),
                      "centered_mgf_tests": 18 * len(cases),
                      "entropy_bound_examples": result["entropy_bound_examples"]}, indent=2))


if __name__ == "__main__":
    main()
