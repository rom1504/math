#!/usr/bin/env python3
"""Replay for the common-Gibbs fourth-order frame comparison.

The trapezoid normalization and balanced fourth-power sums are exact.
Gaussian expectations use deterministic Gauss--Hermite quadrature and are
diagnostics rather than certified integration or theorem premises.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp, roots_hermitenorm
import sympy as sp


def words(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def hadamard(k: int) -> np.ndarray:
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < k:
        h = np.block([[h, h], [h, -h]])
    return h


def exact_trapezoid() -> dict:
    u = sp.symbols("u", positive=True)
    s = sp.symbols("s", real=True)
    for degree in range(9):
        f = s ** degree
        lhs = sp.integrate(f, (s, -u, u)) / (2 * u) - (f.subs(s, -u) + f.subs(s, u)) / 2
        rhs = -sp.integrate((u ** 2 - s ** 2) * sp.diff(f, s, 2), (s, -u, u)) / (4 * u)
        assert sp.simplify(lhs - rhs) == 0
    derivative_kernel_mass = sp.simplify(sp.integrate(u ** 2 - s ** 2, (s, -u, u)) / (8 * u))
    assert derivative_kernel_mass == u ** 2 / 6
    return {"degrees_checked": list(range(9)), "derivative_kernel_mass": str(derivative_kernel_mass),
            "integrated_cumulant_constant": str(sp.Rational(19, 12))}


def exact_balanced_fourth() -> list[dict]:
    records = []
    for k, p, q in [(2, 3, 4), (4, 2, 8), (4, 2, 7), (4, 3, 9)]:
        h, n = hadamard(k), k * p
        counts = np.array([q // k + int(a < q % k) for a in range(k)], dtype=np.int64)
        x = words(n).reshape(-1, k, p)
        coeff = np.einsum("ab,ibp->iap", h, x)
        second = np.einsum("a,iap->i", counts, coeff ** 2)
        fourth = np.einsum("a,iap->i", counts, coeff ** 4)
        assert np.all(fourth <= k * k * second)
        assert np.all(second <= (q + k) * n)
        if q % k == 0:
            assert np.all(second == q * n)
            assert int(fourth.max()) == k * k * q * n
        records.append({"k": k, "p": p, "q": q, "queries_enumerated": len(x),
                        "exact_max_second_sum": int(second.max()),
                        "exact_max_fourth_sum": int(fourth.max()),
                        "claimed_fourth_upper_bound": k * k * (q + k) * n})
    return records


def quadrature_audit(rng: np.random.Generator) -> list[dict]:
    d, count = 3, 7
    coeff = rng.integers(-2, 3, size=(count, d))
    offsets = rng.normal(size=count)
    eps = words(d)
    nodes, weights = roots_hermitenorm(28)
    weights = weights / math.sqrt(2 * math.pi)
    indices = np.array(list(itertools.product(range(len(nodes)), repeat=d)), dtype=np.int64)
    gauss = nodes[indices]
    gauss_weights = np.prod(weights[indices], axis=1)
    b = int(np.abs(coeff).max())
    fourth = int(np.sum(coeff ** 4, axis=1).max())
    records = []
    for tau in [.04, .08, .16, .32]:
        sign_values = logsumexp(tau * (offsets + eps @ coeff.T), axis=1) / tau
        gaussian_values = logsumexp(tau * (offsets + gauss @ coeff.T), axis=1) / tau
        difference = float(sign_values.mean() - gauss_weights @ gaussian_values)
        bound = 19 / 12 * tau ** 3 * math.exp(4 * tau * b) * fourth
        assert abs(difference) <= bound + 1e-10
        records.append({"tau": tau, "B": b, "W4": fourth,
                        "softmax_sign_minus_gaussian_quadrature": difference,
                        "analytic_softmax_comparison_bound": bound,
                        "absolute_gap_over_bound": abs(difference) / bound})
    for _ in range(1000):
        a = rng.normal(size=count)
        probabilities = rng.dirichlet(np.ones(count))
        mean = float(probabilities @ a)
        centered = a - mean
        variance = float(probabilities @ centered ** 2)
        cumulant = float(probabilities @ centered ** 4 - 3 * variance ** 2)
        raw_fourth = float(probabilities @ a ** 4)
        assert abs(cumulant) <= 19 * raw_fourth + 1e-10
    return records


def truncation_audit(rng: np.random.Generator) -> dict:
    constant = 5 + math.sqrt(2) + 19 * math.exp(2) / 96
    assert constant < 8
    maximum_tail_l1_ratio = 0.0
    maximum_tail_l2_ratio = 0.0
    for _ in range(1000):
        count, dimension = 9, 6
        coefficients = rng.standard_t(2, size=(count, dimension))
        fourth = float(np.sum(coefficients ** 4, axis=1).max())
        ll = math.log(2 * count)
        cutoff = (fourth / ll) ** .25
        tail = np.where(np.abs(coefficients) > cutoff, coefficients, 0)
        l1 = float(np.abs(tail).sum(axis=1).max())
        l2 = float(np.sum(tail ** 2, axis=1).max())
        assert l1 <= fourth / cutoff ** 3 + 1e-9
        assert l2 <= fourth / cutoff ** 2 + 1e-9
        maximum_tail_l1_ratio = max(maximum_tail_l1_ratio, l1 / (fourth / cutoff ** 3))
        maximum_tail_l2_ratio = max(maximum_tail_l2_ratio, l2 / (fourth / cutoff ** 2))
    return {"queries_tested": 9000, "explicit_universal_constant": constant,
            "maximum_l1_tail_bound_ratio": maximum_tail_l1_ratio,
            "maximum_l2_squared_tail_bound_ratio": maximum_tail_l2_ratio}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/symmetric_frame_audit.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(2026091705)
    result = {"status": "PASS", "seed": 2026091705,
              "trapezoid": exact_trapezoid(),
              "balanced_fourth_sums": exact_balanced_fourth(),
              "quadrature_diagnostics": quadrature_audit(rng),
              "querywise_truncation": truncation_audit(rng),
              "scope": "Exact algebra plus non-certified Gaussian quadrature diagnostics; proof is analytic."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
