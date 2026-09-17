"""Finite exact-branching checks of biased Gram--Schmidt walk statements.

The branching probabilities are computed in floating point; all final
colorings are exhaustively enumerated for the small input. This is a
falsification test, not a proof of the imported theorem.
"""

from __future__ import annotations

import json

import numpy as np


def projection(columns: np.ndarray) -> np.ndarray:
    if columns.shape[1] == 0:
        return np.zeros((columns.shape[0], columns.shape[0]))
    return columns @ np.linalg.pinv(columns)


def enumerate_walk(v: np.ndarray, initial: np.ndarray):
    """Return (probability, final sign, random phase compensator) leaves."""
    n = len(initial)
    leaves = []

    def visit(z, probability, gamma, pivot=None, phase_start=None):
        active = np.flatnonzero(np.abs(z) < 1.0 - 1e-9)
        if pivot is not None and pivot not in active:
            difference = phase_start - projection(v[:, active])
            difference = (difference + difference.T) / 2.0
            residual = difference @ v[:, pivot]
            gamma = gamma + float(residual @ residual) * difference
            pivot, phase_start = None, None
        if len(active) == 0:
            leaves.append((probability, np.sign(z), gamma))
            return
        if pivot is None:
            pivot = int(active[-1])
            phase_start = projection(v[:, active])
        others = active[active != pivot]
        direction = np.zeros(n)
        direction[pivot] = 1.0
        if len(others):
            direction[others] = -np.linalg.lstsq(v[:, others], v[:, pivot], rcond=None)[0]
        lower, upper = -np.inf, np.inf
        for i in active:
            if abs(direction[i]) < 1e-12:
                continue
            a = (-1.0 - z[i]) / direction[i]
            b = (1.0 - z[i]) / direction[i]
            lower = max(lower, min(a, b))
            upper = min(upper, max(a, b))
        assert lower < 0.0 < upper
        for delta, prob in [(lower, upper / (upper - lower)),
                            (upper, -lower / (upper - lower))]:
            updated = np.clip(z + delta * direction, -1.0, 1.0)
            visit(updated, probability * prob, gamma, pivot, phase_start)

    visit(initial.copy(), 1.0, np.zeros((v.shape[0], v.shape[0])))
    return leaves


def run_case(v, initial, rng):
    leaves = enumerate_walk(v, initial)
    probabilities = np.array([row[0] for row in leaves])
    signs = np.array([row[1] for row in leaves])
    gamma = np.array([row[2] for row in leaves])
    assert abs(probabilities.sum() - 1.0) < 1e-10
    mean_error = float(np.max(np.abs(probabilities @ signs - initial)))
    assert mean_error < 1e-9
    error_vectors = (signs - initial) @ v.T
    gamma_mean = np.einsum("k,kij->ij", probabilities, gamma)
    largest_gamma_eigenvalue = max(float(np.linalg.eigvalsh(g)[-1]) for g in gamma)
    assert largest_gamma_eigenvalue <= 1.0 + 1e-8
    worst_sharp = -np.inf
    worst_self_normalized = -np.inf
    worst_averaged = -np.inf
    counterexample = None
    for _ in range(100):
        direction = rng.normal(size=v.shape[0])
        direction *= rng.uniform(0.1, 10.0) / np.linalg.norm(direction)
        linear = error_vectors @ direction
        quadratic = np.einsum("i,kij,j->k", direction, gamma, direction)
        log_mgf = np.logaddexp.reduce(np.log(probabilities) + linear)
        sharp_difference = log_mgf - float(direction @ direction) / 2.0
        self_normalized = np.logaddexp.reduce(np.log(probabilities) + linear - quadratic / 2.0)
        averaged_difference = log_mgf - float(direction @ gamma_mean @ direction) / 2.0
        worst_sharp = max(worst_sharp, sharp_difference)
        worst_self_normalized = max(worst_self_normalized, self_normalized)
        worst_averaged = max(worst_averaged, averaged_difference)
        if averaged_difference > 1e-7 and counterexample is None:
            counterexample = {
                "v": v.tolist(), "initial": initial.tolist(),
                "theta": direction.tolist(), "log_mgf": float(log_mgf),
                "half_mean_compensator": float(direction @ gamma_mean @ direction) / 2.0,
                "gap": float(averaged_difference),
            }
    assert worst_sharp <= 1e-8
    assert worst_self_normalized <= 1e-8
    return {
        "leaves": len(leaves), "mean_error": mean_error,
        "max_gamma_eigenvalue": largest_gamma_eigenvalue,
        "max_log_mgf_minus_sharp_proxy": float(worst_sharp),
        "max_log_self_normalized_mgf": float(worst_self_normalized),
        "max_log_mgf_minus_averaged_phase_proxy": float(worst_averaged),
        "averaged_proxy_counterexample": counterexample,
    }


def main():
    rng = np.random.default_rng(202609171)
    cases = []
    counterexample = None
    for _ in range(100):
        n = 4
        v = rng.normal(size=(n, n))
        v /= np.linalg.norm(v, axis=0)
        initial = rng.uniform(-0.95, 0.95, size=n)
        case = run_case(v, initial, rng)
        if case["averaged_proxy_counterexample"] is not None:
            counterexample = case["averaged_proxy_counterexample"]
        cases.append(case)
        if counterexample is not None:
            break
    report = {
        "status": "all sharp and random-compensator tests passed",
        "cases_checked": len(cases),
        "largest_mean_error": max(row["mean_error"] for row in cases),
        "worst_sharp_gap": max(row["max_log_mgf_minus_sharp_proxy"] for row in cases),
        "worst_compensated_gap": max(row["max_log_self_normalized_mgf"] for row in cases),
        "averaged_proxy_counterexample": counterexample,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
