#!/usr/bin/env python3
"""Finite checks of the Schur/flat-seed/twin identities; no optimizer claims.

Run with the repository's .venv/bin/python. All matrices, seeds, and test
counts are deterministic. Integer identities are checked exactly; scalar
Gaussian/Cholesky identities are numerical checks of the accompanying proof.
"""

from __future__ import annotations

import json
import math

import numpy as np


def g(t: float, v: np.ndarray) -> np.ndarray:
    z = 2.0 * t * np.asarray(v)
    rho = 2.0 * z / (1.0 + np.sqrt(1.0 + 4.0 * z * z))
    return -t * v * (1.0 - rho) + 0.25 * np.log1p(-rho * rho)


def c(t: float, precision: np.ndarray) -> np.ndarray:
    u = precision / t
    return 0.25 * np.log(u * (2.0 - u))


def all_spins(n: int) -> np.ndarray:
    # Global reversal preserves all tested quadratic quantities.
    ids = np.arange(1 << (n - 1), dtype=np.int64)
    out = np.ones((len(ids), n), dtype=np.int64)
    out[:, 1:] = 1 - 2 * ((ids[:, None] >> np.arange(n - 1)) & 1)
    return out


def cap(a: np.ndarray, spins: np.ndarray) -> int:
    values = np.einsum("bi,ij,bj->b", spins, a, spins)
    assert np.all(values % 2 == 0)
    return int(np.max(np.abs(values)) // 2)


def hollow(rng: np.random.Generator, n: int) -> np.ndarray:
    a = np.triu(rng.choice([-1, 1], size=(n, n)), 1)
    return (a + a.T).astype(np.int64)


def twin(a: np.ndarray, pairs: list[tuple[int, int]], signs: np.ndarray):
    n = len(a)
    representatives = np.arange(n)
    for r, s in pairs:
        representatives[s] = r
    b = a[np.ix_(representatives, representatives)].copy()
    for (r, s), sign in zip(pairs, signs):
        b[r, s] = b[s, r] = sign
    return b, representatives


def main() -> None:
    rng = np.random.default_rng(20260906)
    schur_tests = 0
    least_schur_gap = math.inf
    least_diagonal_gap = math.inf
    for d in range(2, 10):
        for _ in range(30):
            r = rng.normal(size=(d, d))
            sigma = r @ r.T + 0.03 * np.eye(d)
            o, _ = np.linalg.qr(rng.normal(size=(d, d)))
            eigenvalues = np.linalg.eigvalsh(sigma)
            pivots = np.diag(np.linalg.cholesky(sigma)) ** 2
            rotated_diagonal = np.diag(o @ sigma @ o.T)
            log_p = np.sort(np.log(pivots))[::-1]
            log_e = np.sort(np.log(eigenvalues))[::-1]
            assert np.max(np.cumsum(log_p) - np.cumsum(log_e)) < 2e-10
            t = float(np.exp(rng.uniform(-2, 2)))
            gap_one = float(np.sum(g(t, pivots)) - np.sum(g(t, eigenvalues)))
            gap_two = float(np.sum(g(t, eigenvalues)) - np.sum(g(t, rotated_diagonal)))
            assert gap_one > -2e-10 and gap_two > -2e-10
            least_schur_gap = min(least_schur_gap, gap_one)
            least_diagonal_gap = min(least_diagonal_gap, gap_two)
            schur_tests += 1

    twin_results = []
    for n, count in [(6, 1), (10, 2), (12, 3)]:
        a = hollow(rng, n)
        pairs = [(2 * j, 2 * j + 1) for j in range(count)]
        signs = rng.choice([-1, 1], size=count)
        b, representatives = twin(a, pairs, signs)
        spins = all_spins(n)
        jset = [j for pair in pairs for j in pair]
        row_norm = int(np.max(np.sum(np.abs(spins @ a[:, jset]), axis=1)))
        qa, qb = cap(a, spins), cap(b, spins)
        assert qb <= qa + 1.5 * row_norm + count

        # Exact polynomial identity and exact small-eigenvector identity.
        for x in spins:
            y = np.zeros(n, dtype=np.int64)
            np.add.at(y, representatives, x)
            pair_term = sum(int(s) * int(x[r]) * int(x[v])
                            for (r, v), s in zip(pairs, signs))
            assert int(x @ b @ x) == int(y @ a @ y) + 2 * pair_term
        for (r, s), sign in zip(pairs, signs):
            v = np.zeros(n, dtype=np.int64)
            v[r], v[s] = 1, -1
            assert np.array_equal(b @ v, -sign * v)

        tested_completions = 0
        precision_checks = 0
        max_determinant_slack = -math.inf
        for _ in range(60):
            full = b + np.diag(rng.choice([-1, 1], size=n))
            singular = np.linalg.svd(full, compute_uv=False)
            assert singular[-count] <= 2 + 1e-9
            gram = full.T @ full / n
            sign, logdet = np.linalg.slogdet(gram)
            delta = count / n
            bound = delta * math.log(4 / n) - (1 - delta) * math.log(1 - delta)
            if singular[-1] > 1e-7:
                assert sign > 0 and logdet / n <= bound + 1e-9
                max_determinant_slack = max(max_determinant_slack, logdet / n - bound)
                for _ in range(3):
                    t = float(np.exp(rng.uniform(-1, 1)))
                    lam = t * rng.uniform(0.05, 1, size=n)
                    order = rng.permutation(n)
                    m = full[:, order] / math.sqrt(n)
                    p = m.T @ (lam[:, None] * m)
                    pivots = np.diag(np.linalg.cholesky(p)) ** 2
                    assert np.max(pivots) <= t + 1e-8
                    penalty = float(np.sum(c(t, lam)) - np.sum(c(t, pivots)))
                    identity = -0.25 * logdet + 0.25 * np.sum(
                        np.log((2 - lam / t) / (2 - pivots / t)))
                    assert abs(penalty - identity) < 2e-7
                    assert penalty >= -0.25 * logdet - n * math.log(2) / 4 - 2e-7
                    precision_checks += 1
            tested_completions += 1
        twin_results.append({
            "order": n, "twin_pairs": count, "original_cap": qa,
            "twin_cap": qb, "chosen_row_bilinear_norm": row_norm,
            "all_quadratic_identity_checks": len(spins),
            "completion_tests": tested_completions,
            "nonsingular_precision_checks": precision_checks,
            "largest_logdet_minus_bound": max_determinant_slack,
        })

    # Exact sign-flat nonorthogonal reversed butterfly and Gram identity.
    d, q = 6, 4
    seed = hollow(rng, d) + np.diag(rng.choice([-1, 1], size=d))
    h4 = np.ones((q, q), dtype=np.int64) - 2 * np.eye(q, dtype=np.int64)
    children = np.zeros((d * q, d * q), dtype=np.int64)
    for j in range(d):
        child = h4[rng.permutation(q)] * rng.choice([-1, 1], size=(q, 1))
        children[j * q:(j + 1) * q, j * q:(j + 1) * q] = child
    raw = np.kron(seed, np.eye(q, dtype=np.int64)) @ children
    assert np.all(np.abs(raw) == 1)
    assert np.array_equal(raw.T @ raw,
                          children.T @ np.kron(seed.T @ seed, np.eye(q, dtype=np.int64)) @ children)

    print(json.dumps({
        "status": "all finite identities and inequalities passed",
        "random_seed": 20260906,
        "schur_tests": schur_tests,
        "least_cholesky_reward_gap": least_schur_gap,
        "least_rotated_diagonal_reward_gap": least_diagonal_gap,
        "twin_tests": twin_results,
        "exact_flat_base_order": d * q,
        "scope": "finite verification, not evidence of minimax convergence",
    }, indent=2))


if __name__ == "__main__":
    main()
