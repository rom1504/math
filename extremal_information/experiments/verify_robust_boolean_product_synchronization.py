#!/usr/bin/env python3
"""Finite checks for robust Boolean product synchronization.

The proof is analytic.  These checks diagnose its normalizations:

* Fourier Parseval and the selector witness;
* the sharp intrinsic and operator-defect response bounds;
* the PSD diagonal Cauchy--Schwarz bound;
* tensor Gram/Rayleigh identities and intrinsic subadditivity;
* Gram and indefinite-Rayleigh Schur contractions;
* projective-histogram response normalization;
* exact-sign completion Lipschitzness on a small instance.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import numpy as np


def words(p: int):
    return list(product((-1, 1), repeat=p))


def majority_tau(a):
    total = sum(a)
    if total:
        return 1 if total > 0 else -1
    return a[0]


def fourier_tau(p: int):
    cube = words(p)
    active = []
    for mask in range(1 << p):
        value = sum(
            majority_tau(a)
            * np.prod([a[i] for i in range(p) if mask & (1 << i)], dtype=int)
            for a in cube
        )
        coeff = Fraction(int(value), 1 << p)
        if coeff:
            active.append((mask, coeff))
    return active


def active_columns(w: np.ndarray, active):
    n, p = w.shape
    cols = []
    for mask, _ in active:
        col = np.ones(n, dtype=float)
        for i in range(p):
            if mask & (1 << i):
                col *= w[:, i]
        cols.append(col)
    return np.column_stack(cols)


def coefficient_vector(epsilon, active):
    out = []
    for mask, coeff in active:
        sign = np.prod(
            [epsilon[i] for i in range(len(epsilon)) if mask & (1 << i)],
            dtype=int,
        )
        out.append(float(coeff * int(sign)))
    return np.asarray(out)


def state(h: np.ndarray, w: np.ndarray, r: float, active):
    z = active_columns(w, active)
    n = len(h)
    g = z.T @ z / n
    rayleigh = z.T @ h @ z / (r * n)
    return z, g, rayleigh, g - rayleigh


def trust_response(h: np.ndarray, w: np.ndarray, m: int, epsilon) -> float:
    n = len(h)
    field = w @ np.asarray(epsilon, dtype=float)
    best = -float("inf")
    for x_tuple in words(n):
        x = np.asarray(x_tuple, dtype=float)
        child = float(x @ h @ x) / 2
        value = abs(child) + m * float(field @ x)
        best = max(best, value)
    return best


def antipodal_system(p: int, gamma: Fraction):
    cube = words(p)
    w = np.asarray(cube, dtype=float)
    index = {a: i for i, a in enumerate(cube)}
    j = np.zeros((len(cube), len(cube)), dtype=float)
    for i, a in enumerate(cube):
        j[i, index[tuple(-v for v in a)]] = 1
    h = -(1 - float(gamma)) * j
    return h, w


def projective_key(row):
    row = tuple(int(v) for v in row)
    neg = tuple(-v for v in row)
    return min(row, neg)


def histogram(w: np.ndarray):
    hist = {}
    for row in w:
        key = projective_key(row)
        hist[key] = hist.get(key, 0) + 1 / len(w)
    return hist


def response_metric(mu, nu, p: int):
    def expected(measure, epsilon):
        return sum(
            mass * abs(sum(s[i] * epsilon[i] for i in range(p))) / p
            for s, mass in measure.items()
        )

    return max(
        abs(expected(mu, epsilon) - expected(nu, epsilon))
        for epsilon in words(p)
    )


def check_selector_response_and_histogram() -> int:
    p = 3
    gamma = Fraction(1, 4)
    h, w = antipodal_system(p, gamma)
    n = len(h)
    r = 1.0  # a declared norm bound; ||H||=3/4
    m = 1
    active = fourier_tau(p)
    z, g, rayleigh, d = state(h, w, r, active)

    # Odd Boolean characters are orthogonal, so D=gamma I.
    assert np.allclose(g, np.eye(len(active)))
    assert np.allclose(d, float(gamma) * np.eye(len(active)))
    delta = np.linalg.norm(d, 2)
    assert np.isclose(delta, float(gamma))

    diagonal_bound = 0.0
    for (_, coeff), deficit in zip(active, np.diag(d)):
        diagonal_bound += abs(float(coeff)) * np.sqrt(deficit)
    diagonal_bound **= 2
    fourier_l1 = sum(abs(float(coeff)) for _, coeff in active)
    assert np.isclose(diagonal_bound, float(gamma) * fourier_l1**2)

    mu = histogram(w)
    first = next(iter(mu))
    nu = {first: 1.0}
    eta = response_metric(mu, nu, p)
    c = m * p / r

    checks = 0
    for epsilon in words(p):
        a = coefficient_vector(epsilon, active)
        assert np.isclose(a @ a, 1)
        selector = z @ a
        direct = np.asarray(
            [majority_tau(tuple(epsilon[i] * int(row[i]) for i in range(p))) for row in w]
        )
        assert np.allclose(selector, direct)
        assert set(np.rint(selector).astype(int)) <= {-1, 1}

        field = w @ np.asarray(epsilon)
        assert np.isclose(field @ selector, np.abs(field).sum())
        intrinsic = float(a @ d @ a)
        assert intrinsic <= delta + 1e-10
        assert intrinsic <= diagonal_bound + 1e-10

        actual = trust_response(h, w, m, epsilon)
        ideal = r * n / 2 + m * np.abs(field).sum()
        assert -1e-9 <= ideal - actual <= intrinsic * r * n / 2 + 1e-9
        # This benchmark saturates the robust loss exactly.
        assert np.isclose(ideal - actual, float(gamma) * n / 2)

        decoded_support = n * sum(
            mass * abs(sum(s[i] * epsilon[i] for i in range(p)))
            for s, mass in nu.items()
        )
        decoded = r * n / 2 + m * decoded_support
        assert abs(actual - decoded) <= r * n * (intrinsic / 2 + c * eta) + 1e-9
        checks += 8
    return checks


def check_tensor_laws() -> int:
    p = 3
    active = fourier_tau(p)
    h1, w1 = antipodal_system(p, Fraction(1, 5))
    h2, w2 = antipodal_system(p, Fraction(1, 3))
    z1, g1, r1, d1 = state(h1, w1, 1.0, active)
    z2, g2, r2, d2 = state(h2, w2, 1.0, active)

    h12 = np.kron(h1, h2)
    w12 = np.column_stack([np.kron(w1[:, i], w2[:, i]) for i in range(p)])
    z12, g12, r12, d12 = state(h12, w12, 1.0, active)
    assert np.allclose(z12, np.column_stack([np.kron(z1[:, j], z2[:, j]) for j in range(len(active))]))
    assert np.allclose(g12, g1 * g2)
    assert np.allclose(r12, r1 * r2)
    assert np.allclose(d12, d1 * g2 + r1 * d2)
    assert np.allclose(d12, d1 * g2 + g1 * d2 - d1 * d2)

    def intrinsic(d):
        return max(
            coefficient_vector(epsilon, active) @ d @ coefficient_vector(epsilon, active)
            for epsilon in words(p)
        )

    delta1 = np.linalg.norm(d1, 2)
    delta2 = np.linalg.norm(d2, 2)
    delta12 = np.linalg.norm(d12, 2)
    assert intrinsic(d12) <= intrinsic(d1) + intrinsic(d2) + 1e-10
    assert delta12 <= delta1 + delta2 + 1e-10

    kappa1 = max(np.linalg.norm(h1 @ z1[:, j]) / np.sqrt(len(h1)) for j in range(z1.shape[1]))
    kappa2 = max(np.linalg.norm(h2 @ z2[:, j]) / np.sqrt(len(h2)) for j in range(z2.shape[1]))
    assert delta12 <= delta1 + kappa1 * delta2 + 1e-10
    assert delta12 <= kappa2 * delta1 + delta2 + 1e-10
    return 9


def check_schur_contractions(seed: int = 20260817) -> int:
    rng = np.random.default_rng(seed)
    checks = 0
    for q in range(2, 9):
        dim = q + 2
        a = rng.normal(size=(q, dim))
        b = rng.normal(size=(q, dim))
        a /= np.maximum(1.0, np.linalg.norm(a, axis=1))[:, None]
        b /= np.maximum(1.0, np.linalg.norm(b, axis=1))[:, None]
        kernel = a @ b.T
        x = rng.normal(size=(q, q))
        assert np.linalg.norm(kernel * x, 2) <= np.linalg.norm(x, 2) + 1e-10
        checks += 1

        # A Rayleigh kernel can be indefinite and remains contractive.
        u = rng.normal(size=(q, dim))
        u /= np.linalg.norm(u, axis=1)[:, None]
        t = rng.normal(size=(dim, dim))
        t = (t + t.T) / 2
        t /= max(1.0, np.linalg.norm(t, 2))
        rayleigh = u @ t @ u.T
        eig = np.linalg.eigvalsh(rayleigh)
        if eig[0] < -1e-8 and eig[-1] > 1e-8:
            assert np.linalg.norm(rayleigh * x, 2) <= np.linalg.norm(x, 2) + 1e-10
            checks += 1
    return checks


def cap(matrix: np.ndarray) -> float:
    best = 0.0
    for x_tuple in words(len(matrix)):
        x = np.asarray(x_tuple, dtype=float)
        best = max(best, abs(float(x @ matrix @ x) / 2))
    return best


def check_completion_lipschitz() -> int:
    # Small exact-sign old block, two one-vertex shores, and one completion edge.
    h = np.asarray(
        [
            [1, 1, 1, -1],
            [1, 1, -1, 1],
            [1, -1, -1, -1],
            [-1, 1, -1, -1],
        ],
        dtype=float,
    )
    assert np.allclose(h, h.T)
    assert np.trace(h) == 0
    old = h - np.diag(np.diag(h))
    assert set(old[np.triu_indices(4, 1)].astype(int)) <= {-1, 1}
    w = np.asarray([[1, 1], [1, -1], [-1, 1], [-1, -1]], dtype=float)
    bridge = w
    incomplete = np.block(
        [[old, bridge], [bridge.T, np.zeros((2, 2), dtype=float)]]
    )
    completion = np.asarray([[0, -1], [-1, 0]], dtype=float)
    completed = incomplete.copy()
    completed[4:, 4:] = completion

    labelled = max(trust_response(h, w, 1, epsilon) for epsilon in words(2))
    assert np.isclose(cap(incomplete), labelled)
    assert abs(cap(completed) - cap(incomplete)) <= cap(completion) + 1e-10
    return 4


def main() -> None:
    checks = 0
    checks += check_selector_response_and_histogram()
    checks += check_tensor_laws()
    checks += check_schur_contractions()
    checks += check_completion_lipschitz()
    print(f"robust Boolean product-synchronization checks passed: {checks}")


if __name__ == "__main__":
    main()
