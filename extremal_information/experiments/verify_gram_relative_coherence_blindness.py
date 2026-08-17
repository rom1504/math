#!/usr/bin/env python3
"""Exact and numerical checks for Gram-relative coherence blindness.

The analytic proof is in
``drafts/gram_relative_coherence_blindness.md``.  This verifier checks:

* the exact Fourier formula for odd majority;
* orthogonality of the generated odd Walsh pole table;
* equality of all marginal deficit data in the coherent/diagonal pair;
* the exact selector losses and relative defects;
* realization by symmetric positive contractions; and
* the finite-group twirl and commutator bounds.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb

import numpy as np


def cube(p: int):
    return list(product((-1, 1), repeat=p))


def odd_masks(p: int):
    return [mask for mask in range(1 << p) if bin(mask).count("1") % 2 == 1]


def character(word, mask: int) -> int:
    out = 1
    for i, value in enumerate(word):
        if mask & (1 << i):
            out *= value
    return out


def majority_coefficients(p: int):
    assert p % 2 == 1
    words = cube(p)
    result = {}
    for mask in odd_masks(p):
        numerator = 0
        for word in words:
            tau = 1 if sum(word) > 0 else -1
            numerator += tau * character(word, mask)
        result[mask] = Fraction(numerator, 1 << p)
    return result


def predicted_coefficient(p: int, degree: int) -> Fraction:
    assert p % 2 == degree % 2 == 1
    m = (p - 1) // 2
    k = (degree - 1) // 2
    a1 = Fraction(comb(2 * m, m), 1 << (2 * m))
    return ((-1) ** k) * a1 * Fraction(comb(m, k), comb(2 * m, 2 * k))


def verify_exact_fourier():
    checks = 0
    for p in (3, 5, 7, 9):
        coeff = majority_coefficients(p)
        for mask, value in coeff.items():
            assert value == predicted_coefficient(p, bin(mask).count("1"))
            checks += 1
        assert sum(value * value for value in coeff.values()) == 1
        m = (p - 1) // 2
        a1 = Fraction(comb(2 * m, m), 1 << (2 * m))
        assert max(abs(value) for value in coeff.values()) == a1
    return checks


def pole_table(p: int):
    words = cube(p)
    masks = odd_masks(p)
    v = np.asarray(
        [[character(word, mask) for mask in masks] for word in words],
        dtype=float,
    )
    return masks, v


def verify_coherence_pair():
    checks = 0
    for p in (3, 5, 7, 9):
        masks, v = pole_table(p)
        q = len(masks)
        base_n = 1 << p
        g = v.T @ v / base_n
        assert np.allclose(g, np.eye(q), atol=1e-12)

        coeff_exact = majority_coefficients(p)
        a = np.asarray([float(coeff_exact[mask]) for mask in masks])
        singleton = np.asarray([bin(mask).count("1") == 1 for mask in masks])
        high = ~singleton
        rho = float(a[high] @ a[high])
        assert rho > 0
        u = np.zeros(q)
        u[high] = a[high] / np.sqrt(rho)

        d_coh = np.outer(u, u)
        d_diag = np.diag(u * u)
        r_coh = np.eye(q) - d_coh
        r_diag = np.eye(q) - d_diag

        assert np.allclose(np.diag(d_coh), np.diag(d_diag), atol=1e-12)
        assert np.allclose(d_coh[singleton, :], 0, atol=1e-12)
        assert np.allclose(d_diag[singleton, :], 0, atol=1e-12)
        assert np.isclose(np.trace(d_coh), 1)
        assert np.isclose(np.trace(d_diag), 1)
        assert np.isclose(np.linalg.eigvalsh(d_coh)[-1], 1)
        assert np.isclose(
            np.linalg.eigvalsh(d_diag)[-1], np.max(u * u)
        )

        m = (p - 1) // 2
        a1 = float(Fraction(comb(2 * m, m), 1 << (2 * m)))
        assert np.max(u * u) <= a1 * a1 / rho + 1e-12
        assert np.isclose(a @ d_coh @ a, rho, atol=1e-12)
        assert a @ d_diag @ a <= a1 * a1 + 1e-12

        # Repetition of every row type leaves normalized Gram data unchanged,
        # so the smaller cube table suffices to verify the N=2^(2p) theorem.
        e = v / np.sqrt(base_n)
        for rayleigh in (r_coh, r_diag):
            t = e @ rayleigh @ e.T
            eig = np.linalg.eigvalsh(t)
            assert eig[0] >= -1e-10
            assert eig[-1] <= 1 + 1e-10
            assert np.allclose(e.T @ t @ e, rayleigh, atol=1e-10)
        checks += q * q + 10
    return checks


def group_representations(p: int, masks):
    """Character-diagonal representation matrices for the Boolean group."""
    reps = []
    for word in cube(p):
        diagonal = [character(word, mask) for mask in masks]
        reps.append(np.diag(diagonal))
    return reps


def verify_twirling():
    rng = np.random.default_rng(20260817)
    checks = 0
    for p in (2, 3, 4):
        # Use all characters here, which are pairwise distinct.
        masks = list(range(1 << p))
        q = len(masks)
        reps = group_representations(p, masks)

        raw = rng.normal(size=(q, q))
        raw = (raw + raw.T) / 2
        # Scale and shift so A is self-adjoint and A <= I.
        raw /= max(1.0, np.linalg.norm(raw, 2))
        a = 0.7 * raw
        largest = np.linalg.eigvalsh(a)[-1]
        if largest > 0.9:
            a -= (largest - 0.9) * np.eye(q)

        twirl = sum(rep.T @ a @ rep for rep in reps) / len(reps)
        assert np.allclose(twirl, np.diag(np.diag(a)), atol=1e-12)

        eta = np.linalg.norm(a - twirl, 2)
        d = float(np.max(1 - np.diag(a)))
        defect = np.eye(q) - a
        assert np.linalg.eigvalsh(defect)[0] >= -1e-10
        assert np.linalg.eigvalsh(defect)[-1] <= d + eta + 1e-10

        average_comm = sum(
            np.linalg.norm(a @ rep - rep @ a, 2) for rep in reps
        ) / len(reps)
        assert eta <= average_comm + 1e-10

        generators = []
        for i in range(p):
            word = [1] * p
            word[i] = -1
            diagonal = [character(tuple(word), mask) for mask in masks]
            generators.append(np.diag(diagonal))
        gamma = max(np.linalg.norm(a @ rep - rep @ a, 2) for rep in generators)
        assert eta <= p * gamma / 2 + 1e-10
        checks += len(reps) + p + 4
    return checks


def main():
    total = 0
    total += verify_exact_fourier()
    total += verify_coherence_pair()
    total += verify_twirling()
    print(f"verified {total} Gram-relative coherence identities and bounds")


if __name__ == "__main__":
    main()
