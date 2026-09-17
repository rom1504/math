"""Finite replay of Paley pinning, exact sign bridges, and certificate floor."""

from __future__ import annotations

import itertools
import json
import math

import numpy as np


def entropy(value):
    if value <= 0 or value >= 1:
        return 0.0
    return -value*math.log(value)-(1-value)*math.log1p(-value)


def paley_pinned(prime, target, rng):
    chi = np.zeros(prime, dtype=np.int64)
    for value in range(1, prime):
        chi[value] = 1 if pow(value, (prime-1)//2, prime) == 1 else -1
    assert chi[-1] == 1 and chi.sum() == 0
    base_eigenvalues = np.fft.fft(chi)
    assert np.max(np.abs(base_eigenvalues.imag)) < 1e-10
    assert np.max(np.abs(np.abs(base_eigenvalues[1:])-math.sqrt(prime))) < 1e-10
    pairs = [value for value in range(1, (prime+1)//2) if chi[value] == -1]
    count = int(math.floor(target*math.sqrt(prime)/4))
    assert 0 < count < len(pairs)
    for attempt in range(1000):
        selected = rng.choice(pairs, count, replace=False)
        sequence = chi.copy()
        for value in selected:
            sequence[value] = sequence[-value] = 1
        eigenvalues = np.fft.fft(sequence).real
        spectral_remainder = float(np.max(np.abs(eigenvalues[1:])))
        degree = 4*count
        if degree > spectral_remainder:
            break
    else:
        raise AssertionError("No numerically pinned sample found.")
    matrix = sequence[(np.arange(prime)[:, None]-np.arange(prime)[None, :]) % prime]
    assert np.array_equal(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0)
    assert set(matrix[np.triu_indices(prime, 1)]) <= {-1, 1}
    assert np.all(matrix.sum(axis=1) == degree)
    perturbation = eigenvalues-base_eigenvalues.real
    for frequency in range(1, prime):
        population = np.asarray([4*math.cos(2*math.pi*frequency*value/prime)
                                 for value in pairs])
        predicted_mean = -4*count*(1+base_eigenvalues[frequency].real)/(prime-1)
        assert abs(count*population.mean()-predicted_mean) < 1e-10
        assert abs(perturbation[frequency]-predicted_mean) <= math.sqrt(
            512*count*math.log(prime))+1e-10
    return matrix, degree, spectral_remainder, selected.tolist()


def exhaustive_checks(matrix, degree, remainder, rng):
    n = len(matrix)
    indices = np.arange(2**(n-1), dtype=np.int64)
    words = np.ones((len(indices), n), dtype=np.int64)
    words[:, 1:] = 1-2*((indices[:, None] >> np.arange(n-1)) & 1)
    energies = np.sum((words @ matrix)*words, axis=1)//2
    cap = n*degree//2
    assert int(np.max(np.abs(energies))) == cap
    radii = np.minimum(np.sum(words < 0, axis=1), np.sum(words > 0, axis=1))
    rho = radii/n
    assert np.all(cap-energies >= 2*(degree-remainder)*n*rho*(1-rho)-1e-9)
    assert np.max(-energies) <= remainder*n/2+1e-9
    assert np.sum(np.abs(energies) == cap) == 1  # One representative of ±1.

    q = 3
    bridge = np.empty((n, q), dtype=np.int64)
    drivers = rng.choice((-1, 1), size=(n//2, q))
    bridge[0:2*(n//2):2] = drivers
    bridge[1:2*(n//2):2] = -drivers
    if n % 2:
        bridge[-1] = rng.choice((-1, 1), size=q)
    assert np.all(np.abs(bridge.sum(axis=0)) == n % 2)
    responses = words @ bridge
    absolute_response = np.sum(np.abs(responses), axis=1)
    for radius in range(n//2+1):
        mask = radii == radius
        if not np.any(mask):
            continue
        u = math.log(2*math.comb(n, radius))+3*math.log(n+1)
        bound = q+2*q*math.sqrt(radius)+math.sqrt(8*q*radius*u)
        assert np.max(absolute_response[mask]) <= bound+1e-9

    new_words = np.asarray(list(itertools.product((-1, 1), repeat=q)), dtype=np.int64)
    child = np.eye(q, dtype=np.int64)-np.ones((q, q), dtype=np.int64)
    child_energies = np.sum((new_words @ child)*new_words, axis=1)//2
    parent_values = energies[:, None]+child_energies[None, :]+responses @ new_words.T
    parent_cap = int(np.max(np.abs(parent_values)))
    direct_triangle = int(np.max(np.abs(energies)+absolute_response)+np.max(np.abs(child_energies)))
    assert parent_cap <= direct_triangle
    return {"spin_pairs_checked": len(words), "exact_old_cap": cap,
            "exact_parent_cap_q3": parent_cap,
            "old_normalized_cap": cap/n**1.5,
            "parent_normalized_cap_q3": parent_cap/(n+q)**1.5}


def main():
    rng = np.random.default_rng(2026091723)
    reports = []
    for prime in (13, 17, 29, 101, 257):
        matrix, degree, remainder, selected = paley_pinned(prime, 3.0, rng)
        report = {"order": prime, "row_sum": degree, "orthogonal_norm_numeric": remainder,
                  "selected_negative_pairs": selected}
        if prime <= 17:
            report.update(exhaustive_checks(matrix, degree, remainder, rng))
        reports.append(report)
    for child_constant in (0.0, 0.4333221116640807, 0.493608094, 0.5):
        for epsilon in np.logspace(-8, 8, 1001):
            lower = (0.5+math.sqrt(2)*epsilon+2*math.sqrt(epsilon*math.log(2))
                     +child_constant*epsilon**1.5)/(1+epsilon)**1.5
            assert lower > child_constant
    print(json.dumps({"status": "all finite Paley, bridge, and floor checks passed",
                      "spectral_note": "Fourier norms numerical; asymptotic proof is analytic.",
                      "cases": reports, "certificate_grid_checks": 4004}, indent=2))


if __name__ == "__main__":
    main()

