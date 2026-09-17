"""Finite diagnostics for the proved critical-scale realization theorem.

No randomized diagnostic replaces the uniform proof. Outputs JSON to stdout.
Run from repository root with .venv/bin/python.
"""

import itertools
import json
import math

import numpy as np
from scipy.linalg import expm, hadamard


SEED = 2026091709
rng = np.random.default_rng(SEED)
C0 = 16 * math.sqrt(math.e * math.pi / 2)


def signs(d):
    return np.array(list(itertools.product((-1.0, 1.0), repeat=d)))


def law(k):
    if k & (k - 1):
        atoms = signs(k)
        return atoms, np.full(len(atoms), 1 / len(atoms))
    h = hadamard(k).astype(float)
    h2 = h[:, rng.permutation(k)] * rng.choice((-1.0, 1.0), size=k)
    atoms = np.concatenate((h, h2, -h))
    weights = np.concatenate((np.full(k, .2 / k),
                              np.full(k, .5 / k),
                              np.full(k, .3 / k)))
    return atoms, weights


def nuclear(t):
    return np.linalg.svd(t, compute_uv=False).sum()


def response(t, atoms, weights):
    return weights @ np.linalg.norm(atoms @ t.T, axis=1)


cases = []
mgf_error = 0.0
nuclear_tests = 0
gaussian_increment_tests = 0
absolute_mgf_tests = 0
max_relative_slack_ratio = 0.0
for k in (2, 3, 4, 8):
    atoms, weights = law(k)
    cov = atoms.T @ (weights[:, None] * atoms)
    assert np.max(np.abs(cov - np.eye(k))) < 1e-12
    for theta in (-.2, .03, .2):
        actual = sum(w * expm(theta * np.outer(h, h))
                     for h, w in zip(atoms, weights))
        wanted = (1 + np.expm1(theta * k) / k) * np.eye(k)
        mgf_error = max(mgf_error, np.max(np.abs(actual - wanted)))
        assert np.max(np.abs(actual - wanted)) < 2e-12
    q = math.ceil(120 * k * math.log(4 * k))
    delta = math.sqrt(3 * k * math.log(4 * k) / q)
    labels = atoms[rng.choice(len(atoms), q, p=weights)]
    gram = labels.T @ labels
    eig = np.linalg.eigvalsh(gram / q)
    assert eig[-1] <= 1 + delta
    max_nuclear_ratio = 0.0
    for _ in range(150):
        p = int(rng.integers(1, k + 3))
        t = rng.normal(size=(p, k))
        if rng.random() < .5:
            t = np.outer(rng.normal(size=p), rng.normal(size=k))
        phi = response(t, atoms, weights)
        nu = nuclear(t)
        max_nuclear_ratio = max(max_nuclear_ratio, nu / (math.sqrt(k)*phi))
        assert nu <= math.sqrt(k) * phi + 1e-10
        empirical = np.linalg.norm(labels @ t.T, axis=1).mean()
        rhs = phi + C0 * math.sqrt(k / q) * nu
        assert empirical <= rhs + 1e-10
        max_relative_slack_ratio = max(max_relative_slack_ratio,
            max(0, empirical - phi) / (C0 * math.sqrt(k/q) * nu))
        nuclear_tests += 1
        u = rng.normal(size=(p, k))
        inc_norm = np.linalg.norm(labels @ t.T, axis=1) - np.linalg.norm(labels @ u.T, axis=1)
        lhs = inc_norm @ inc_norm
        rhs = np.linalg.norm(labels @ (t-u).T)**2
        assert lhs <= rhs + 1e-8
        gaussian_increment_tests += 1
    # Exact single-column laws, multiplied across independent columns.
    p, ell = 2, 1
    n = k*p + ell
    scalar_atoms = signs(p + ell)
    for _ in range(8):
        w = rng.normal(size=n)
        core, left = w[:k*p].reshape(k,p), w[k*p:]
        coeff = np.concatenate((labels @ core,
                                np.tile(left, (q,1))), axis=1)
        variance = np.sum(coeff**2)
        assert variance <= q*(1+delta)*(w@w) + 1e-8
        column_values = np.abs(coeff @ scalar_atoms.T)
        means = column_values.mean(axis=1)
        for lam in (-.3, -.05, .05, .3):
            log_mgf = np.log(np.exp(lam*(column_values-means[:,None])).mean(axis=1)).sum()
            assert log_mgf <= lam**2 * variance / 2 + 1e-9
            absolute_mgf_tests += 1
    cases.append(dict(k=k, q=q, delta=delta,
                      gram_eigenvalues=eig.tolist(),
                      max_nuclear_ratio=max_nuclear_ratio))

# Exact fixed-cardinality flip-energy identity and counting obstruction.
flip_tests = 0
for n in range(4, 10):
    words = signs(n)
    a = rng.choice((-1.0, 1.0), size=(n,n))
    a = np.triu(a, 1)
    a = a + a.T
    energies = np.einsum('bi,ij,bj->b', words, a, words)/2
    idx = np.argmax(np.abs(energies))
    x = words[idx]
    sigma = np.sign(energies[idx])
    cap = abs(energies[idx])
    distances = np.sum(words != x, axis=1)
    deficits = cap - sigma*energies
    assert deficits.min() >= -1e-10
    for r in range(n+1):
        avg = deficits[distances == r].mean()
        expected = 4*r*(n-r)*cap/(n*(n-1))
        assert abs(avg-expected) < 1e-10
        flip_tests += 1

print(json.dumps(dict(status='PASS', seed=SEED, C0=C0,
    cases=cases, scalar_matrix_mgf_max_error=float(mgf_error),
    nuclear_tests=nuclear_tests, gaussian_increment_tests=gaussian_increment_tests,
    centered_absolute_mgf_tests=absolute_mgf_tests,
    max_sampled_empirical_error_fraction_of_bound=max_relative_slack_ratio,
    exact_flip_identity_tests=flip_tests,
    scope='Finite diagnostics, not proof of uniform claims or optimizer geometry'), indent=2))
