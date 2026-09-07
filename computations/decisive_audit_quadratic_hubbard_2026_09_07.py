"""Finite-state diagnostics for cap-controlled quadratic fluctuations.

The proof is deterministic. Floating eigenvalues and Gaussian quadrature
here check its normalizations, including the quenched covariance order.
"""
from itertools import combinations, product
import math
import numpy as np
from scipy.special import ndtr


def cube(n):
    return np.asarray(list(product((-1.0, 1.0), repeat=n)))


def gaussian_abs(mean, variance):
    if variance == 0:
        return np.abs(mean)
    sigma = math.sqrt(variance)
    z = mean / sigma
    return sigma * math.sqrt(2 / math.pi) * np.exp(-z*z/2) + mean*(2*ndtr(z)-1)


def check(j, branch_field):
    n = len(j)
    x0 = cube(n)
    energies0 = np.einsum('ij,ij->i', x0, x0 @ j)/2
    cap = float(np.max(np.abs(energies0)))
    q = cap/n
    x = np.concatenate((x0, x0))
    s = np.concatenate((np.ones(len(x0)), -np.ones(len(x0))))
    energy = np.concatenate((energies0, energies0))
    fields = x @ j
    d = float(np.linalg.norm(j, 2))
    y_given_x = gaussian_abs(d*x+s[:, None]*fields, d)
    edges = list(combinations(range(n), 2))
    statistics = np.stack([x[:, i]*x[:, k] for i, k in edges], axis=1)
    signed_statistics = statistics*s[:, None]
    if branch_field:
        nodes, weights = np.polynomial.hermite.hermgauss(48)
        nodes *= math.sqrt(2)
        weights /= math.sqrt(math.pi)
    else:
        nodes, weights = np.asarray([0.0]), np.asarray([1.0])
    y = np.zeros(n)
    field_norm = 0.0
    signed_energy = 0.0
    covariances = [np.zeros((len(edges), len(edges))) for _ in range(2)]
    for g, weight in zip(nodes, weights):
        logits = s*(energy+branch_field*g)
        p = np.exp(logits-logits.max())
        p /= p.sum()
        y += weight*(p @ y_given_x)
        field_norm += weight*float(p @ np.abs(fields).sum(axis=1))
        signed_energy += weight*float(p @ (s*energy))
        for covariance, statistic in zip(covariances, (statistics, signed_statistics)):
            mean = p @ statistic
            covariance += weight*(statistic.T @ (p[:, None]*statistic)-np.outer(mean, mean))
    assert field_norm <= min(4*cap, 2*cap+n)+2e-10
    assert field_norm <= 2*signed_energy+n+2e-10
    budget = min(4*q, 2*q+1)+d+math.sqrt(2*d/math.pi)
    assert y.sum() <= n*budget+2e-10
    diagonal = np.asarray([math.exp(-2*(y[i]+y[k])) for i, k in edges])
    for covariance in covariances:
        minimum = np.linalg.eigvalsh(covariance-np.diag(diagonal)).min()
        assert minimum > -2e-10, (n, branch_field, minimum)
        sign, logdet = np.linalg.slogdet(covariance)
        assert sign > 0
        assert logdet >= -2*(n-1)*y.sum()-2e-8
    # A whole flat bridge has eta=N/min(m,N-m), not twice that number.
    for m in range(1, n):
        coeff = np.asarray([j[i, k] != 0 if (i < m <= k) else 0 for i, k in edges], dtype=float)
        # For the sign matrices below every cross entry is nonzero.
        v = float(coeff @ coeff)
        eta = n/min(m, n-m)
        for covariance in covariances:
            assert coeff @ covariance @ coeff >= v*math.exp(-2*eta*budget)-2e-10
    return 2, 2*(n-1)


def main():
    checks = bridges = systems = 0
    for n in (3, 4, 5):
        free = list(combinations(range(1, n), 2))
        for signs in product((-1.0, 1.0), repeat=len(free)):
            a = np.ones((n, n))-np.eye(n)
            for (i, k), sign in zip(free, signs):
                a[i, k] = a[k, i] = sign
            for beta in (0.2, 0.9, 2.0):
                for h in (0.0, 1.3):
                    c, b = check(beta*a/math.sqrt(n), h)
                    checks += c
                    bridges += b
                    systems += 1
    rng = np.random.RandomState(20260907)
    for n in (6, 7, 8):
        for _ in range(8):
            upper = np.triu(rng.choice((-1.0, 1.0), size=(n, n)), 1)
            a = upper+upper.T
            for beta in (0.4, 1.5):
                c, b = check(beta*a/math.sqrt(n), 0.7)
                checks += c
                bridges += b
                systems += 1
    print('PASS: %d Gibbs systems, %d full quadratic covariance/determinant tests, %d bridges'
          % (systems, checks, bridges))


if __name__ == '__main__':
    main()
