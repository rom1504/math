"""Exact conference/sign identity checks and finite covariance diagnostics."""
import argparse
import json
from pathlib import Path

import numpy as np


def conference(q):
    assert q % 4 == 1
    assert all(q % d for d in range(2, int(q**0.5)+1))
    chi = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        chi[x] = 1 if pow(x, (q-1)//2, q) == 1 else -1
    matrix = np.ones((q+1, q+1), dtype=np.int64)
    matrix[0, 0] = 0
    for a in range(q):
        for b in range(q):
            matrix[a+1, b+1] = chi[(a-b) % q]
    assert np.array_equal(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0)
    assert np.array_equal(matrix@matrix, q*np.eye(q+1, dtype=np.int64))
    return matrix


def run(q, samples):
    bmat = conference(q)/np.sqrt(q)
    n = q+1
    kappa = h = 0.25
    first = kappa*np.exp(-0.5)
    t = np.pi/(2*first)
    tau2 = kappa*kappa*np.exp(-1)*(np.sinh(1)-1)
    variance = t*t*tau2
    mean_cos = np.exp(-variance/2)
    exhaustive = n <= 18
    total = 2**n if exhaustive else samples
    rng = np.random.default_rng(9062026+q)
    sums = np.zeros((2, 2))
    identity_error = 0.0
    for start in range(0, total, 1024):
        count = min(1024, total-start)
        if exhaustive:
            index = np.arange(start, start+count, dtype=np.uint64)
            seeds = 2*((index[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(float)-1
        else:
            seeds = rng.choice((-1.0, 1.0), size=(count, n))
        g = seeds@bmat
        f = kappa*np.sin(g)
        noise = (f-first*g)@bmat
        actual_response = h*np.sin(t*(f@bmat))
        parity_response = h*seeds*np.cos(t*noise)
        identity_error = max(identity_error, float(np.abs(actual_response-parity_response).max()))
        c0 = h*mean_cos*seeds
        residual = actual_response-c0
        values = np.stack((np.mean(residual**2, axis=1),
                           np.mean(residual*c0, axis=1)))
        sums[:, 0] += values.sum(axis=1)
        sums[:, 1] += (values**2).sum(axis=1)
    means = sums[:, 0]/total
    se = np.zeros(2) if exhaustive else np.sqrt(np.maximum(sums[:, 1]/total-means**2, 0)/total)
    assert identity_error < 1e-12
    return {'q': q, 'n': n, 'integer_conference_identity': True,
        'exhaustive_boolean_average': exhaustive, 'seed_count': total,
        'maximum_trigonometric_identity_error': identity_error,
        'actual_average_eta_variance': float(means[0]),
        'eta_variance_sampling_se': float(se[0]),
        'actual_signed_average_rho': float(means[1]),
        'rho_sampling_se': float(se[1]),
        'proved_limit_even_noise_variance': float(h*h*(1-np.exp(-variance))**2/2),
        'incorrect_zero_first_sine_prediction': float(h*h*np.exp(-variance)*(np.sinh(variance)-variance)),
        'proved_positive_limiting_difference': float(h*h*np.exp(-variance)*(variance-1+np.exp(-variance))),
        'scope': 'Integer signing identity exact; transcendental and finite averages floating; asymptotic conclusion proved separately.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--primes', type=int, nargs='+', default=[5, 13, 17, 29])
    parser.add_argument('--samples', type=int, default=16000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records = [run(q, args.samples) for q in args.primes]
    args.output.write_text(json.dumps(records, indent=2)+'\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
