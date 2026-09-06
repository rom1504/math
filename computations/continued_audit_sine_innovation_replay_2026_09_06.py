"""Actual-signing finite diagnostics for the proved sine innovation floor.

The interval constant is certified. Matrix inequalities are floating
diagnostics of exact formulas; Monte Carlo values are not asymptotic proofs.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpmath import iv

from continued_director_steiner_signing_2026_09_06 import make_signing


def constant_interval():
    iv.dps = 70
    one = iv.mpf(1)
    sinh_one = (iv.exp(one)-iv.exp(-one))/2
    value = iv.sqrt(2 / iv.pi) * iv.exp(-one) / (48 * iv.sqrt(sinh_one - one))
    assert value.a > 0
    return str(value)


def run(m, samples, batch, seed):
    signing, metadata = make_signing(m)
    n = len(signing)
    b = signing.astype(float) / np.sqrt(n - 1)
    q = b @ b
    cap = 2.0
    actual_op = max(abs(metadata['positive_eigenvalue']),
                    abs(metadata['negative_eigenvalue'])) / np.sqrt(n - 1)
    assert actual_op <= cap
    alpha = np.exp(1.5) / 2
    kappa = 1 / (4 * (1 + alpha))

    def source(x):
        return kappa * (np.sin(x) - alpha * np.sin(2 * x))

    def gaussian_source_cov(correlation):
        return kappa**2 * (np.exp(-1) * np.sinh(correlation)
            - 2 * alpha * np.exp(-2.5) * np.sinh(2 * correlation)
            + alpha**2 * np.exp(-4) * np.sinh(4 * correlation))

    tau2 = float(gaussian_source_cov(1.0))
    tau = np.sqrt(tau2)
    t = 1 / (cap * tau)
    h = 0.25
    rf = gaussian_source_cov(q)
    covariance = b @ rf @ b
    diagonal = np.diag(covariance)
    weights = np.exp(-t*t*diagonal/2)
    coefficient = h*t*weights
    residual_cov = h*h*weights[:, None]*weights[None, :] * (
        np.sinh(t*t*covariance) - t*t*covariance)
    ideal_eta_cov = b @ residual_cov @ b
    ideal_v = np.diag(ideal_eta_cov)
    a_floor = h*h*t**6*tau**6*np.exp(-1)/(6*cap**2)
    row_cap = h*h*cap**2*(np.sinh(1)-1)
    psd_margin = float(np.linalg.eigvalsh(cap**2*tau2*q-covariance).min())
    assert psd_margin >= -1e-10
    assert diagonal.max() <= cap**2*tau2+1e-10
    assert diagonal.mean() >= tau2-1e-10
    assert ideal_v.mean() >= a_floor-1e-10
    assert ideal_v.max() <= row_cap+1e-10

    rng = np.random.default_rng(seed+m)
    moments = np.zeros((3, 2))
    row_v = np.zeros(n)
    row_rho = np.zeros(n)
    row_rho2 = np.zeros(n)
    done = 0
    while done < samples:
        count = min(batch, samples-done)
        spins = rng.choice((-1.0, 1.0), size=(count, n))
        g = spins @ b
        f = source(g)
        z = f @ b
        c = h*np.sin(t*z)
        coherent = (z*coefficient) @ b
        eta = c @ b-coherent
        endpoint = f+c
        gain = np.mean((1-np.abs(endpoint))*np.abs(z+c@b), axis=1)
        rho_samples = eta*coherent
        values = np.stack((np.mean(eta**2, axis=1),
                           np.mean(rho_samples, axis=1), gain))
        moments[:, 0] += values.sum(axis=1)
        moments[:, 1] += (values**2).sum(axis=1)
        row_v += (eta**2).sum(axis=0)
        row_rho += rho_samples.sum(axis=0)
        row_rho2 += (rho_samples**2).sum(axis=0)
        done += count
    means = moments[:, 0]/samples
    errors = np.sqrt(np.maximum(moments[:, 1]/samples-means**2, 0)/samples)
    row_v /= samples
    row_rho /= samples
    row_se = np.sqrt(np.maximum(row_rho2/samples-row_rho**2, 0)/samples)
    return {
        'n': n, 'steiner_parameter': m, 'actual_operator_norm': actual_op,
        'fixed_operator_cap': cap, 'tau_squared': tau2, 'response_t': t,
        'exact_formula_float_checks': {
            'min_eigenvalue_L2tau2Q_minus_T': psd_margin,
            'mean_old_variance': float(diagonal.mean()),
            'min_old_variance': float(diagonal.min()),
            'mean_ideal_eta_variance': float(ideal_v.mean()),
            'max_ideal_eta_variance': float(ideal_v.max()),
            'proved_variance_floor': a_floor, 'proved_ideal_row_cap': row_cap},
        'monte_carlo': {
            'samples': samples, 'seed': seed+m,
            'mean_eta_variance': float(means[0]), 'eta_variance_se': float(errors[0]),
            'signed_mean_rho': float(means[1]), 'signed_mean_rho_se': float(errors[1]),
            'mean_absolute_estimated_row_rho': float(np.abs(row_rho).mean()),
            'mean_row_rho_standard_error': float(row_se.mean()),
            'mean_absolute_estimated_variance_error': float(np.abs(row_v-ideal_v).mean()),
            'actual_next_gain': float(means[2]), 'next_gain_se': float(errors[2])},
        'scope': 'Finite diagnostics, not a convergence-rate or theorem certification.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--parameters', type=int, nargs='+', default=[3, 4, 5])
    parser.add_argument('--samples', type=int, default=4000)
    parser.add_argument('--batch', type=int, default=100)
    parser.add_argument('--seed', type=int, default=6092026)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records = []
    for m in args.parameters:
        result = run(m, args.samples, args.batch, args.seed)
        records.append(result)
        print(json.dumps(result), flush=True)
    result = {'certified_gain_coefficient_interval': constant_interval(),
              'records': records}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['certified_gain_coefficient_interval'], flush=True)


if __name__ == '__main__':
    main()
