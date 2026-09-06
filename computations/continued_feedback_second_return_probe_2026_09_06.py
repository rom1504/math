#!/usr/bin/env python3
"""Focused falsifier of a fresh-noise law for the SECOND feedback return.

This tests η=B(C-c0-aZ) against L=B(c0+aZ), never replacing the actual
coherent return. Regression/pilot and quadrature uncertainties are NOT
included in the conditional Monte Carlo standard errors.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
from numpy.polynomial.hermite import hermgauss

from continued_director_steiner_signing_2026_09_06 import make_signing
from continued_feedback_coherent_energy_projection_2026_09_06 import signing
from continued_feedback_hidden_third_return_2026_09_06 import construct
from continued_feedback_variance_falsifier_2026_09_06 import SEED


def run(kind, size, samples, batch, seed, angle, response_scale):
    rng = np.random.default_rng(seed)
    if kind == "hidden":
        a = construct(size, seed).astype(float)
    elif kind == "steiner_seed":
        a = np.kron(SEED, make_signing(size)[0]).astype(float)
        a[a == 0] = 1
        np.fill_diagonal(a, 0)
    else:
        a = signing(kind, size, rng)
    n = len(a)
    b = a / np.sqrt(n - 1)
    q = b @ b
    first = np.exp(-0.5)
    sigma = np.sqrt(np.maximum(np.diag(b @ (np.exp(-1) * (np.sinh(q) - q)) @ b), 0))
    nodes, weights = hermgauss(48)
    nodes *= np.sqrt(2)
    weights /= np.sqrt(np.pi)

    def fields(count):
        s = rng.choice((-1.0, 1.0), size=(count, n))
        g = s @ b
        d = s * (g**2 - 1) / np.sqrt(2)
        old = d @ b
        score = np.cos(angle) * g + np.sin(angle) * old
        coherent = first * (np.cos(angle) * (s @ q) + np.sin(angle) * (old @ b))
        f = np.sin(score)
        h = 1 - np.abs(f)
        feedback = f @ b
        z = feedback - coherent
        c = h * np.tanh(response_scale * feedback)
        c0 = np.zeros_like(c)
        derivative = np.zeros_like(c)
        for node, weight in zip(nodes, weights):
            response = np.tanh(response_scale * (coherent + node * sigma))
            c0 += weight * h * response
            derivative += weight * h * response_scale * (1 - response**2)
        return f, feedback, z, c, c0, derivative

    coefficient = np.zeros(n)
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        coefficient += fields(count)[5].sum(axis=0)
        done += count
    coefficient /= samples

    regression_moments = np.zeros((3, n))
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        _, _, z, c, c0, _ = fields(count)
        linear_return = (c0 + z * coefficient) @ b
        innovation = c @ b - linear_return
        regression_moments[0] += (innovation**2).sum(axis=0)
        regression_moments[1] += (linear_return**2).sum(axis=0)
        regression_moments[2] += (innovation * linear_return).sum(axis=0)
        done += count
    regression_moments /= samples
    variance_eta, variance_l, covariance = regression_moments
    kept = variance_eta > max(float(variance_eta.mean()) * 1e-3, 1e-12)
    beta = np.divide(covariance, variance_eta, out=np.zeros(n), where=kept)
    variance_r = variance_l - 2 * beta * covariance + beta**2 * variance_eta
    normalizer = float(np.mean(variance_eta[kept] * variance_r[kept]))
    fourth_normalizer = float(np.mean(variance_eta[kept] ** 2))

    accum = np.zeros((7, 2))
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        f, feedback, z, c, c0, _ = fields(count)
        bc = c @ b
        linear_return = (c0 + z * coefficient) @ b
        eta = bc - linear_return
        residual = linear_return - beta * eta
        eta_keep = eta[:, kept]
        r_keep = residual[:, kept]
        square_covariance = np.mean(
            (eta_keep**2 - variance_eta[kept]) * (r_keep**2 - variance_r[kept]), axis=1
        )
        fourth_cumulant = np.mean(
            eta_keep**4 - 6 * variance_eta[kept] * eta_keep**2 + 3 * variance_eta[kept] ** 2,
            axis=1,
        )
        endpoint = f + c
        next_field = feedback + bc
        next_gain = np.mean((1 - np.abs(endpoint)) * np.abs(next_field), axis=1)
        data = np.stack((
            square_covariance,
            np.mean(eta_keep * r_keep**3, axis=1),
            np.mean(eta_keep**3 * r_keep, axis=1),
            fourth_cumulant,
            np.mean(eta_keep * r_keep, axis=1),
            next_gain,
            np.mean(eta_keep**2, axis=1),
        ))
        accum[:, 0] += data.sum(axis=1)
        accum[:, 1] += (data**2).sum(axis=1)
        done += count
    means = accum[:, 0] / samples
    errors = np.sqrt(np.maximum(accum[:, 1] / samples - means**2, 0) / samples)
    return {
        "status": "Falsifier only; errors conditional on BOTH pilots and 48-node quadrature",
        "family": kind, "size_parameter": size, "n": n, "seed": seed,
        "samples_each_of_three_stages": samples, "angle": angle, "response_scale": response_scale,
        "kept_fraction": float(kept.mean()), "mean_eta_variance": float(means[6]),
        "square_covariance": float(means[0]), "square_covariance_se": float(errors[0]),
        "normalized_square_covariance": float(means[0] / normalizer),
        "normalized_square_covariance_se": float(errors[0] / normalizer),
        "eta_residual_cubic": float(means[1]), "eta_residual_cubic_se": float(errors[1]),
        "eta_cubic_residual": float(means[2]), "eta_cubic_residual_se": float(errors[2]),
        "normalized_eta_fourth_cumulant": float(means[3] / fourth_normalizer),
        "normalized_eta_fourth_cumulant_se": float(errors[3] / fourth_normalizer),
        "residual_linear_covariance": float(means[4]),
        "next_actual_endpoint_gain": float(means[5]), "next_gain_se": float(errors[5]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("steiner", "steiner_seed", "hidden", "random", "seed_lift"), default="steiner")
    parser.add_argument("--sizes", type=int, nargs="+", default=[4, 5])
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--batch", type=int, default=200)
    parser.add_argument("--seed", type=int, default=10092026)
    parser.add_argument("--angle", type=float, default=float(np.pi / 4))
    parser.add_argument("--response-scale", type=float, default=2)
    args = parser.parse_args()
    for size in args.sizes:
        print(json.dumps(run(args.kind, size, args.samples, args.batch, args.seed + size,
                             args.angle, args.response_scale)), flush=True)


if __name__ == "__main__":
    main()
