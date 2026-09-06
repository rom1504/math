"""Exact finite-product identities and MC diagnostics for actual Boolean feedback.

No Gaussian substitution is made for Q S. Small complete seed enumerations
check the characteristic-function algebra, not the asymptotic energy theorem.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def inserted_products(angles):
    cosine = np.cos(angles)
    left = np.concatenate((np.ones((len(angles), 1)), np.cumprod(cosine[:, :-1], axis=1)), axis=1)
    right = np.concatenate((np.cumprod(cosine[:, :0:-1], axis=1)[:, ::-1], np.ones((len(angles), 1))), axis=1)
    return np.prod(cosine, axis=1), np.sin(angles) * left * right


def certificate(a, kappa=0.6, frequency=1.0):
    n = len(a)
    bmat = a.astype(float) / np.sqrt(n - 1)
    q = bmat @ bmat
    b3 = q @ bmat
    first = kappa * np.exp(-0.5)
    height = 1 - kappa
    residual_cov = first**2 * (np.sinh(q) - q)
    transported = bmat @ residual_cov @ bmat
    sigma2 = np.diag(transported)
    damping = height * np.exp(-0.5 * frequency**2 * sigma2)
    angles = frequency * first * q
    products, inserted = inserted_products(angles)
    projection = damping[:, None] * inserted
    derivative = frequency * damping * products
    kernel = np.empty((n, n))
    for i in range(n):
        kernel[i] = kappa * damping[i] / 2 * (
            np.prod(np.cos(angles[i][None, :] - bmat), axis=1)
            - np.prod(np.cos(angles[i][None, :] + bmat), axis=1)
        )
    kernel -= first * (projection @ bmat)
    # Tr(B M M^T) is sum((B M) * M), not Tr(M B M^T).
    coherent_energy = float(np.sum((bmat @ projection) * projection) / (2 * n))
    cross_energy = float(np.sum(((bmat * derivative[None, :]) @ bmat) * kernel) / n)
    noise_energy = float(np.sum(bmat * (derivative[:, None] * transported * derivative[None, :])) / (2 * n))
    gain = float(np.mean(damping * (first * np.sum(q * inserted, axis=1) + frequency * sigma2 * products)))
    old_energy = first**2 * float(np.trace(b3)) / (2 * n)
    return {"B": bmat, "Q": q, "sigma2": sigma2, "damping": damping,
            "M": projection, "K": kernel, "a": derivative, "first": first,
            "terms": {"coherent": coherent_energy, "bare_star": cross_energy,
                      "linear_noise": noise_energy, "gain": gain, "old": old_energy,
                      "feedback": coherent_energy + cross_energy + noise_energy,
                      "endpoint": abs(gain) + abs(old_energy + coherent_energy + cross_energy + noise_energy)}}


def enumerate_identities(a, kappa, frequency):
    data = certificate(a, kappa, frequency)
    n = len(a)
    seed = 1 - 2 * ((np.arange(2**n)[:, None] >> np.arange(n)) & 1)
    g = seed @ data["B"]
    v = seed @ data["Q"]
    c0 = data["damping"][None, :] * np.sin(frequency * data["first"] * v)
    m = c0.T @ seed / len(seed)
    residual = kappa * np.sin(g) - data["first"] * g
    kernel = c0.T @ residual / len(seed)
    derivative = frequency * data["damping"] * np.mean(np.cos(frequency * data["first"] * v), axis=0)
    errors = {"M": float(np.max(np.abs(m - data["M"]))),
              "K": float(np.max(np.abs(kernel - data["K"]))),
              "a": float(np.max(np.abs(derivative - data["a"])))}
    assert max(errors.values()) < 1e-11
    return {"n": n, "seed_count": len(seed), "max_errors": errors,
            "status": "complete finite seed enumeration, floating identity check"}


def diagnose(a, samples, batch, seed_value, kappa, frequency):
    data = certificate(a, kappa, frequency)
    n = len(a)
    rng = np.random.default_rng(seed_value)
    energy, gain = [], []
    for start in range(0, samples, batch):
        count = min(batch, samples - start)
        seed = rng.choice([-1.0, 1.0], size=(count, n))
        old = kappa * np.sin(seed @ data["B"])
        returned = old @ data["B"]
        response = (1 - kappa) * np.sin(frequency * returned)
        energy.extend(np.sum(response * (response @ data["B"]), axis=1) / (2 * n))
        gain.extend(np.sum(response * returned, axis=1) / n)

    def report(values, target):
        values = np.asarray(values)
        return {"mean": float(np.mean(values)), "standard_error": float(np.std(values, ddof=1) / np.sqrt(samples)),
                "projected": target, "finite_difference": float(np.mean(values)) - target}

    return {"n": n, "samples": samples, "seed": seed_value, "terms": data["terms"],
            "feedback": report(energy, data["terms"]["feedback"]), "gain": report(gain, data["terms"]["gain"]),
            "status": "Monte Carlo diagnostics, not an asymptotic certificate"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steiner-m", nargs="+", type=int, default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=12000)
    parser.add_argument("--batch", type=int, default=300)
    parser.add_argument("--seed", type=int, default=9060851)
    parser.add_argument("--kappa", type=float, default=0.6)
    parser.add_argument("--frequency", type=float, default=1.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    small, _ = make_signing(2)
    rng = np.random.default_rng(args.seed)
    upper = np.triu(rng.choice([-1, 1], size=(10, 10)), 1)
    checks = [enumerate_identities(a, args.kappa, args.frequency) for a in [small, upper + upper.T]]
    print(json.dumps({"enumeration": checks}), flush=True)
    records = []
    for m in args.steiner_m:
        a, _ = make_signing(m)
        record = diagnose(a, args.samples, args.batch, args.seed + m, args.kappa, args.frequency)
        records.append(record)
        print(json.dumps(record), flush=True)
    if args.output:
        args.output.write_text(json.dumps({"kappa": args.kappa, "frequency": args.frequency,
                                         "enumeration": checks, "records": records}, indent=2) + "\n")


if __name__ == "__main__":
    main()
