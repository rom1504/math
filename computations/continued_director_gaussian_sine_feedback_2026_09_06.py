"""Reproducible finite diagnostics for an explicit Gaussian-seed energy law.

The deterministic projected value is evaluated without numerical integration.
Monte Carlo values and standard errors are diagnostics, not proof certificates.
Use repository-local tmp for any output not intended as a committed artifact.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def evaluate(a, samples, batch, seed, kappa, frequency):
    n = len(a)
    bmat = a.astype(float) / np.sqrt(n - 1)
    q = bmat @ bmat
    b3 = q @ bmat
    first = kappa * np.exp(-0.5)
    height = 1 - kappa
    total_cov = first**2 * (bmat @ np.sinh(q) @ bmat)
    variances = np.diag(total_cov)
    derivative = height * frequency * np.exp(-0.5 * frequency**2 * variances)
    weighted_b = bmat * derivative[None, :]
    linear_energy = float(np.sum(bmat * (derivative[:, None] * total_cov * derivative[None, :])) / (2 * n))
    cross_kernel = (first / frequency) * derivative[:, None] * (
        np.sinh(frequency * first * b3) - frequency * first * b3
    )
    cross_energy = float(np.sum((weighted_b @ bmat) * cross_kernel) / n)
    predicted_energy = linear_energy + cross_energy
    predicted_gain = float(np.mean(derivative * variances))
    old_energy = first**2 * float(np.trace(b3)) / (2 * n)
    rng = np.random.default_rng(seed)
    energies, gains, original = [], [], []
    for start in range(0, samples, batch):
        count = min(batch, samples - start)
        xi = rng.normal(size=(count, n))
        field = xi @ bmat
        old = kappa * np.sin(field)
        returned = old @ bmat
        response = height * np.sin(frequency * returned)
        energies.extend(np.sum(response * (response @ bmat), axis=1) / (2 * n))
        gains.extend(np.sum(response * returned, axis=1) / n)
        original.extend(np.sum(old * returned, axis=1) / (2 * n))

    def report(values, prediction):
        values = np.asarray(values)
        mean = float(np.mean(values))
        return {"mean": mean, "standard_error": float(np.std(values, ddof=1) / np.sqrt(len(values))),
                "predicted_asymptotic": prediction, "finite_difference": mean - prediction}

    return {"n": n, "samples": samples, "seed": seed, "kappa": kappa,
            "frequency": frequency, "linear_feedback_term": linear_energy,
            "retained_bare_star": cross_energy, "feedback": report(energies, predicted_energy),
            "gain": report(gains, predicted_gain), "old_energy": report(original, old_energy),
            "asymptotic_endpoint_certificate": predicted_gain + abs(old_energy + predicted_energy),
            "evidence": "Monte Carlo diagnostics; finite-size asymptotic error is not a sampling error"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steiner-m", nargs="+", type=int, default=[3, 4, 5])
    parser.add_argument("--random-orders", nargs="*", type=int, default=[64, 128, 256])
    parser.add_argument("--samples", type=int, default=12000)
    parser.add_argument("--batch", type=int, default=300)
    parser.add_argument("--seed", type=int, default=9060837)
    parser.add_argument("--kappa", type=float, default=0.6)
    parser.add_argument("--frequency", type=float, default=1.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for m in args.steiner_m:
        a, _ = make_signing(m)
        record = evaluate(a, args.samples, args.batch, args.seed + m, args.kappa, args.frequency)
        record["family"] = "exact Steiner"
        records.append(record)
        print(json.dumps(record), flush=True)
    rng = np.random.default_rng(args.seed)
    for n in args.random_orders:
        upper = np.triu(rng.choice([-1, 1], size=(n, n)), 1)
        a = upper + upper.T
        record = evaluate(a, args.samples, args.batch, args.seed + n, args.kappa, args.frequency)
        record["family"] = "random signing, not a minimizer"
        records.append(record)
        print(json.dumps(record), flush=True)
    if args.output:
        args.output.write_text(json.dumps({"records": records}, indent=2) + "\n")


if __name__ == "__main__":
    main()
