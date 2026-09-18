#!/usr/bin/env python3
"""Numerical falsification checks for the control/affine phase theorem.

The proof is analytic in the director artifact. These finite floating-point
checks are diagnostics, not optimality or asymptotic certificates.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def fwht(values):
    result = np.array(values, dtype=complex, copy=True)
    width = 1
    while width < result.size:
        rows = result.reshape(-1, 2 * width)
        left, right = rows[:, :width].copy(), rows[:, width:].copy()
        rows[:, :width], rows[:, width:] = left + right, left - right
        width *= 2
    return result / result.size


def entropy(probabilities):
    values = np.asarray(probabilities)
    positive = values[values > 0]
    return float(-np.sum(positive * np.log(positive)))


def one_case(rng, k, n, degree, time):
    z = np.arange(1 << k)
    characters = np.array(
        [[(-1) ** bin(int(v) & mask).count("1") for mask in range(1 << k)]
         for v in z], dtype=float)
    coefficients = rng.normal(size=(n + 1, 1 << k))
    for row in range(n + 1):
        budget = degree if row == 0 else degree - 1
        coefficients[row, [mask for mask in range(1 << k)
                           if bin(mask).count("1") > budget]] = 0
    fields = coefficients @ characters.T
    cap = float(np.max(np.abs(fields[0]) + np.sum(np.abs(fields[1:]), axis=0)))
    fields /= cap
    coefficients /= cap
    x = np.array([[(-1) ** ((int(v) >> i) & 1) for i in range(n)]
                  for v in range(1 << n)])
    energies = fields[0][None, :] + x @ fields[1:]
    full_coefficients = fwht(np.exp(1j * time * energies).reshape(-1))
    joint = np.abs(full_coefficients.reshape(1 << n, 1 << k)) ** 2
    assert abs(float(joint.sum()) - 1) < 1e-11
    marginal = joint.sum(axis=1)
    amplitudes = np.ones((1 << n, 1 << k), dtype=complex)
    for mask in range(1 << n):
        for i in range(n):
            angle = time * fields[i + 1]
            amplitudes[mask] *= (1j * np.sin(angle) if (mask >> i) & 1
                                 else np.cos(angle))
    conditional = np.abs(amplitudes) ** 2
    assert np.max(np.abs(conditional.sum(axis=0) - 1)) < 1e-11
    assert np.max(np.abs(marginal - conditional.mean(axis=1))) < 1e-11
    conditional_entropy = sum(entropy(conditional[:, i]) for i in z) / len(z)
    information = entropy(marginal) - conditional_entropy
    derivative_energy = 0.0
    field_influence = 0.0
    for j in range(k):
        flips = z ^ (1 << j)
        derivative_energy += float(np.sum(np.abs(
            (amplitudes - amplitudes[:, flips]) / 2) ** 2) / len(z))
        field_influence += float(np.sum(
            ((fields[1:] - fields[1:, flips]) / 2) ** 2) / len(z))
    total_entropy = entropy(joint)
    bound = k * np.log(2) + 2 * abs(time) + 2 * time ** 2 * (degree - 1)
    linear_bound = (6 + 2 * k * np.log(2)) * abs(time)
    tolerance = 1e-9
    assert conditional_entropy <= 2 * abs(time) + tolerance
    assert information <= 2 * derivative_energy + tolerance
    assert derivative_energy <= time ** 2 * field_influence + tolerance
    assert field_influence <= degree - 1 + tolerance
    assert total_entropy <= entropy(marginal) + k * np.log(2) + tolerance
    assert total_entropy <= bound + tolerance
    assert total_entropy <= linear_bound + tolerance
    if abs(time) <= 1:
        assert np.max(1 - conditional[0]) <= time ** 2 + tolerance
        control_mass = float(joint[:, 1:].sum())
        assert control_mass <= time ** 2 + tolerance
        assert information <= 2 * abs(time) + time ** 2 * k * np.log(2) + tolerance
    return {"control_bits": k, "data_bits": n, "degree_budget": degree,
            "time_cap": time, "entropy": total_entropy, "bound": float(bound),
            "linear_bound": float(linear_bound),
            "mutual_information": information,
            "log_sobolev_bound": 2 * derivative_energy,
            "field_influence": field_influence}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=2026091806)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    cases = [one_case(rng, k, n, d, t)
             for k, n, d in [(0, 8, 1), (2, 6, 2), (3, 7, 3), (4, 8, 4)]
             for t in [0.0, 0.001, 0.1, 1.0, 3.0]]
    report = {"status": "PASS", "evidence": "floating-point diagnostics only",
              "seed": args.seed, "cases": cases}
    output = json.dumps(report, indent=2)
    print(output)
    if args.output:
        args.output.write_text(output + "\n")


if __name__ == "__main__":
    main()
