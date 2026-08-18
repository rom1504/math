#!/usr/bin/env python3
"""Finite checks of the reverse row-refresh de Bruijn identity.

This verifier uses a non-product random law against a nonuniform product
base, integrates the exact row-refresh Fisher functional, and compares it
with ``D(r||q)``.  It also checks the two closed-form benchmark endpoints:
a product target (zero interaction lifetime) and full parity (inverse-order
cost).  The experiment is an implementation audit, not evidence for an
asymptotic actual-child estimate.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[2]


def average_axis(values: np.ndarray, probability: np.ndarray, axis: int) -> np.ndarray:
    """Average ``axis`` and broadcast the result back to the full shape."""

    shape = [1] * values.ndim
    shape[axis] = len(probability)
    averaged = np.sum(values * probability.reshape(shape), axis=axis, keepdims=True)
    return np.broadcast_to(averaged, values.shape)


def refresh(values: np.ndarray, bases: list[np.ndarray], time: float) -> np.ndarray:
    retained = math.exp(-time)
    result = values.copy()
    for axis, probability in enumerate(bases):
        result = retained * result + (1.0 - retained) * average_axis(
            result, probability, axis
        )
    return result


def product_mass(bases: list[np.ndarray]) -> np.ndarray:
    result = np.ones(tuple(len(probability) for probability in bases))
    for axis, probability in enumerate(bases):
        shape = [1] * len(bases)
        shape[axis] = len(probability)
        result *= probability.reshape(shape)
    return result


def fisher(density: np.ndarray, bases: list[np.ndarray], base_mass: np.ndarray) -> float:
    total = 0.0
    reciprocal = 1.0 / density
    for axis, probability in enumerate(bases):
        mean = average_axis(density, probability, axis)
        inverse_mean = average_axis(reciprocal, probability, axis)
        # The integrand is constant along the averaged axis.  Averaging with
        # the full product mass is therefore exactly E_(r_-i).
        total += float(np.sum(base_mass * (mean * inverse_mean - 1.0)))
    return total


def random_case() -> dict:
    rng = np.random.default_rng(20260818)
    bases = [np.asarray([0.31, 0.69]), np.asarray([0.57, 0.43]), np.asarray([0.22, 0.78])]
    base = product_mass(bases)
    raw = np.exp(rng.normal(size=base.shape))
    target = raw / np.sum(raw)
    density = target / base
    exact = float(np.sum(base * np.log(base / target)))

    def integrand(time: float) -> float:
        return fisher(refresh(density, bases, time), bases, base)

    integral, error = quad(integrand, 0.0, 40.0, epsabs=2e-12, epsrel=2e-12)
    tail_bound = abs(integrand(40.0))
    return {
        "D_r_parallel_q": exact,
        "integrated_row_refresh_fisher": integral,
        "quadrature_error_estimate": error,
        "integrand_at_time_40": tail_bound,
        "absolute_identity_residual": abs(exact - integral),
    }


def product_case() -> dict:
    bases = [np.asarray([0.4, 0.6]), np.asarray([0.7, 0.3])]
    base = product_mass(bases)
    target_factors = [np.asarray([0.8, 0.2]), np.asarray([0.25, 0.75])]
    target = product_mass(target_factors)
    # Choosing r=q is the canonical base for an already row-product target.
    density = target / target
    value = fisher(density, target_factors, target)
    return {
        "maximum_density_deviation_from_one": float(np.max(np.abs(density - 1.0))),
        "row_refresh_fisher": value,
    }


def parity_case(bit_count: int = 12, amplitude: float = 0.37) -> dict:
    indices = np.arange(1 << bit_count, dtype=np.uint64)
    parity = 1.0 - 2.0 * (np.bitwise_count(indices) & 1).astype(float)
    theta = math.tanh(amplitude)
    density = (1.0 - theta * parity).reshape((2,) * bit_count)
    bases = [np.asarray([0.5, 0.5]) for _ in range(bit_count)]
    base = product_mass(bases)

    def integrand(time: float) -> float:
        return fisher(refresh(density, bases, time), bases, base)

    integral, error = quad(integrand, 0.0, 8.0, epsabs=2e-11, epsrel=2e-11)
    exact = math.log(math.cosh(amplitude))
    time_zero = bit_count * math.sinh(amplitude) ** 2
    return {
        "bit_count": bit_count,
        "amplitude": amplitude,
        "time_zero_raw_fisher": time_zero,
        "integrated_row_refresh_fisher": integral,
        "closed_form_log_cosh": exact,
        "quadrature_error_estimate": error,
        "absolute_identity_residual": abs(exact - integral),
    }


def main() -> None:
    payload = {
        "schema": "reverse-row-noise-lifetime-verification-v1",
        "classification": "finite numerical verification of exact identities",
        "random_nonproduct_case": random_case(),
        "product_case": product_case(),
        "parity_case": parity_case(),
    }
    output = ROOT / "computations/results/reverse_row_noise_lifetime_verify.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
