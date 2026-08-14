#!/usr/bin/env python3
"""Exhaustive small-order audit of the joint reverse-KL bridge identity.

The normalized partition function is

    Zbar_n(A,t) = 2^{-n} sum_x cosh(t H_A(x)).

Switching symmetry is fixed by taking every edge incident with vertex 0 to
be +1.  This gives one representative of every switching orbit.  The audit
is exhaustive over those representatives, both relative child orientations,
and every bridge for total order at most ``--max-n``.

Floating-point evaluations are exploratory, not formal certificates.  The
enumerations themselves are complete at the requested orders.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def spins(n: int) -> np.ndarray:
    return np.asarray(list(itertools.product((-1, 1), repeat=n)), dtype=np.int8)


def canonical_signings(n: int) -> list[np.ndarray]:
    """One signing per switching orbit, with a_0j=+1."""
    if n <= 1:
        return [np.zeros((n, n), dtype=np.int8)]
    free_edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    result: list[np.ndarray] = []
    for mask in range(1 << len(free_edges)):
        a = np.zeros((n, n), dtype=np.int8)
        a[0, 1:] = 1
        a[1:, 0] = 1
        for bit, (i, j) in enumerate(free_edges):
            value = -1 if (mask >> bit) & 1 else 1
            a[i, j] = value
            a[j, i] = value
        result.append(a)
    return result


def energies(a: np.ndarray, x: np.ndarray) -> np.ndarray:
    if a.shape[0] <= 1:
        return np.zeros(x.shape[0], dtype=np.int64)
    return np.einsum("bi,ij,bj->b", x, a, x, dtype=np.int64) // 2


def log_mean_cosh(values: np.ndarray) -> float:
    magnitudes = np.abs(np.asarray(values, dtype=np.float64))
    peak = float(np.max(magnitudes))
    scaled = 0.5 * (
        np.exp(magnitudes - peak) + np.exp(-magnitudes - peak)
    )
    return peak + math.log(float(np.mean(scaled)))


def pressure(energy: np.ndarray, raw_t: float) -> float:
    return log_mean_cosh(raw_t * energy)


def minimizers(
    all_energies: list[np.ndarray], raw_t: float, tolerance: float = 1e-11
) -> tuple[float, list[int]]:
    values = [pressure(energy, raw_t) for energy in all_energies]
    optimum = min(values)
    indices = [i for i, value in enumerate(values) if value <= optimum + tolerance]
    return optimum, indices


def bridge_matrices(m: int, n: int):
    for mask in range(1 << (m * n)):
        entries = np.ones(m * n, dtype=np.int8)
        for bit in range(m * n):
            if (mask >> bit) & 1:
                entries[bit] = -1
        yield entries.reshape((m, n))


def mean_parent_log_pressure(
    a_energy: np.ndarray,
    d_energy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    raw_t: float,
) -> float:
    total = 0.0
    count = 0
    internal_a = a_energy[:, None]
    internal_d = d_energy[None, :]
    for bridge in bridge_matrices(x.shape[1], y.shape[1]):
        cross = x.astype(np.int64) @ bridge.astype(np.int64) @ y.astype(np.int64).T
        for orientation in (-1, 1):
            parent_energy = internal_a + orientation * internal_d + cross
            total += pressure(parent_energy.reshape(-1), raw_t)
            count += 1
    return total / count


def audit(max_n: int, betas: list[float]) -> dict:
    spin_cache = {n: spins(n) for n in range(1, max_n + 1)}
    signing_cache = {n: canonical_signings(n) for n in range(1, max_n + 1)}
    energy_cache = {
        n: [energies(a, spin_cache[n]) for a in signing_cache[n]]
        for n in range(1, max_n + 1)
    }

    payload: dict = {
        "status": "exhaustive enumeration with floating-point evaluation",
        "normalization": "Zbar=2^{-n} sum_x cosh(t H_A(x))",
        "max_total_order": max_n,
        "betas": betas,
        "records": [],
    }

    for beta in betas:
        high_optima = {
            n: minimizers(energy_cache[n], beta / math.sqrt(n))[0]
            for n in range(1, max_n + 1)
        }
        for total_n in range(2, max_n + 1):
            raw_t = beta / math.sqrt(total_n)
            parent_optimum = high_optima[total_n]
            for m in range(1, total_n // 2 + 1):
                n = total_n - m
                low_m, minimizer_m = minimizers(energy_cache[m], raw_t)
                low_n, minimizer_n = minimizers(energy_cache[n], raw_t)
                best_average = math.inf
                best_pair = None
                for i in minimizer_m:
                    for j in minimizer_n:
                        average = mean_parent_log_pressure(
                            energy_cache[m][i],
                            energy_cache[n][j],
                            spin_cache[m],
                            spin_cache[n],
                            raw_t,
                        )
                        if average < best_average:
                            best_average = average
                            best_pair = [i, j]

                bridge_annealed = m * n * math.log(math.cosh(raw_t))
                reverse_kl = low_m + low_n + bridge_annealed - best_average
                thermal_gap = (
                    high_optima[m] - low_m + high_optima[n] - low_n
                )
                compensation = thermal_gap + reverse_kl
                interface_defect = bridge_annealed - compensation
                actual_same_beta_defect = (
                    parent_optimum - high_optima[m] - high_optima[n]
                )
                payload["records"].append(
                    {
                        "beta": beta,
                        "N": total_n,
                        "split": [m, n],
                        "raw_t": raw_t,
                        "low_child_minimizer_counts": [
                            len(minimizer_m),
                            len(minimizer_n),
                        ],
                        "best_low_child_pair_indices": best_pair,
                        "bridge_annealed": bridge_annealed,
                        "thermal_gap": thermal_gap,
                        "reverse_kl": reverse_kl,
                        "joint_compensation": compensation,
                        "interface_defect": interface_defect,
                        "interface_defect_per_N": interface_defect / total_n,
                        "actual_same_beta_defect": actual_same_beta_defect,
                        "actual_same_beta_defect_per_N": (
                            actual_same_beta_defect / total_n
                        ),
                    }
                )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=6)
    parser.add_argument(
        "--betas", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit(args.max_n, args.betas)
    rendered = json.dumps(payload, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
