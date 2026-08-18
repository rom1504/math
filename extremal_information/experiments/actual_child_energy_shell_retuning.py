#!/usr/bin/env python3
"""Exact finite audit of the actual latent energy-shell retuning split.

The program exhaustively selects contracted-temperature minimizing children,
enumerates the complete bridge cube, and evaluates the averaged ordinary
forward-channel posterior under an inverse-disorder bridge tilt.  It then
checks the exact shell/geometric KL chain rule of Theorem 37.62.

Only actual minimizing children are used.  Integer child energies, bridge
words, and shell labels are exact; Gibbs probabilities and KL quantities are
floating-point evaluations.  The default balanced child order four keeps the
complete calculation small enough to reproduce quickly.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import actual_child_bridge_law_exact as exact  # noqa: E402


def logsumexp(values: np.ndarray) -> float:
    peak = float(np.max(values))
    return peak + math.log(float(np.sum(np.exp(values - peak))))


def normalize_log_weights(log_weights: np.ndarray) -> np.ndarray:
    shifted = np.asarray(log_weights, dtype=np.float64) - float(
        np.max(log_weights)
    )
    weights = np.exp(shifted)
    return weights / float(np.sum(weights))


def kl(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    positive = left > 0.0
    if np.any(right[positive] <= 0.0):
        return math.inf
    return float(np.dot(left[positive], np.log(left[positive] / right[positive])))


def latent_prior(
    left: np.ndarray,
    right: np.ndarray,
    beta: float,
    epsilon: int,
) -> dict:
    m, n = len(left), len(right)
    total_n = m + n
    t = beta / math.sqrt(total_n)
    x = exact.projective_spins(m).astype(np.int16)
    y = exact.projective_spins(n).astype(np.int16)
    ex = exact.energies_for_matrix(left, x).astype(np.int64)
    ey = exact.energies_for_matrix(right, y).astype(np.int64)
    d = m * n
    full_mask = (1 << d) - 1

    projective_patterns: list[int] = []
    combined_energies: list[int] = []
    projective_weights: list[float] = []
    for xi, exi in zip(x, ex):
        for yj, eyj in zip(y, ey):
            signs = (xi[:, None] * yj[None, :]).reshape(-1)
            pattern = 0
            for bit, sign in enumerate(signs):
                if sign < 0:
                    pattern |= 1 << bit
            projective_patterns.append(pattern)
            energy = int(exi + epsilon * eyj)
            combined_energies.append(energy)
            projective_weights.append(math.cosh(t * energy))

    if len(set(projective_patterns)) != len(projective_patterns):
        raise AssertionError("projective rank-one words collided")
    projective_probability = np.asarray(projective_weights, dtype=np.float64)
    projective_probability /= float(np.sum(projective_probability))

    patterns = np.empty(2 * len(projective_patterns), dtype=np.uint64)
    probability = np.empty(len(patterns), dtype=np.float64)
    shell = np.empty(len(patterns), dtype=np.int64)
    for index, (pattern, energy, weight) in enumerate(
        zip(projective_patterns, combined_energies, projective_probability)
    ):
        patterns[2 * index] = pattern
        patterns[2 * index + 1] = pattern ^ full_mask
        probability[2 * index : 2 * index + 2] = 0.5 * weight
        # Signed antipodes have the same quadratic combined-energy label.
        shell[2 * index : 2 * index + 2] = energy

    if len(set(int(value) for value in patterns)) != len(patterns):
        raise AssertionError("signed rank-one words collided")
    if abs(float(np.sum(probability)) - 1.0) > 2e-14:
        raise AssertionError(float(np.sum(probability)))
    return {
        "patterns": patterns,
        "probability": probability,
        "shell": shell,
        "raw_t": t,
    }


def shell_index(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, inverse = np.unique(np.asarray(labels, dtype=np.int64), return_inverse=True)
    return values, inverse


def aggregate(probability: np.ndarray, inverse: np.ndarray, count: int) -> np.ndarray:
    return np.bincount(inverse, weights=probability, minlength=count).astype(
        np.float64
    )


def audit_pair(
    left: np.ndarray,
    right: np.ndarray,
    beta: float,
    epsilon: int,
    inverse_disorder: float,
    chunk_size: int,
) -> dict:
    prior_record = latent_prior(left, right, beta, epsilon)
    patterns = prior_record["patterns"]
    prior = prior_record["probability"]
    labels = prior_record["shell"]
    t = float(prior_record["raw_t"])
    m, n = len(left), len(right)
    d = m * n
    bridge_count = 1 << d
    bridge_masks = np.arange(bridge_count, dtype=np.uint64)
    shell_values, inverse = shell_index(labels)
    shell_count = len(shell_values)
    prior_shell = aggregate(prior, inverse, shell_count)
    shell_indicator = np.zeros((len(prior), shell_count), dtype=np.float64)
    shell_indicator[np.arange(len(prior)), inverse] = 1.0

    # First pass: exact channel output likelihood and inverse bridge tilt.
    output = np.empty(bridge_count, dtype=np.float64)
    normalizer = math.cosh(t) ** d
    for start in range(0, bridge_count, chunk_size):
        stop = min(bridge_count, start + chunk_size)
        xor = np.bitwise_xor(bridge_masks[start:stop, None], patterns[None, :])
        dot = d - 2 * np.bitwise_count(xor).astype(np.int16)
        channel = np.exp(t * dot.astype(np.float64)) / normalizer
        output[start:stop] = channel @ prior
    if abs(float(np.mean(output)) - 1.0) > 2e-12:
        raise AssertionError(("output normalization", float(np.mean(output))))
    bridge_law = normalize_log_weights(-inverse_disorder * np.log(output))

    # Second pass: averaged posterior and posterior information terms.
    posterior_average = np.zeros(len(prior), dtype=np.float64)
    for start in range(0, bridge_count, chunk_size):
        stop = min(bridge_count, start + chunk_size)
        xor = np.bitwise_xor(bridge_masks[start:stop, None], patterns[None, :])
        dot = d - 2 * np.bitwise_count(xor).astype(np.int16)
        channel = np.exp(t * dot.astype(np.float64)) / normalizer
        posterior = channel * prior[None, :] / output[start:stop, None]
        posterior_average += bridge_law[start:stop] @ posterior
    posterior_average /= float(np.sum(posterior_average))
    posterior_shell = aggregate(posterior_average, inverse, shell_count)

    total_retuning = kl(posterior_average, prior)
    shell_retuning = kl(posterior_shell, prior_shell)
    within_shell = 0.0
    for shell_id in range(shell_count):
        members = inverse == shell_id
        mass = float(posterior_shell[shell_id])
        if mass <= 0.0:
            continue
        conditional = posterior_average[members] / mass
        uniform = np.full(np.count_nonzero(members), 1.0 / np.count_nonzero(members))
        within_shell += mass * kl(conditional, uniform)

    total_information = 0.0
    shell_information = 0.0
    posterior_work = 0.0
    for start in range(0, bridge_count, chunk_size):
        stop = min(bridge_count, start + chunk_size)
        xor = np.bitwise_xor(bridge_masks[start:stop, None], patterns[None, :])
        dot = d - 2 * np.bitwise_count(xor).astype(np.int16)
        channel = np.exp(t * dot.astype(np.float64)) / normalizer
        posterior = channel * prior[None, :] / output[start:stop, None]
        posterior_e = posterior @ shell_indicator
        q_chunk = bridge_law[start:stop]
        total_information += float(
            np.sum(
                q_chunk[:, None]
                * posterior
                * np.log(posterior / posterior_average[None, :])
            )
        )
        shell_information += float(
            np.sum(
                q_chunk[:, None]
                * posterior_e
                * np.log(posterior_e / posterior_shell[None, :])
            )
        )
        posterior_work += float(
            np.sum(
                q_chunk[:, None]
                * posterior
                * np.log(posterior / prior[None, :])
            )
        )

    return {
        "left_order": m,
        "right_order": n,
        "parent_order": m + n,
        "bridge_edge_count": d,
        "bridge_cube_size": bridge_count,
        "rank_one_word_count": len(prior),
        "combined_energy_shell_count": shell_count,
        "combined_energy_shell_values": shell_values.astype(int).tolist(),
        "beta": beta,
        "raw_t": t,
        "orientation": epsilon,
        "inverse_disorder": inverse_disorder,
        "prior_maximum_atom": float(np.max(prior)),
        "prior_minimum_shell_mass": float(np.min(prior_shell)),
        "posterior_maximum_atom": float(np.max(posterior_average)),
        "total_latent_retuning_KL": total_retuning,
        "energy_shell_retuning_KL": shell_retuning,
        "within_shell_geometric_KL": within_shell,
        "shell_chain_rule_residual": total_retuning
        - shell_retuning
        - within_shell,
        "posterior_mutual_information_total": total_information,
        "posterior_mutual_information_energy_shell": shell_information,
        "posterior_work_Eq_D_nuB_parallel_prior": posterior_work,
        "posterior_budget_residual": posterior_work
        - total_information
        - total_retuning,
        "fractions_when_total_retuning_positive": {
            "energy_shell": shell_retuning / total_retuning
            if total_retuning > 0.0
            else 0.0,
            "within_shell": within_shell / total_retuning
            if total_retuning > 0.0
            else 0.0,
        },
        "classification": (
            "exact signing/energy/bridge enumeration; floating-point Gibbs and KL"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child-order", type=int, default=4)
    parser.add_argument("--left-order", type=int)
    parser.add_argument("--right-order", type=int)
    parser.add_argument("--betas", nargs="+", type=float, default=[1.0, 2.0, 4.0])
    parser.add_argument(
        "--lambdas", nargs="+", type=float, default=[1.0, 5.382104195764755]
    )
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_energy_shell_retuning.json",
    )
    args = parser.parse_args()
    left_order = args.left_order if args.left_order is not None else args.child_order
    right_order = (
        args.right_order if args.right_order is not None else args.child_order
    )
    if min(left_order, right_order) < 2:
        parser.error("child orders must be at least two")
    if left_order * right_order > 24:
        parser.error("complete bridge cube is intentionally limited to 24 edges")
    if any(value <= 0 for value in args.betas + args.lambdas):
        parser.error("beta and lambda values must be positive")
    mp.mp.dps = args.mp_dps

    started = time.time()
    total_n = left_order + right_order
    left_space = exact.build_signing_space(left_order)
    right_space = (
        left_space
        if right_order == left_order
        else exact.build_signing_space(right_order)
    )
    records: list[dict] = []
    certificates: list[dict] = []
    for beta in args.betas:
        left_classes, left_certificate = exact.thermal_minimizer_classes(
            left_space, format(beta, ".17g"), total_n
        )
        right_classes, right_certificate = exact.thermal_minimizer_classes(
            right_space, format(beta, ".17g"), total_n
        )
        certificates.extend([left_certificate, right_certificate])
        for left_class in left_classes:
            left = np.asarray(left_class["representative_matrix"], dtype=np.int8)
            for right_class in right_classes:
                right = np.asarray(right_class["representative_matrix"], dtype=np.int8)
                for epsilon in (-1, 1):
                    for inverse_disorder in args.lambdas:
                        record = audit_pair(
                            left,
                            right,
                            beta,
                            epsilon,
                            inverse_disorder,
                            args.chunk_size,
                        )
                        record["left_class_id"] = left_class["class_id"]
                        record["right_class_id"] = right_class["class_id"]
                        records.append(record)
                        print(
                            f"m+n={left_order}+{right_order} beta={beta:g} "
                            f"eps={epsilon:+d} "
                            f"lambda={inverse_disorder:g} "
                            f"D={record['total_latent_retuning_KL']:.8g} "
                            f"shell={record['energy_shell_retuning_KL']:.8g} "
                            f"geo={record['within_shell_geometric_KL']:.8g}",
                            flush=True,
                        )

    payload = {
        "schema": "actual-child-energy-shell-retuning-v1",
        "parameters": {
            "child_order": args.child_order,
            "left_order": left_order,
            "right_order": right_order,
            "parent_order": total_n,
            "betas": args.betas,
            "lambdas": args.lambdas,
            "mp_dps": args.mp_dps,
            "chunk_size": args.chunk_size,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "mpmath": mp.__version__,
            "elapsed_seconds": time.time() - started,
        },
        "child_minimizer_certificates": certificates,
        "records": records,
        "evidentiary_status": (
            "finite exact-enumeration experiment; no asymptotic inference"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
