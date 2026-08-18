#!/usr/bin/env python3
"""Exhaustive finite audit of projective Gibbs concentration in actual children.

For each child order and parent-scaled temperature, this program enumerates
every root-gauged signing, selects every contracted-temperature pressure
minimizer, and evaluates the projective and augmented-projective Gibbs laws.
It then evaluates every pair of minimizing classes and both relative
orientations.  No conference, Paley, or other surrogate signing is loaded.

Integer energies and signing enumeration are exact.  Selection between
transcendental pressure histograms uses mpmath at the requested precision;
reported probabilities and entropies are floating-point evaluations.
"""

from __future__ import annotations

import argparse
import gc
import itertools
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
    maximum = float(np.max(values))
    return maximum + math.log(float(np.sum(np.exp(values - maximum))))


def log_cosh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    absolute = np.abs(values)
    return absolute + np.log1p(np.exp(-2.0 * absolute)) - math.log(2.0)


def probabilities_from_logits(logits: np.ndarray) -> np.ndarray:
    flat = np.asarray(logits, dtype=np.float64).reshape(-1)
    shifted = flat - float(np.max(flat))
    probability = np.exp(shifted)
    probability /= float(np.sum(probability))
    return probability


def distribution_summary(probability: np.ndarray) -> dict:
    probability = np.asarray(probability, dtype=np.float64).reshape(-1)
    if np.any(probability <= 0.0):
        raise AssertionError("all finite-temperature atoms must be positive")
    if abs(float(np.sum(probability)) - 1.0) > 2e-12:
        raise AssertionError(float(np.sum(probability)))
    log_probability = np.log(probability)
    shannon = -float(np.dot(probability, log_probability))
    collision = float(np.dot(probability, probability))
    renyi2 = -math.log(collision)
    count = len(probability)
    maximum = float(np.max(probability))
    descending = np.sort(probability)[::-1]
    return {
        "atom_count": count,
        "maximum_atom": maximum,
        "minimum_atom": float(np.min(probability)),
        "maximum_tie_count": int(
            np.count_nonzero(
                np.isclose(probability, maximum, rtol=2e-13, atol=2e-15)
            )
        ),
        "min_entropy_nats": -math.log(maximum),
        "shannon_entropy_nats": shannon,
        "shannon_entropy_fraction_of_full": shannon / math.log(count)
        if count > 1
        else 1.0,
        "shannon_entropy_deficit_nats": math.log(count) - shannon,
        "effective_shannon_atom_count": math.exp(shannon),
        "collision_probability": collision,
        "renyi2_entropy_nats": renyi2,
        "renyi2_entropy_fraction_of_full": renyi2 / math.log(count)
        if count > 1
        else 1.0,
        "renyi2_entropy_deficit_nats": math.log(count) - renyi2,
        "effective_collision_atom_count": 1.0 / collision,
        "effective_collision_fraction_of_full": 1.0 / (collision * count),
        "top_mass": {
            str(k): float(np.sum(descending[: min(k, count)]))
            for k in (1, 2, 4, 8, 16)
        },
    }


def absolute_histogram(energies: np.ndarray) -> dict[str, int]:
    absolute = np.abs(np.asarray(energies, dtype=np.int64))
    return {
        str(int(value)): int(np.count_nonzero(absolute == value))
        for value in np.unique(absolute)
    }


def shell_atom_lower_bound(energies: np.ndarray, t: float) -> dict:
    """A certified lower bound using only the top absolute shell and its gap."""

    absolute = np.abs(np.asarray(energies, dtype=np.int64).reshape(-1))
    levels = np.unique(absolute)
    cap = int(levels[-1])
    top_count = int(np.count_nonzero(absolute == cap))
    count = len(absolute)
    if len(levels) == 1:
        second = None
        ratio = 0.0
        lower = 1.0 / top_count
    else:
        second = int(levels[-2])
        ratio = math.cosh(t * second) / math.cosh(t * cap)
        lower = 1.0 / (top_count + (count - top_count) * ratio)
    return {
        "absolute_cap": cap,
        "top_absolute_shell_count": top_count,
        "second_absolute_level": second,
        "top_shell_gap": None if second is None else cap - second,
        "non_top_to_top_weight_ratio_bound": ratio,
        "projective_max_atom_shell_lower_bound": lower,
    }


def child_record(matrix: np.ndarray, beta: float, parent_order: int) -> dict:
    order = len(matrix)
    spins = exact.projective_spins(order)
    energies = exact.energies_for_matrix(matrix, spins).astype(np.int64)
    t = beta / math.sqrt(parent_order)
    projective_logits = log_cosh(t * energies)
    projective = probabilities_from_logits(projective_logits)
    augmented_logits = np.stack((t * energies, -t * energies), axis=0)
    augmented = probabilities_from_logits(augmented_logits)
    projective_summary = distribution_summary(projective)
    augmented_summary = distribution_summary(augmented)
    cap = int(np.max(np.abs(energies)))
    exact_projective_max = math.cosh(t * cap) / float(
        np.sum(np.cosh(t * energies))
    )
    exact_augmented_max = math.exp(t * cap) / float(
        2.0 * np.sum(np.cosh(t * energies))
    )
    relation = exact_projective_max / (1.0 + math.exp(-2.0 * t * cap))
    if abs(exact_projective_max - projective_summary["maximum_atom"]) > 2e-13:
        raise AssertionError("child projective atom identity failed")
    if abs(exact_augmented_max - augmented_summary["maximum_atom"]) > 2e-13:
        raise AssertionError("child augmented atom identity failed")
    if abs(exact_augmented_max - relation) > 2e-13:
        raise AssertionError("child augmented/projective relation failed")
    return {
        "order": order,
        "parent_order": parent_order,
        "beta": beta,
        "raw_t": t,
        "energy_histogram_projective": absolute_histogram(energies),
        "projective_cosh_gibbs": projective_summary,
        "augmented_projective_gibbs_s_xclass": augmented_summary,
        "atom_identities": {
            "projective_max_cosh_tK_over_sum_cosh": exact_projective_max,
            "augmented_max_exp_tK_over_2_sum_cosh": exact_augmented_max,
            "augmented_max_from_projective_max": relation,
            "maximum_identity_residual": max(
                abs(exact_projective_max - projective_summary["maximum_atom"]),
                abs(exact_augmented_max - augmented_summary["maximum_atom"]),
            ),
        },
        "shell_bound": shell_atom_lower_bound(energies, t),
        "energies": energies,
    }


def pair_record(
    left: dict,
    right: dict,
    epsilon: int,
    beta: float,
    parent_order: int,
) -> dict:
    t = beta / math.sqrt(parent_order)
    left_energies = np.asarray(left["energies"], dtype=np.float64)
    right_energies = np.asarray(right["energies"], dtype=np.float64)
    energies = (
        left_energies[:, None] + epsilon * right_energies[None, :]
    ).astype(np.int64)
    projective = probabilities_from_logits(log_cosh(t * energies))
    augmented = probabilities_from_logits(
        np.stack((t * energies, -t * energies), axis=0)
    )
    projective_summary = distribution_summary(projective)
    augmented_summary = distribution_summary(augmented)
    cap = int(np.max(np.abs(energies)))
    denominator = float(np.sum(np.cosh(t * energies)))
    exact_projective_max = math.cosh(t * cap) / denominator
    exact_augmented_max = math.exp(t * cap) / (2.0 * denominator)
    relation = exact_projective_max / (1.0 + math.exp(-2.0 * t * cap))
    if abs(exact_projective_max - projective_summary["maximum_atom"]) > 2e-13:
        raise AssertionError("pair projective atom identity failed")
    if abs(exact_augmented_max - augmented_summary["maximum_atom"]) > 2e-13:
        raise AssertionError("pair augmented atom identity failed")
    if abs(exact_augmented_max - relation) > 2e-13:
        raise AssertionError("pair augmented/projective relation failed")

    # Independently reconstruct the laws from the exact sector factorization
    # pi_s mu_(A,s) tensor mu_(D,epsilon*s).  This guards the normalization
    # and confirms that summing the two sectors gives the cosh law above.
    sector_joints = []
    sector_weights = []
    for sector in (1, -1):
        left_logits = t * sector * left_energies
        right_logits = t * epsilon * sector * right_energies
        left_probability = probabilities_from_logits(left_logits)
        right_probability = probabilities_from_logits(right_logits)
        left_partition = math.exp(logsumexp(left_logits) - math.log(len(left_logits)))
        right_partition = math.exp(
            logsumexp(right_logits) - math.log(len(right_logits))
        )
        sector_weights.append(left_partition * right_partition)
        sector_joints.append(
            left_probability[:, None] * right_probability[None, :]
        )
    sector_weights = np.asarray(sector_weights, dtype=np.float64)
    sector_weights /= float(np.sum(sector_weights))
    augmented_from_sectors = np.concatenate(
        [
            (sector_weights[index] * sector_joints[index]).reshape(-1)
            for index in range(2)
        ]
    )
    projective_from_sectors = (
        sector_weights[0] * sector_joints[0]
        + sector_weights[1] * sector_joints[1]
    ).reshape(-1)
    sector_residual = max(
        float(np.max(np.abs(augmented - augmented_from_sectors))),
        float(np.max(np.abs(projective - projective_from_sectors))),
    )
    if sector_residual > 3e-13:
        raise AssertionError(("sector/projective identity failed", sector_residual))
    signed_collision = 0.5 * projective_summary["collision_probability"]
    return {
        "epsilon": epsilon,
        "projective_rank_one_qclass": projective_summary,
        "augmented_projective_gibbs_s_xclass_yclass": augmented_summary,
        "signed_rank_one_Q": {
            "atom_count": 2 * projective_summary["atom_count"],
            "maximum_atom": 0.5 * projective_summary["maximum_atom"],
            "collision_probability": signed_collision,
            "renyi2_entropy_nats": -math.log(signed_collision),
            "effective_collision_atom_count": 1.0 / signed_collision,
            "shannon_entropy_nats": (
                projective_summary["shannon_entropy_nats"] + math.log(2.0)
            ),
        },
        "two_word_common_sign_mass": projective_summary["maximum_atom"],
        "atom_identities": {
            "projective_max_cosh_tK_over_sum_cosh": exact_projective_max,
            "augmented_max_exp_tK_over_2_sum_cosh": exact_augmented_max,
            "augmented_max_from_projective_max": relation,
            "signed_max_half_projective_max": (
                0.5 * projective_summary["maximum_atom"]
            ),
            "signed_collision_half_projective_collision": signed_collision,
            "sector_factorization_maximum_residual": sector_residual,
            "maximum_identity_residual": max(
                abs(exact_projective_max - projective_summary["maximum_atom"]),
                abs(exact_augmented_max - augmented_summary["maximum_atom"]),
            ),
        },
        "shell_bound": shell_atom_lower_bound(energies, t),
        "combined_energy_histogram_projective_pairs": absolute_histogram(energies),
    }


def strip_internal(record: dict) -> dict:
    return {key: value for key, value in record.items() if key != "energies"}


def extrema(records: list[dict], path: tuple[str, ...]) -> dict:
    def extract(record: dict) -> float:
        value = record
        for key in path:
            value = value[key]
        return float(value)

    values = [extract(record) for record in records]
    return {"minimum": min(values), "maximum": max(values)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", type=int, default=list(range(2, 9)))
    parser.add_argument("--betas", nargs="+", type=float, default=[1.0, 2.0, 4.0])
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_projective_atom_entropy.json",
    )
    args = parser.parse_args()
    if args.mp_dps <= 20:
        parser.error("--mp-dps must exceed the selector's 20 guard digits")
    if any(order < 2 for order in args.orders):
        parser.error("orders must be at least two")
    mp.mp.dps = args.mp_dps

    started = time.time()
    records = []
    certificates = []
    for order in args.orders:
        build_started = time.time()
        space = exact.build_signing_space(order)
        print(
            f"built order {order}: {len(space.energies)} root-gauged signings "
            f"in {time.time()-build_started:.2f}s",
            flush=True,
        )
        parent_order = 2 * order
        for beta in args.betas:
            classes, certificate = exact.thermal_minimizer_classes(
                space, format(beta, ".17g"), parent_order
            )
            certificates.append(certificate)
            children = []
            for item in classes:
                matrix = np.asarray(item["representative_matrix"], dtype=np.int16)
                child = child_record(matrix, beta, parent_order)
                child.update(
                    {
                        "class_id": item["class_id"],
                        "representative_sha256": item["representative_sha256"],
                        "root_gauged_member_count": item["root_gauged_member_count"],
                    }
                )
                children.append(child)
            pairs = []
            for left, right, epsilon in itertools.product(children, children, (-1, 1)):
                pair = pair_record(left, right, epsilon, beta, parent_order)
                pair.update(
                    {
                        "left_class_id": left["class_id"],
                        "right_class_id": right["class_id"],
                        "left_sha256": left["representative_sha256"],
                        "right_sha256": right["representative_sha256"],
                    }
                )
                pairs.append(pair)
            aggregate = {
                "child_projective_max_atom": extrema(
                    children, ("projective_cosh_gibbs", "maximum_atom")
                ),
                "child_augmented_projective_max_atom": extrema(
                    children,
                    ("augmented_projective_gibbs_s_xclass", "maximum_atom"),
                ),
                "child_projective_effective_collision_count": extrema(
                    children,
                    (
                        "projective_cosh_gibbs",
                        "effective_collision_atom_count",
                    ),
                ),
                "pair_projective_max_atom": extrema(
                    pairs, ("projective_rank_one_qclass", "maximum_atom")
                ),
                "pair_augmented_projective_max_atom": extrema(
                    pairs,
                    (
                        "augmented_projective_gibbs_s_xclass_yclass",
                        "maximum_atom",
                    ),
                ),
                "pair_projective_effective_collision_count": extrema(
                    pairs,
                    (
                        "projective_rank_one_qclass",
                        "effective_collision_atom_count",
                    ),
                ),
                "pair_projective_shannon_entropy_nats": extrema(
                    pairs,
                    ("projective_rank_one_qclass", "shannon_entropy_nats"),
                ),
                "pair_projective_renyi2_entropy_nats": extrema(
                    pairs,
                    ("projective_rank_one_qclass", "renyi2_entropy_nats"),
                ),
            }
            records.append(
                {
                    "child_order": order,
                    "parent_order": parent_order,
                    "beta": beta,
                    "raw_t": beta / math.sqrt(parent_order),
                    "thermal_minimizer_class_count": len(classes),
                    "aggregate": aggregate,
                    "children": [strip_internal(child) for child in children],
                    "pairs": pairs,
                }
            )
            print(
                f"m={order} beta={beta:g} classes={len(classes)} "
                f"pair-max={aggregate['pair_projective_max_atom']['maximum']:.8g} "
                f"pair-Neff2-min={aggregate['pair_projective_effective_collision_count']['minimum']:.8g}",
                flush=True,
            )
        del space
        gc.collect()

    payload = {
        "schema": "actual-child-projective-atom-entropy-v1",
        "classification": (
            "finite exhaustive child-signing selection; exact integer energies; "
            "high-precision numerical pressure comparison; floating-point Gibbs summaries"
        ),
        "definitions": {
            "child_projective": "[x] with weight proportional to cosh(t H_A(x))",
            "child_augmented_projective": (
                "(s,[x]) with weight proportional to exp(t s H_A(x))"
            ),
            "pair_projective_rank_one": (
                "[Q]=[xy^T], equivalently ([x],[y]), conditional on epsilon, "
                "with weight proportional to cosh(t(H_A(x)+epsilon H_D(y)))"
            ),
            "pair_augmented_projective": (
                "(s,[x],[y]) with weight proportional to "
                "exp(t s(H_A(x)+epsilon H_D(y)))"
            ),
            "signed_Q_fibre": (
                "conditional on [Q], Q and -Q have equal mass; signed max is "
                "half the projective max and signed collision is half the projective collision"
            ),
            "effective_collision_atom_count": "1/sum_atom p(atom)^2",
        },
        "parameters": {
            "orders": args.orders,
            "betas": args.betas,
            "parent_order_rule": "N=2m",
            "raw_temperature_rule": "t=beta/sqrt(N)",
            "mp_dps": args.mp_dps,
            "all_thermal_minimizer_classes": True,
            "all_class_pairs": True,
            "all_relative_orientations": True,
            "surrogate_signings": False,
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "elapsed_seconds": time.time() - started,
        },
        "child_minimizer_certificates": certificates,
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
