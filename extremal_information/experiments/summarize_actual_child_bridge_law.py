#!/usr/bin/env python3
"""Create a compact, reproducible summary of the actual-child exact audits."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def selected_law(record: dict, lam: float) -> dict:
    return next(law for law in record["laws"] if law["lambda"] == lam)


def audit(exact: dict, shadow: dict, threshold_shadow: dict | None = None) -> dict:
    exact_laws = [
        (record, law)
        for record in exact["records"]
        for law in record["laws"]
    ]
    shadow_laws = [
        (record, law)
        for record in shadow["records"]
        for law in record["laws"]
    ]

    sanity = {
        "maximum_direct_Walsh_vs_spin_log_pressure_error": max(
            record["bridge_pressure_audit"]["maximum_direct_log_pressure_error"]
            for record in exact["records"]
        ),
        "maximum_KL_escort_identity_residual": max(
            abs(law["KL_identity_residual"]) for _, law in exact_laws
        ),
        "minimum_Renyi2_minus_KL": min(
            law["Renyi2_q_parallel_U"] - law["KL_q_parallel_U"]
            for _, law in exact_laws
        ),
        "minimum_KL_minus_row_TC": min(
            law["KL_q_parallel_U"] - law["row_total_correlation"]
            for _, law in exact_laws
        ),
        "minimum_KL_minus_column_TC": min(
            law["KL_q_parallel_U"] - law["column_total_correlation"]
            for _, law in exact_laws
        ),
        "minimum_mean_gain_minus_soft_gain": min(
            law["mean_pressure_gain"] - law["negative_moment_pressure_gain"]
            for _, law in exact_laws
        ),
        "maximum_conditional_D2_fraction_of_AC10_bound": max(
            max(
                filtration["aggregate"]["maximum"]
                for filtration in law["row_filtrations"]
            )
            /
            (
                law["lambda"] ** 2
                * record["beta"] ** 2
                * record["split"][1]
                / record["N"]
            )
            for record, law in exact_laws
        ),
        "maximum_shadow_objective_increase_in_best_trace": max(
            max(
                later - earlier
                for earlier, later in zip(
                    law["best_run"]["objective_trace"],
                    law["best_run"]["objective_trace"][1:],
                )
            )
            for _, law in shadow_laws
        ),
        "shadow_best_runs_hitting_100_sweep_cap": sum(
            not law["best_run"]["converged_by_objective_tolerance"]
            for _, law in shadow_laws
        ),
        "shadow_terminal_objective_multiplicity": dict(
            sorted(
                Counter(
                    law["distinct_terminal_objective_count_1e-9"]
                    for _, law in shadow_laws
                ).items()
            )
        ),
    }

    certificates = exact["child_minimizer_certificates"]
    finite_selector = {
        "queried_child_temperature_contexts": len(certificates),
        "maximum_signed_permutation_global_sign_class_count": max(
            certificate["signed_permutation_global_sign_class_count"]
            for certificate in certificates
        ),
        "minimum_high_precision_gap_to_next_distinct_histogram": min(
            float(certificate["mp_gap_to_next_histogram"])
            for certificate in certificates
            if certificate["mp_gap_to_next_histogram"] is not None
        ),
        "maximum_root_gauged_signing_count_enumerated": max(
            certificate["root_gauged_signing_count"]
            for certificate in certificates
        ),
    }

    target_feasibility = []
    for total_n in range(4, 10):
        for beta in (0.25, 0.5, 1.0, 2.0, 4.0):
            candidates = [
                record
                for record in exact["records"]
                if record["N"] == total_n and record["beta"] == beta
            ]
            best = min(
                candidates,
                key=lambda record: record["bridge_pressure_minimum"]
                - record["same_temperature_minimum_child_target"],
            )
            target_feasibility.append(
                {
                    "N": total_n,
                    "beta": beta,
                    "best_split": best["split"],
                    "best_orientation": best["relative_child_orientation"],
                    "minimum_bridge_gap": best["bridge_pressure_minimum"]
                    - best["same_temperature_minimum_child_target"],
                    "minimum_bridge_gap_per_N": (
                        best["bridge_pressure_minimum"]
                        - best["same_temperature_minimum_child_target"]
                    )
                    / total_n,
                }
            )

    threshold_laws = []
    for record in exact["records"]:
        threshold = record["target_threshold"]
        if "law" not in threshold:
            continue
        law = threshold["law"]
        threshold_laws.append(
            {
                "N": record["N"],
                "split": record["split"],
                "beta": record["beta"],
                "orientation": record["relative_child_orientation"],
                "threshold_lambda": threshold["lambda"],
                "KL_per_N": law["KL_per_parent_vertex"],
                "row_TC_per_N": law["row_total_correlation"] / record["N"],
                "best_one_row_latent_residual_TC_per_N": law[
                    "best_latent_residual_tc"
                ]["value"]
                / record["N"],
                "best_filtration_maximum_step_mean_D2": law[
                    "best_row_filtration"
                ]["maximum_step_weighted_mean"],
                "chain_support_proxy_fraction": law[
                    "chain_support_proxy_fraction"
                ],
            }
        )

    inverse_proposal = []
    for total_n in range(4, 10):
        m = total_n // 2
        n = total_n - m
        for beta in (1.0, 2.0, 4.0):
            candidates = [
                record
                for record in exact["records"]
                if record["N"] == total_n
                and record["split"] == [m, n]
                and record["beta"] == beta
            ]
            record = min(
                candidates,
                key=lambda candidate: selected_law(candidate, 1.0)[
                    "KL_q_parallel_positive_output_law"
                ],
            )
            law = selected_law(record, 1.0)
            edge_count = m * n
            raw_t = beta / math.sqrt(total_n)
            log_k = (total_n - 1) * math.log(2.0)
            im2_bound = edge_count * math.log(math.cosh(raw_t)) - raw_t * math.sqrt(
                2 * edge_count * log_k
            )
            inverse_proposal.append(
                {
                    "N": total_n,
                    "split": [m, n],
                    "beta": beta,
                    "orientation": record["relative_child_orientation"],
                    "KL_U_parallel_p_per_N": law[
                        "KL_U_parallel_positive_output_law_per_parent_vertex"
                    ],
                    "KL_q_parallel_p_per_N": law[
                        "KL_q_parallel_positive_output_law_per_parent_vertex"
                    ],
                    "TV_q_p": law["TV_q_positive_output_law"],
                    "affinity_q_p": law["affinity_q_positive_output_law"],
                    "IM2_support_bound_per_N": im2_bound / total_n,
                }
            )

    shadow_balanced = []
    for total_n in range(4, 10):
        m = total_n // 2
        n = total_n - m
        for beta in (1.0, 2.0, 4.0):
            for record in shadow["records"]:
                if not (
                    record["N"] == total_n
                    and record["split"] == [m, n]
                    and record["beta"] == beta
                ):
                    continue
                law = selected_law(record, 1.0)
                shadow_balanced.append(
                    {
                        "N": total_n,
                        "split": [m, n],
                        "beta": beta,
                        "orientation": record["relative_child_orientation"],
                        "captured_gain_fraction": law["captured_gain_fraction"],
                        "candidate_product_gain_per_N": law[
                            "rigorous_candidate_row_product_gain"
                        ]
                        / total_n,
                        "candidate_reverse_projection_upper_bound_per_N": law[
                            "candidate_reverse_projection_upper_bound"
                        ]
                        / total_n,
                        "best_response_residual_L1": law["best_run"][
                            "simultaneous_best_response_maximum_l1"
                        ],
                        "terminal_objective_count_1e-9": law[
                            "distinct_terminal_objective_count_1e-9"
                        ],
                    }
                )

    threshold_shadow_summary = []
    if threshold_shadow is not None:
        for record in threshold_shadow["records"]:
            if record["split"] != [4, 4] or record["beta"] != 4.0:
                continue
            law = record["laws"][0]
            threshold_shadow_summary.append(
                {
                    "N": record["N"],
                    "split": record["split"],
                    "beta": record["beta"],
                    "orientation": record["relative_child_orientation"],
                    "lambda": law["lambda"],
                    "captured_gain_fraction": law["captured_gain_fraction"],
                    "candidate_product_gain_per_N": law[
                        "rigorous_candidate_row_product_gain"
                    ]
                    / record["N"],
                    "candidate_reverse_projection_upper_bound_per_N": law[
                        "candidate_reverse_projection_upper_bound"
                    ]
                    / record["N"],
                    "best_response_residual_L1": law["best_run"][
                        "simultaneous_best_response_maximum_l1"
                    ],
                    "terminal_objective_count_1e-9": law[
                        "distinct_terminal_objective_count_1e-9"
                    ],
                }
            )

    return {
        "schema": "actual-child-bridge-law-compact-summary-v1",
        "classification": (
            "derived exact-enumeration numerical summary; no asymptotic inference; "
            "row-product candidates are feasible but not globally certified"
        ),
        "sources": {
            "exact": "computations/results/actual_child_bridge_law_exact.json",
            "row_product_shadow": "computations/results/actual_child_row_product_shadow.json",
            "target_threshold_row_product_shadow": (
                "computations/results/actual_child_row_product_shadow_target_threshold_n8.json"
                if threshold_shadow is not None
                else None
            ),
        },
        "sanity_checks": sanity,
        "finite_thermal_child_selector": finite_selector,
        "best_all_split_target_feasibility": target_feasibility,
        "finite_target_threshold_laws": threshold_laws,
        "balanced_lambda1_inverse_proposal_comparison": inverse_proposal,
        "balanced_lambda1_row_product_shadow": shadow_balanced,
        "target_threshold_row_product_shadow": threshold_shadow_summary,
        "global_extrema": {
            "maximum_KL_per_N": max(
                law["KL_per_parent_vertex"] for _, law in exact_laws
            ),
            "maximum_row_TC_per_row": max(
                law["row_total_correlation_per_row"] for _, law in exact_laws
            ),
            "minimum_chain_support_proxy_fraction": min(
                law["chain_support_proxy_fraction"] for _, law in exact_laws
            ),
            "maximum_chain_support_proxy_fraction": max(
                law["chain_support_proxy_fraction"] for _, law in exact_laws
            ),
            "maximum_lambda1_TV_q_p": max(
                selected_law(record, 1.0)["TV_q_positive_output_law"]
                for record in exact["records"]
            ),
            "minimum_shadow_captured_gain_fraction": min(
                law["captured_gain_fraction"] for _, law in shadow_laws
            ),
            "maximum_shadow_captured_gain_fraction": max(
                law["captured_gain_fraction"] for _, law in shadow_laws
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--exact",
        type=Path,
        default=ROOT / "computations/results/actual_child_bridge_law_exact.json",
    )
    parser.add_argument(
        "--threshold-shadow",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_row_product_shadow_target_threshold_n8.json",
    )
    parser.add_argument(
        "--shadow",
        type=Path,
        default=ROOT / "computations/results/actual_child_row_product_shadow.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_bridge_law_summary.json",
    )
    args = parser.parse_args()
    payload = audit(
        json.loads(args.exact.read_text()),
        json.loads(args.shadow.read_text()),
        (
            json.loads(args.threshold_shadow.read_text())
            if args.threshold_shadow.exists()
            else None
        ),
    )
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
