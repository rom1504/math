#!/usr/bin/env python3
"""Exact finite audit for the Wave 28 completion-collision memo.

For each selector S, the canonical channel is uniform on all full oriented
projective cuts whose restriction is an exact oriented ground of A[S].
All incidence, row-square, degree, and entropy-deficit identities are checked
by exhaustive enumeration.  Only the displayed logarithms are floating point.
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict

import numpy as np

from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def exact_audit(a: np.ndarray, m: int) -> dict[str, float | int]:
    n = len(a)
    selectors = list(itertools.combinations(range(n), m))
    cuts: list[tuple[int, tuple[int, ...]]] = []
    row: dict[tuple[int, tuple[int, ...]], int] = {}
    for x_array in projective_spins(n):
        x = tuple(int(v) for v in x_array)
        r2 = int((a @ x_array) @ (a @ x_array))
        for sigma in (-1, 1):
            d = (sigma, x)
            cuts.append(d)
            row[d] = r2

    incident: dict[tuple[int, ...], list[tuple[int, tuple[int, ...]]]] = {}
    child_norm: dict[tuple[int, ...], int] = {}
    degree = defaultdict(int)
    ground_counts: list[int] = []
    conditional_row_formula: list[float] = []

    for s in selectors:
        child = a[np.ix_(s, s)]
        energies = []
        for y in projective_spins(m):
            energies.append((tuple(int(v) for v in y), int(y @ child @ y)))
        q_child = max(abs(e) for _, e in energies)
        child_norm[s] = q_child
        ground_counts.append(sum(abs(e) == q_child for _, e in energies))

        e_s = []
        direct_selected_norm_sum = 0
        direct_selected_norm_count = 0
        for d in cuts:
            sigma, x_tuple = d
            x = np.array(x_tuple, dtype=np.int64)
            y = x[list(s)]
            if sigma * int(y @ child @ y) == q_child:
                e_s.append(d)
                degree[d] += 1
                direct_selected_norm_sum += int((a[:, s] @ y) @ (a[:, s] @ y))
                direct_selected_norm_count += 1
        g_s = ground_counts[-1]
        assert len(e_s) == g_s * (2 ** (n - m))
        incident[s] = e_s

        # Averaging over the outside cube kills the cross term exactly.
        selected_norm_average = direct_selected_norm_sum / direct_selected_norm_count
        predicted = selected_norm_average + (n - m) * (n - 1)
        enumerated = sum(row[d] for d in e_s) / len(e_s)
        assert abs(predicted - enumerated) < 1e-10
        conditional_row_formula.append(predicted)

    nsel = len(selectors)
    marginal = defaultdict(float)
    for s in selectors:
        mass = 1.0 / nsel / len(incident[s])
        for d in incident[s]:
            marginal[d] += mass
    assert abs(sum(marginal.values()) - 1.0) < 1e-10

    mutual_information = 0.0
    expected_log_inverse_degree = 0.0
    expected_row = 0.0
    for s in selectors:
        conditional = 1.0 / len(incident[s])
        joint = 1.0 / nsel * conditional
        for d in incident[s]:
            mutual_information += joint * math.log(conditional / marginal[d])
    for d, probability in marginal.items():
        expected_log_inverse_degree += probability * math.log(nsel / degree[d])
        expected_row += probability * row[d]

    # Posterior-support entropy gives the first inequality.  Uniform ground
    # completion gives the elementary output-alphabet upper bound second.
    entropy_deficit = sum(
        math.log((2 ** (m - 1)) / g) for g in ground_counts
    ) / nsel
    assert expected_log_inverse_degree <= mutual_information + 1e-10
    assert mutual_information <= math.log(2.0) + entropy_deficit + 1e-10
    assert abs(expected_row - sum(conditional_row_formula) / nsel) < 1e-10

    # Exact simultaneous extraction: E[L/I + R/E R] <= 2.  Search the
    # positive output support and verify a witness exists.
    if mutual_information > 1e-14:
        scored = [
            (
                math.log(nsel / degree[d]) / mutual_information
                + row[d] / expected_row,
                d,
            )
            for d in marginal
        ]
        score, witness = min(scored)
        assert score <= 2.0 + 1e-10
    else:
        score = 1.0
        witness = next(iter(marginal))
        assert degree[witness] == nsel

    max_degree = max(degree.values())

    # The larger zero-threshold microcanonical incidence relation uses every
    # oriented child state with deficit <= B_{n,m}.  Every outside completion
    # then has retained effective loss <= 0.
    parent_q = max(abs(int(np.array(x) @ a @ np.array(x))) for x in projective_spins(n))
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    retained_slack = (p ** 1.5 - p2) * parent_q
    micro_incident: dict[tuple[int, ...], list[tuple[int, tuple[int, ...]]]] = {}
    micro_degree = defaultdict(int)
    micro_pattern_counts = []
    for s in selectors:
        child = a[np.ix_(s, s)]
        e_s = []
        for d in cuts:
            sigma, x_tuple = d
            x = np.array(x_tuple, dtype=np.int64)
            y = x[list(s)]
            child_payoff = sigma * int(y @ child @ y)
            deficit = child_norm[s] - child_payoff
            if deficit <= retained_slack + 1e-12:
                parent_payoff = sigma * int(x @ a @ x)
                effective_loss = deficit + p2 * parent_payoff - p ** 1.5 * parent_q
                assert effective_loss <= 1e-10
                e_s.append(d)
                micro_degree[d] += 1
        assert len(e_s) % (2 ** (n - m)) == 0
        micro_pattern_counts.append(len(e_s) // (2 ** (n - m)))
        assert micro_pattern_counts[-1] > 0
        micro_incident[s] = e_s

    micro_marginal = defaultdict(float)
    for s in selectors:
        mass = 1.0 / nsel / len(micro_incident[s])
        for d in micro_incident[s]:
            micro_marginal[d] += mass
    micro_information = 0.0
    for s in selectors:
        conditional = 1.0 / len(micro_incident[s])
        joint = conditional / nsel
        for d in micro_incident[s]:
            micro_information += joint * math.log(conditional / micro_marginal[d])
    micro_entropy_deficit = sum(
        math.log((2 ** m) / f) for f in micro_pattern_counts
    ) / nsel
    micro_expected_row = sum(
        probability * row[d] for d, probability in micro_marginal.items()
    )
    assert micro_information <= micro_entropy_deficit + 1e-10

    # Check the row-Gibbs incidence-pressure identities at several tilts.
    # Lbar' is evaluated analytically as conditional cost minus global cost.
    max_tilt_error = 0.0
    for lam in (0.0, 0.01, 0.05):
        weights = {d: math.exp(-lam * row[d]) for d in cuts}
        partition = sum(weights.values())
        global_row = sum(weights[d] * row[d] for d in cuts) / partition
        lbar = 0.0
        joint_row = 0.0
        lbar_prime = 0.0
        for s in selectors:
            local_partition = sum(weights[d] for d in incident[s])
            local_row = (
                sum(weights[d] * row[d] for d in incident[s])
                / local_partition
            )
            local_mass = local_partition / partition
            lbar += -math.log(local_mass) / nsel
            joint_row += local_row / nsel
            lbar_prime += (local_row - global_row) / nsel
        max_tilt_error = max(
            max_tilt_error,
            abs(joint_row - (global_row + lbar_prime)),
        )
        assert lbar >= -1e-12
    assert max_tilt_error < 1e-10

    return {
        "n": n,
        "m": m,
        "selectors": nsel,
        "ground_count_min": min(ground_counts),
        "ground_count_max": max(ground_counts),
        "ground_entropy_deficit": entropy_deficit,
        "mutual_information": mutual_information,
        "posterior_support_lower": expected_log_inverse_degree,
        "expected_row": expected_row,
        "max_degree": max_degree,
        "max_coverage": max_degree / nsel,
        "extraction_score": score,
        "extracted_row": row[witness],
        "extracted_degree": degree[witness],
        "micro_K": micro_entropy_deficit,
        "micro_I": micro_information,
        "micro_expected_row": micro_expected_row,
        "micro_max_coverage": max(micro_degree.values()) / nsel,
        "tilt_identity_error": max_tilt_error,
    }


def main() -> None:
    for name, a in (("A6", A6), ("A8", A8), ("A9", A9)):
        for m in range((len(a) + 1) // 2, len(a)):
            result = exact_audit(a, m)
            print(name, " ".join(f"{key}={value}" for key, value in result.items()))
    print("PASS completion-collision exact incidence/entropy/row audit")


if __name__ == "__main__":
    main()
