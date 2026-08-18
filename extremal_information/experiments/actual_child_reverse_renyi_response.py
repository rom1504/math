#!/usr/bin/env python3
"""Exact finite verifier for the reverse-Renyi response identity.

The child signings are complete-enumeration thermal minimizers.  Every bridge
and every leave-one-coordinate channel response is then enumerated exactly;
floating-point transcendental evaluation is numerical.  The script verifies
RR.7 by centered finite differences and records the response density and
leave-one-out dependence of the actual inverse escort.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_row_product_shadow as shadow  # noqa: E402
import verify_actual_child_row_anova as anova  # noqa: E402


def likelihood(
    left: np.ndarray,
    right: np.ndarray,
    internal_t: float,
    bridge_u: float,
    epsilon: int,
) -> np.ndarray:
    pressure = anova.generalized_bridge_pressures(
        left, right, internal_t, bridge_u, epsilon
    )
    shifted = pressure - float(np.max(pressure))
    p = np.exp(shifted)
    p /= float(np.mean(p))
    return p


def response_metrics(p: np.ndarray, rows: int, columns: int, u: float, lam: float) -> dict:
    d = rows * columns
    rho = math.tanh(u)
    logp = np.log(p)
    logq = -lam * logp
    logq -= float(np.max(logq))
    q = np.exp(logq)
    q /= float(np.sum(q))

    indices = np.arange(len(p), dtype=np.uint64)
    total_r2 = np.zeros(len(p), dtype=np.float64)
    total_j = np.zeros(len(p), dtype=np.float64)
    total_mi = np.zeros(len(p), dtype=np.float64)
    total_s2 = np.zeros(len(p), dtype=np.float64)
    maximum_formula_residual = 0.0

    if rho == 0.0:
        return {
            "u": u,
            "rho": rho,
            "reverse_renyi": 0.0,
            "analytic_derivative": 0.0,
            "S": 0.0,
            "S_density": 0.0,
            "leave_one_out_mutual_information_sum": 0.0,
            "leave_one_out_bias_square_sum": 0.0,
            "maximum_conditional_bias_formula_residual": 0.0,
        }

    for e in range(d):
        bit = np.uint64(1 << e)
        zero = indices[(indices & bit) == 0]
        one = zero | bit
        # Our mask convention uses bit zero for bridge sign +1 and bit one
        # for bridge sign -1.  The sign choice only changes r, not r^2 or J.
        p_plus = p[zero]
        p_minus = p[one]
        r = (p_plus - p_minus) / (rho * (p_plus + p_minus))
        r = np.clip(r, -1.0, 1.0)
        a = np.clip(rho * r, -1.0 + 1e-15, 1.0 - 1e-15)
        s_formula = -np.tanh(lam * np.arctanh(a))

        q_plus = q[zero]
        q_minus = q[one]
        pair_mass = q_plus + q_minus
        s_direct = (q_plus - q_minus) / pair_mass
        # Because mask zero is B_e=+1, this is the q conditional mean.
        maximum_formula_residual = max(
            maximum_formula_residual,
            float(np.max(np.abs(s_direct - s_formula))),
        )

        j = (
            (1.0 - rho * rho)
            * r
            * (rho * r + np.tanh(lam * np.arctanh(a)))
            / (1.0 - a * a)
        )
        mi = 0.5 * (
            (1.0 + s_direct) * np.log1p(s_direct)
            + (1.0 - s_direct) * np.log1p(-s_direct)
        )

        # r, j, and the conditional information are constant on each pair.
        total_r2[zero] += r * r
        total_r2[one] += r * r
        total_j[zero] += j
        total_j[one] += j
        total_mi[zero] += mi
        total_mi[one] += mi
        total_s2[zero] += s_direct * s_direct
        total_s2[one] += s_direct * s_direct

    reverse_renyi = math.log(float(np.mean(p ** (-lam)))) / lam
    reverse_kl = -float(np.mean(logp))
    row_size = 1 << columns
    tensor = shadow.pressure_tensor(p, rows, row_size)
    row_likelihoods = []
    for row in range(rows):
        axes = tuple(axis for axis in range(rows) if axis != row)
        row_likelihoods.append(np.mean(tensor, axis=axes))
    row_works = [
        math.log(float(np.mean(row_p ** (-lam)))) / lam
        for row_p in row_likelihoods
    ]
    row_escorts = [
        row_p ** (-lam) / float(np.sum(row_p ** (-lam)))
        for row_p in row_likelihoods
    ]
    product_likelihood = np.ones_like(tensor)
    product_escort = np.ones_like(tensor)
    for row in range(rows):
        shape = [1] * rows
        shape[row] = row_size
        product_likelihood *= row_likelihoods[row].reshape(shape)
        product_escort *= row_escorts[row].reshape(shape)
    interaction = np.log(tensor / product_likelihood)
    interaction_mean = float(np.sum(product_escort * interaction))
    centered_cumulant = math.log(
        float(
            np.sum(
                product_escort
                * np.exp(-lam * (interaction - interaction_mean))
            )
        )
    )
    return {
        "u": u,
        "rho": rho,
        "reverse_renyi": reverse_renyi,
        "reverse_KL": reverse_kl,
        "centered_negative_gain": reverse_renyi - reverse_kl,
        "canonical_row_erased_sum_work": sum(row_works),
        "canonical_row_interaction_mean": interaction_mean,
        "canonical_row_centered_interaction_cumulant": centered_cumulant,
        "canonical_row_product_inverse_work": sum(row_works) - interaction_mean,
        "canonical_row_certificate_reverse_KL": centered_cumulant,
        "canonical_decomposition_residual": (
            reverse_renyi
            - sum(row_works)
            + interaction_mean
            - centered_cumulant / lam
        ),
        "analytic_derivative": float(np.dot(q, total_j)),
        "S": float(np.dot(q, total_r2)),
        "S_density": float(np.dot(q, total_r2)) / d,
        "leave_one_out_mutual_information_sum": float(np.dot(q, total_mi)),
        "leave_one_out_bias_square_sum": float(np.dot(q, total_s2)),
        "maximum_conditional_bias_formula_residual": maximum_formula_residual,
    }


def run(args: argparse.Namespace) -> dict:
    mp.mp.dps = 80
    N = args.total_order
    m = args.left_order
    n = N - m
    internal_t = args.beta / math.sqrt(N)
    spaces = {k: exact.build_signing_space(k) for k in {m, n}}
    children = {}
    for k in (m, n):
        classes = exact.thermal_minimizer_classes(
            spaces[k], format(args.beta, ".12g"), N
        )[0]
        children[k] = np.asarray(classes[0]["representative_matrix"], dtype=np.int8)
    left = children[m]
    right = children[n]

    records = []
    for epsilon in args.epsilons:
        for lam in args.lambdas:
            laws = []
            for fraction in args.amplitude_fractions:
                u = fraction * internal_t
                p = likelihood(left, right, internal_t, u, epsilon)
                item = response_metrics(p, m, n, u, lam)
                step = args.finite_difference_step
                if u > step:
                    p_minus = likelihood(left, right, internal_t, u - step, epsilon)
                    p_plus = likelihood(left, right, internal_t, u + step, epsilon)
                    r_minus = math.log(float(np.mean(p_minus ** (-lam)))) / lam
                    r_plus = math.log(float(np.mean(p_plus ** (-lam)))) / lam
                    item["centered_finite_difference_derivative"] = (
                        r_plus - r_minus
                    ) / (2.0 * step)
                    item["derivative_residual"] = (
                        item["analytic_derivative"]
                        - item["centered_finite_difference_derivative"]
                    )
                laws.append(item)
            records.append(
                {
                    "epsilon": epsilon,
                    "lambda": lam,
                    "amplitudes": laws,
                }
            )

    return {
        "schema": "actual-child-reverse-renyi-response-v1",
        "classification": (
            "complete finite child and bridge enumeration; numerical "
            "transcendental evaluation; no asymptotic inference"
        ),
        "theorem": "extremal_information/drafts/actual_child_reverse_renyi_response_identity.md",
        "N": N,
        "split": [m, n],
        "beta": args.beta,
        "internal_t": internal_t,
        "left_sha": exact.matrix_sha(left),
        "right_sha": exact.matrix_sha(right),
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--total-order", type=int, default=8)
    parser.add_argument("--left-order", type=int, default=4)
    parser.add_argument("--beta", type=float, default=4.0)
    parser.add_argument("--epsilons", type=int, nargs="+", default=[-1, 1])
    parser.add_argument("--lambdas", type=float, nargs="+", default=[1.0, 5.382104])
    parser.add_argument(
        "--amplitude-fractions",
        type=float,
        nargs="+",
        default=[0.25, 0.5, 0.75, 1.0],
    )
    parser.add_argument("--finite-difference-step", type=float, default=1e-5)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "computations/results/actual_child_reverse_renyi_response.json",
    )
    args = parser.parse_args()
    payload = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
