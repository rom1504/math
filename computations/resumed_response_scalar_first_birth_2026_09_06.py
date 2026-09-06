#!/usr/bin/env python3
"""Floating diagnostic: one exact unrestricted gradient birth and line search.

All expectations reduce to one outer Gaussian integral. This is NOT an
exact lower certificate; tails and quadrature need separate auditing.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import ndtr

from resumed_response_scalar_full_optimizer_2026_09_06 import ScalarObjective, phi, hermites


def gaussian_absolute(k, variance):
    tau = np.sqrt(variance)
    z = k / tau
    return 2 * tau * phi(z) + k * (2 * ndtr(z) - 1)


def run(parameters, nodes, cutoff):
    degree = 2 * (len(parameters) - 1)
    alpha = parameters[0]
    objective = ScalarObjective(degree, max(160, nodes // 2))
    value, info = objective.evaluate(parameters)
    coefficients = np.array(info["coefficients_even"])
    base_nodes, base_weights = np.polynomial.legendre.leggauss(nodes)
    center_v = alpha * (base_nodes + 1) / 2
    center_weights = alpha * base_weights * phi(center_v)
    tail_v = alpha + (cutoff - alpha) * (base_nodes + 1) / 2
    tail_weights = (cutoff - alpha) * base_weights * phi(tail_v)
    v = np.r_[center_v, tail_v]
    weights = np.r_[center_weights, tail_weights]
    H = np.r_[np.ones(nodes), np.zeros(nodes)]
    g = hermites(v, degree)[:, ::2] @ coefficients
    lam, gam, variance = info["lambda"], info["gamma"], info["tau"] ** 2
    K = lam * g + gam * H
    psi = gaussian_absolute(K, variance)
    s = 2 * ndtr(K / np.sqrt(variance)) - 1
    B = float(weights @ (H * phi(K / np.sqrt(variance)))) / np.sqrt(variance)
    u = H * (s - 2 * B * gam)
    ug = float(weights @ (u * g))
    normu2 = float(weights @ (u * u))
    A1 = ug - 2 * B * lam
    nu = np.sqrt(max(0.0, normu2 - ug * ug))
    if nu <= 0:
        raise ValueError("Degenerate gradient noise")
    hU = (u - ug * g) / nu
    mean = A1 * v
    threshold = np.maximum(psi - B, 0)
    plus = ndtr((mean - threshold) / nu)
    minus = ndtr((-mean - threshold) / nu)
    q_new = plus + minus
    f_mean = plus - minus
    root_moment = phi((threshold - mean) / nu) + phi((threshold + mean) / nu)
    p_new = float(weights @ q_new)
    aV_new = float(weights @ (v * f_mean))
    aU_new = float(weights @ root_moment)
    K_new = aV_new * g + aU_new * hU
    q_old = 1 - info["p"]
    normK_old2 = q_old - variance
    HhU = (float(weights @ u) - ug * info["c"]) / nu
    cross_K = (aV_new * (lam + info["c"] * gam)
               + aU_new * gam * HhU)
    normK_new2 = aV_new * aV_new + aU_new * aU_new

    def line(theta):
        support = (1 - theta) * q_old + theta * p_new
        normK2 = ((1 - theta) ** 2 * normK_old2
                  + 2 * theta * (1 - theta) * cross_K
                  + theta**2 * normK_new2)
        residual = support - normK2
        if residual <= 0:
            raise ValueError("Invalid line residual")
        Ktheta = (1 - theta) * K + theta * K_new
        Htheta = (1 - theta) * H + theta * (1 - q_new)
        return float(weights @ (Htheta * gaussian_absolute(Ktheta, residual)))

    opt = minimize_scalar(lambda theta: -line(theta), bounds=(0, 1), method="bounded",
                          options={"xatol": 1e-12})
    candidates = [(0.0, line(0)), (1.0, line(1)), (float(opt.x), -float(opt.fun))]
    theta_best, value_best = max(candidates, key=lambda item: item[1])
    old_inner = (info["lambda"] * info["c"] + info["gamma"] * info["p"])
    # Direct full gap, integrating the Gaussian positive part analytically.
    b = psi - B
    nonnegative_b = np.maximum(b, 0)
    positive_part = (
        nu * (phi((nonnegative_b - mean) / nu) + phi((nonnegative_b + mean) / nu))
        + (mean - nonnegative_b) * ndtr((mean - nonnegative_b) / nu)
        + (-mean - nonnegative_b) * ndtr((-mean - nonnegative_b) / nu)
    )
    positive_part = np.where(b >= 0, positive_part, gaussian_absolute(mean, nu * nu) - b)
    gap = float(weights @ (psi + positive_part)) + q_old * B - 2 * line(0)
    finite_difference_gap = (line(1e-5) - line(0)) / 1e-5
    return {
        "status": "floating_diagnostic_not_certificate",
        "degree": degree,
        "nodes_each_interval": nodes,
        "cutoff": cutoff,
        "alpha": alpha,
        "derivative_energy": float(parameters[1:] @ parameters[1:]),
        "old_scalar_value": value,
        "old_integrated_value": line(0),
        "old_marked_value": old_inner,
        "B": B,
        "gradient_A1": A1,
        "gradient_nu": nu,
        "full_gap": gap,
        "line_derivative_check": finite_difference_gap,
        "new_support": p_new,
        "new_aV": aV_new,
        "new_aU": aU_new,
        "full_step_value": line(1),
        "best_theta": theta_best,
        "best_value": value_best,
        "gain": value_best - line(0),
        "line_values": {str(theta): line(theta) for theta in np.linspace(0, 1, 11)},
        "parameters": parameters.tolist(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--nodes", type=int, default=384)
    parser.add_argument("--cutoff", type=float, default=12.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    candidates = json.loads(args.candidate.read_text())["runs"]
    candidate = max(candidates, key=lambda item: item["value"])
    parameters = np.array(candidate["parameters"])
    # Keep a visible floating feasibility margin; this is still not an exact audit.
    parameters[1:] *= min(1.0, (1 - 1e-8) / np.linalg.norm(parameters[1:]))
    first = run(parameters, args.nodes, args.cutoff)
    replay = run(parameters, 2 * args.nodes, args.cutoff)
    result = {"first": first, "double_node_replay": replay,
              "best_value_difference": replay["best_value"] - first["best_value"]}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({key: replay[key] for key in [
        "old_scalar_value", "full_gap", "line_derivative_check", "best_theta",
        "best_value", "gain", "full_step_value"]}))
    print("double-node difference", result["best_value_difference"])


if __name__ == "__main__":
    main()
