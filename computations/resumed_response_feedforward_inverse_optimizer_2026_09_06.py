#!/usr/bin/env python3
"""Floating 1D shape search for a causal inverse above a fixed rich core.

No Sobolev/derivative constraint is imposed on the new inverse: it is a
feedforward function of the already actual Gaussian V. The score is the
conditional-V Jensen LOWER functional, not the full rich-core expectation.
All values here remain diagnostics until exact norm/moment/bin auditing.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr

from resumed_response_scalar_full_optimizer_2026_09_06 import hermites, phi, PHI0


class FeedforwardObjective:
    def __init__(self, coefficients, nodes=192, step=True):
        self.coefficients = coefficients
        self.degree = len(coefficients) - 1
        self.norm_m2 = float(coefficients @ coefficients)
        self.nodes, self.weights = np.polynomial.legendre.leggauss(nodes)
        self.step = step

    def evaluate(self, pars, details=False):
        alpha, y = float(pars[0]), np.asarray(pars[1:])
        x = np.r_[np.sqrt(1 - y @ y), y]
        p = 2 * ndtr(alpha) - 1
        p_small = 2 * ndtr(alpha / 2) - 1
        v0 = alpha * (self.nodes + 1) / 4
        v1 = alpha / 2 + alpha * (self.nodes + 1) / 4
        v = np.r_[v0, v1]
        weights = np.r_[alpha / 2 * self.weights * phi(v0),
                        alpha / 2 * self.weights * phi(v1)]
        small = np.r_[np.ones(len(v0)), np.zeros(len(v1))]
        m = hermites(v, self.degree) @ self.coefficients
        c_H = float(weights @ m)
        moment2_H = float(weights @ (m * m))
        var_center = moment2_H - c_H * c_H / p
        var_tail = self.norm_m2 - moment2_H
        e0 = np.ones(len(v)) / np.sqrt(p)
        e1 = (m - c_H / p) / np.sqrt(var_center)
        correlation = np.array([c_H / np.sqrt(p), np.sqrt(var_center), np.sqrt(var_tail)])
        u = x[0] * e0 + x[1] * e1
        step_norm = step_on_e0 = step_on_e1 = 0.0
        if self.step:
            step_on_e0 = p_small / np.sqrt(p)
            step_on_e1 = float(weights @ (small * e1))
            step_norm = np.sqrt(p_small - step_on_e0**2 - step_on_e1**2)
            e3 = (small - step_on_e0 * e0 - step_on_e1 * e1) / step_norm
            u += x[3] * e3
        c = float(correlation @ x[:3])
        d = np.sqrt(1 - c * c)
        lam = 2 * phi(alpha) * (2 * ndtr(c * alpha / d) - 1)
        gam = 4 * PHI0 * ndtr(-alpha / d)
        residual = 1 - p - lam * lam - 2 * c * lam * gam - gam * gam
        tau = np.sqrt(residual)
        K = lam * m + gam * u
        value = float(weights @ (2 * tau * phi(K / tau) + K * (2 * ndtr(K / tau) - 1)))
        if not details:
            return value
        return {
            "value": value,
            "alpha": alpha,
            "unit_basis_coordinates": x.tolist(),
            "basis_correlations": correlation.tolist(),
            "conditional_inverse_norm_squared": self.norm_m2,
            "center_mass": p,
            "center_inverse_mean_integral": c_H,
            "center_inverse_second_integral": moment2_H,
            "center_variance": var_center,
            "tail_second_integral": var_tail,
            "step_mass": p_small,
            "step_on_e0": step_on_e0,
            "step_on_e1": step_on_e1,
            "step_residual_norm": step_norm,
            "output_covariance": c,
            "lambda": lam,
            "gamma": gam,
            "residual_variance": residual,
            "parameters": pars.tolist(),
        }

    def objective_gradient(self, pars):
        step = 1e-5
        gradient = np.empty(len(pars))
        for j in range(len(pars)):
            left, right = pars.copy(), pars.copy()
            left[j] -= step
            right[j] += step
            gradient[j] = (self.evaluate(right) - self.evaluate(left)) / (2 * step)
        return -self.evaluate(pars), -gradient


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection", type=Path, default=Path(
        "computations/results/resumed_response_conditional_v_certificate_2026_09_06.json"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--nodes", type=int, default=128)
    parser.add_argument("--fixed-alpha", action="store_true")
    args = parser.parse_args()
    source = json.loads(args.projection.read_text())
    coefficients = np.array([float(z["lower"]) for z in source["conditional_projection_coefficients"]])
    records = []
    for use_step in [False, True]:
        objective = FeedforwardObjective(coefficients, args.nodes, step=use_step)
        start = np.r_[0.7246, np.zeros(3 if use_step else 2)]
        bounds = ([(0.7246, 0.7246)] if args.fixed_alpha else [(0.55, 0.9)])
        bounds += [(-0.3, 0.3)] * (len(start) - 1)
        result = minimize(objective.objective_gradient, start, jac=True,
                          method="L-BFGS-B", bounds=bounds,
                          options={"ftol": 2e-15, "gtol": 2e-10, "maxiter": 150})
        info = objective.evaluate(result.x, details=True)
        replay = FeedforwardObjective(coefficients, 3 * args.nodes, step=use_step).evaluate(result.x)
        info.update({"step_basis": use_step, "success": bool(result.success),
                     "message": str(result.message), "iterations": int(result.nit),
                     "replay_value": replay, "replay_difference": replay - info["value"]})
        records.append(info)
        print(json.dumps({key: info[key] for key in ["step_basis", "value", "alpha", "unit_basis_coordinates", "success", "replay_difference"]}), flush=True)
    output = {"status": "floating_diagnostic_not_certificate",
              "fixed_rich_core_projection": str(args.projection), "runs": records}
    if args.output:
        args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
