#!/usr/bin/env python3
"""Floating diagnostic optimization, NOT a certified lower bound.

Optimize the exact scalar full-center objective over even Hermite inverses
of norm one and derivative energy <=1. The nonconstant parameter x_r is
sqrt(r)*g_r, so the derivative constraint is the Euclidean unit ball.
All reported candidates require separate exact feasibility/value audit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr


PHI0 = 1 / np.sqrt(2 * np.pi)


def phi(x):
    return PHI0 * np.exp(-np.asarray(x) ** 2 / 2)


def hermites(v, degree):
    ans = np.empty((len(v), degree + 1), dtype=float)
    ans[:, 0] = 1
    if degree:
        ans[:, 1] = v
    for k in range(1, degree):
        ans[:, k + 1] = (v * ans[:, k] - np.sqrt(k) * ans[:, k - 1]) / np.sqrt(k + 1)
    return ans


class ScalarObjective:
    def __init__(self, degree, nodes):
        self.degree = degree
        self.degrees = np.arange(2, degree + 1, 2, dtype=float)
        self.nodes, self.weights = np.polynomial.legendre.leggauss(nodes)

    def evaluate(self, pars, gradient=False):
        alpha, x = float(pars[0]), np.asarray(pars[1:])
        nonconstant_norm = np.sum(x * x / self.degrees)
        if nonconstant_norm >= 0.999999:
            penalty = 1e3 + 1e3 * nonconstant_norm
            if gradient:
                return -penalty, np.r_[0.0, -2e3 * x / self.degrees], {}
            return -penalty, {}
        b0 = np.sqrt(1 - nonconstant_norm)
        coefficients = np.r_[b0, x / np.sqrt(self.degrees)]
        endpoint = hermites(np.array([alpha]), self.degree)[0]
        p = 2 * ndtr(alpha) - 1
        beta = np.r_[p, -2 * phi(alpha) * endpoint[1::2] / np.sqrt(self.degrees)]
        c = float(coefficients @ beta)
        d2 = p - c * c
        if d2 <= 0:
            raise ValueError("Degenerate covariance in numerical trial")
        d = np.sqrt(d2)
        z = c * alpha / d
        lam = 2 * phi(alpha) * (2 * ndtr(z) - 1)
        gam = 4 * PHI0 / np.sqrt(p) * ndtr(-alpha * np.sqrt(p) / d)
        variance = 1 - p - lam * lam - 2 * c * lam * gam - p * gam * gam
        if variance <= 0:
            raise ValueError("Nonpositive residual variance in numerical trial")
        tau = np.sqrt(variance)
        v = alpha * self.nodes
        weights = alpha * self.weights * phi(v)
        basis = hermites(v, self.degree)[:, ::2]
        g = basis @ coefficients
        k = lam * g + gam
        s = 2 * ndtr(k / tau) - 1
        psi = 2 * tau * phi(k / tau) + k * s
        value = float(weights @ psi)
        information = {
            "alpha": alpha,
            "degree": self.degree,
            "value": value,
            "derivative_energy": float(x @ x),
            "p": float(p),
            "c": c,
            "d": float(d),
            "lambda": float(lam),
            "gamma": float(gam),
            "tau": float(tau),
            "marked_value": float(c * lam + p * gam),
            "coefficients_even": coefficients.tolist(),
        }
        if not gradient:
            return value, information
        B = float(weights @ phi(k / tau)) / tau
        L = 4 * alpha * phi(alpha) * phi(z) / d**3
        lam_c, gam_c = L * p, -L * c
        variance_c = -2 * lam * (gam + L * d2)
        A = lam_c * float(weights @ (g * s)) + gam_c * float(weights @ s) + B * variance_c
        function_gradient = lam * s + A
        parameter_basis = basis[:, 1:] / np.sqrt(self.degrees) - x / (self.degrees * b0)
        grad_x = (weights * function_gradient) @ parameter_basis
        step = 2e-5
        left, right = pars.copy(), pars.copy()
        left[0] -= step
        right[0] += step
        grad_alpha = (self.evaluate(right)[0] - self.evaluate(left)[0]) / (2 * step)
        return value, np.r_[grad_alpha, grad_x], information

    def resolvent_start(self, alpha=0.74, resolvent=9.7):
        endpoint = hermites(np.array([alpha]), self.degree)[0]
        p = 2 * ndtr(alpha) - 1
        beta = np.r_[p, -2 * phi(alpha) * endpoint[1::2] / np.sqrt(self.degrees)]
        coeff = beta / (resolvent + np.r_[0.0, self.degrees])
        coeff /= np.linalg.norm(coeff)
        x = coeff[1:] * np.sqrt(self.degrees)
        if np.linalg.norm(x) > 1:
            x *= 0.999 / np.linalg.norm(x)
        return np.r_[alpha, x]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--degrees", default="20,40,80")
    parser.add_argument("--nodes", type=int, default=160)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--maxiter", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260906)
    parser.add_argument("--broad-starts", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    records = []
    previous = None
    for degree in map(int, args.degrees.split(",")):
        if degree < 2 or degree % 2:
            raise ValueError("Degrees must be positive even integers")
        objective = ScalarObjective(degree, args.nodes)
        starts = [objective.resolvent_start()]
        if previous is not None:
            extended = np.zeros(1 + degree // 2)
            extended[:min(len(previous), len(extended))] = previous[:len(extended)]
            starts[0] = extended
        for _ in range(args.starts - 1):
            if args.broad_starts:
                trial = np.r_[rng.uniform(0.3, 1.2), rng.normal(size=degree // 2)]
                trial[1:] *= rng.uniform(0.4, 0.995) / np.linalg.norm(trial[1:])
            else:
                trial = objective.resolvent_start(alpha=rng.uniform(0.55, 0.95))
                trial[1:] += rng.normal(scale=0.09, size=degree // 2) / np.sqrt(objective.degrees)
                trial[1:] *= min(1.0, 0.995 / np.linalg.norm(trial[1:]))
            starts.append(trial)
        best = None
        for index, start in enumerate(starts):
            def fun(pars):
                val, grad, _ = objective.evaluate(pars, gradient=True)
                return -val, -grad

            result = minimize(
                fun, start, jac=True, method="SLSQP",
                bounds=[(0.25, 1.35)] + [(-1.0, 1.0)] * (degree // 2),
                constraints=[{
                    "type": "ineq",
                    "fun": lambda pars: 1 - pars[1:] @ pars[1:],
                    "jac": lambda pars: np.r_[0.0, -2 * pars[1:]],
                }],
                options={"maxiter": args.maxiter, "ftol": 2e-12, "disp": False},
            )
            val, grad, info = objective.evaluate(result.x, gradient=True)
            replay = ScalarObjective(degree, max(384, 2 * args.nodes)).evaluate(result.x)[0]
            info.update({
                "start_index": index,
                "success": bool(result.success),
                "message": str(result.message),
                "iterations": int(result.nit),
                "high_node_value": replay,
                "quadrature_difference": replay - val,
                "alpha_derivative": float(grad[0]),
                "radial_gradient": float(grad[1:] @ result.x[1:]),
                "parameters": result.x.tolist(),
            })
            records.append(info)
            print(json.dumps({k: info[k] for k in ["degree", "start_index", "value", "derivative_energy", "alpha", "success", "iterations", "quadrature_difference"]}), flush=True)
            if best is None or val > best[0]:
                best = (val, result.x.copy())
        previous = best[1]
    check_objective = ScalarObjective(20, args.nodes)
    check_point = check_objective.resolvent_start()
    _, check_gradient, _ = check_objective.evaluate(check_point, gradient=True)
    errors = []
    for _ in range(6):
        direction = rng.normal(size=len(check_point))
        direction /= np.linalg.norm(direction)
        step = 1e-6
        difference = (check_objective.evaluate(check_point + step * direction)[0]
                      - check_objective.evaluate(check_point - step * direction)[0]) / (2 * step)
        errors.append(float(abs(difference - check_gradient @ direction)))
    output = {
        "status": "floating_diagnostic_not_certificate",
        "gradient_directional_check_max_error": max(errors),
        "runs": records,
    }
    if args.output:
        args.output.write_text(json.dumps(output, indent=2) + "\n")
    else:
        print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
