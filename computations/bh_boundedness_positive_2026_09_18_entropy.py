"""Diagnostic search for unit-circle Fourier entropy / rotation speed.

This is not a certificate and writes no files.  The phase vector is
unrestricted, so no low-degree claim is made for the optimized targets.
"""

import argparse
import json
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize


def fwht(vector):
    out = np.array(vector, dtype=complex, copy=True)
    half = 1
    while half < len(out):
        for start in range(0, len(out), 2 * half):
            left = out[start:start + half].copy()
            right = out[start + half:start + 2 * half].copy()
            out[start:start + half] = left + right
            out[start + half:start + 2 * half] = left - right
        half *= 2
    return out


def adjacency(n):
    size = 2 ** n
    out = np.zeros((size, size))
    for word in range(size):
        for bit in range(n):
            out[word, word ^ (1 << bit)] = .5
    return out


def objective(theta, jx):
    size = len(theta)
    phase = np.exp(1j * theta)
    coefficients = fwht(phase) / size
    probabilities = abs(coefficients) ** 2
    logs = np.log(np.maximum(probabilities, 1e-300))
    entropy = -float(np.dot(probabilities, logs))
    entropy_gradient = 2 * np.imag(phase * np.conjugate(fwht((logs + 1) * coefficients))) / size
    rotated = np.conjugate(phase)[:, None] * jx * phase[None, :]
    eigenvalues, eigenvectors = eigh(rotated - jx, subset_by_index=(size - 1, size - 1))
    speed = float(eigenvalues[0])
    vector = eigenvectors[:, 0]
    speed_gradient = 2 * np.imag(np.conjugate(vector) * (rotated @ vector))
    if speed < 1e-14 or entropy < 1e-14:
        return 1e10, np.zeros(size)
    loss = -entropy / speed
    gradient = -entropy_gradient / speed + entropy * speed_gradient / speed ** 2
    return loss, gradient


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260918)
    args = parser.parse_args()
    jx = adjacency(args.n)
    rng = np.random.default_rng(args.seed)
    test = rng.normal(size=2 ** args.n)
    loss, gradient = objective(test, jx)
    direction = rng.normal(size=len(test)); direction /= np.linalg.norm(direction)
    step = 1e-6
    finite_difference = (objective(test + step * direction, jx)[0]
                         - objective(test - step * direction, jx)[0]) / (2 * step)
    assert abs(finite_difference - gradient @ direction) < 1e-5
    print(json.dumps({"diagnostic_only": True, "n": args.n,
                      "gradient_check_error": float(abs(finite_difference - gradient @ direction))}), flush=True)
    for trial in range(args.starts):
        theta = rng.uniform(-np.pi, np.pi, size=2 ** args.n)
        result = minimize(objective, theta, args=(jx,), jac=True, method="L-BFGS-B",
                          options={"maxiter": args.iterations, "ftol": 1e-11, "gtol": 1e-8})
        print(json.dumps({"trial": trial, "entropy_over_speed": float(-result.fun),
                          "iterations": result.nit, "success": bool(result.success)}), flush=True)


if __name__ == "__main__":
    main()
