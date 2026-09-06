"""Deterministic Gaussian quadrature probe, NOT a rigorous certificate.

Conditional ternary Frank--Wolfe ascent in the finite star bank.
All construction parameters and iteration counts are independent of n.
No files are written; JSON output records the finite trial value.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import ndtr, roots_hermitenorm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stars", type=int, default=6)
    ap.add_argument("--nodes", type=int, default=256)
    ap.add_argument("--steps", type=int, default=160)
    args = ap.parse_args()
    x, weights = roots_hermitenorm(args.nodes)
    weights /= math.sqrt(2 * math.pi)
    hermites = [np.ones_like(x), x.copy()]
    for k in range(2, 2 * args.stars + 1):
        hermites.append((x * hermites[-1] - math.sqrt(k - 1) * hermites[-2]) / math.sqrt(k))
    features = np.stack([hermites[2 * k] for k in range(args.stars + 1)])
    phi = lambda z: np.exp(-np.minimum(z * z, 1500) / 2) / math.sqrt(2 * math.pi)

    def evaluate(q, a):
        p = weights @ q
        t = p - a @ a
        if t <= 0:
            return -math.inf, None
        tau = math.sqrt(t)
        K = a @ features
        z = K / tau
        psi = K * (2 * ndtr(z) - 1) + 2 * tau * phi(z)
        return float(weights @ ((1 - q) * psi)), (p, t, tau, K, z, psi)

    q = np.full_like(x, 0.5)
    a = np.zeros(args.stars + 1)
    a[0] = 0.5 * math.sqrt(2 / math.pi)
    trajectory = []
    for step in range(args.steps):
        value, (p, t, tau, K, z, psi) = evaluate(q, a)
        B = float(weights @ ((1 - q) * phi(z))) / tau
        A = features @ (weights * (1 - q) * (2 * ndtr(z) - 1)) - 2 * B * a
        sigma = float(np.linalg.norm(A[1:]))
        mean = A[0] * x
        threshold_raw = psi - B
        threshold = np.maximum(threshold_raw, 0)
        if sigma > 1e-12:
            zp = (threshold - mean) / sigma
            zm = (-threshold - mean) / sigma
            qnew = ndtr(-zp) + ndtr(zm)
            fmean = ndtr(-zp) - ndtr(zm)
            density = phi(zp) + phi(zm)
            astar = A[1:] / sigma * (weights @ density)
            gain_integrand = mean * fmean + sigma * density - threshold_raw * qnew
        else:
            qnew = (np.abs(mean) > threshold).astype(float)
            fmean = np.sign(mean) * qnew
            astar = np.zeros_like(A[1:])
            gain_integrand = np.maximum(np.abs(mean) - threshold_raw, 0)
        anew = np.concatenate(([weights @ (x * fmean)], astar))
        gap = float(weights @ gain_integrand - A @ a - weights @ ((B - psi) * q))

        def line(theta):
            return evaluate((1 - theta) * q + theta * qnew, (1 - theta) * a + theta * anew)[0]

        grid = np.linspace(0, 1, 33)
        vals = [line(v) for v in grid]
        winner = int(np.argmax(vals))
        lo, hi = grid[max(0, winner - 1)], grid[min(32, winner + 1)]
        res = minimize_scalar(lambda s: -line(s), bounds=(lo, hi), method="bounded", options={"xatol": 1e-12})
        choices = [(float(grid[winner]), vals[winner]), (float(res.x), -float(res.fun)), (0.0, value)]
        theta, newvalue = max(choices, key=lambda pair: pair[1])
        trajectory.append({"step": step, "value": value, "gap": gap, "theta": theta})
        if step % 20 == 0:
            print(json.dumps(trajectory[-1]), flush=True)
        if theta == 0 or newvalue - value < 1e-12:
            break
        q = (1 - theta) * q + theta * qnew
        a = (1 - theta) * a + theta * anew
    value, (p, t, *_rest) = evaluate(q, a)
    print(json.dumps({"status": "quadrature heuristic only", "stars": args.stars, "nodes": args.nodes,
                      "steps": len(trajectory), "value": value, "p": float(p), "tau2": float(t),
                      "a": a.tolist(), "A_last": A.tolist(), "B_last": B,
                      "last": trajectory[-1]}, indent=2))


if __name__ == "__main__":
    main()
